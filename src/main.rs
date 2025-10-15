use clap::Parser;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::collections::{HashSet, HashMap};
use sysinfo::{System};

mod args;
mod fasta_from_sites;
mod genome_array;
mod logger;
mod read_fasta;
mod read_genome_array_summary;
mod read_tab;
mod read_vcf;
mod utils;

use args::Args;
use logger::Logger;

fn main() {

    let args = Args::parse();
    let logger = Logger;

    // Settings
    //let het_code = match args.heterozygosity_encoding.as_str() {
    //    "iupac" => "h-i",
    //    _ => "h-a",
    //};
    let rustatools_summary_string = format!("m-{}-e-{}-z-{}-h-{}",
        args.min_read_depth,
        args.exclude_contig,
        args.restrict_contig,
        args.heterozygosity_encoding
    );

    let rustatools_analysis_string = format!("{}-s-{}", rustatools_summary_string, args.settings);

    //let rustatools_settings_string = format!(
    //    "m-{}-s-{}-e-{}-z-{}-{}",
    //    args.min_read_depth,
    //    args.settings,
    //    args.exclude_contig,
    //    args.restrict_contig,
    //    het_code
    //);

    // Base output folder (phylorust_output by default)
    fs::create_dir_all(&args.output_dir).unwrap_or_else(|error|{
        logger.error(&format!("Error with output directory: {}", error));
        process::exit(1);
    });

    // VCF summaries folder (shared across runs)
    let vcf_summary_dir = format!("{}/vcf_summaries/{}", args.output_dir, rustatools_summary_string);
    fs::create_dir_all(&vcf_summary_dir).unwrap_or_else(|error| {
        logger.error(&format!("Error with VCF summaries directory: {}", error));
        process::exit(1);
    });

    // Run-specific folder (named after the tab file)
    let tabfile_basename = Path::new(&args.name_type_location_filename).to_str().unwrap();
    let run_dir = format!("{}/{}", args.output_dir, tabfile_basename);
    fs::create_dir_all(&run_dir).unwrap_or_else(|error| {
        logger.error(&format!("Error with run-specific directory: {}", error));
        process::exit(1);
    });

    logger.information(&format!("Output directories: vcf_summaries={}, run_dir={}", vcf_summary_dir, run_dir));

    // Read Name Type Location file
    let name_type_locations = read_tab::read_name_type_location_file(&args.name_type_location_filename, &logger);

    // Make VCF path → tab file sample name
    let mut rename_map: HashMap<String, String> = HashMap::new();
    for nt in &name_type_locations {
        // Keyed by the file path from the tab file
        rename_map.insert(nt.location.clone(), nt.name.clone());
    }

    // Read FASTA file to memory
    let fasta = read_fasta::read_fasta(args.fasta_filename.clone(), &logger);

    // go through each VCF converting to summary files of reference and variant positions and bases
    let mut global_sample_names: HashSet<String> = HashSet::new();
    for name_type_location in &name_type_locations {
        logger.information("──────────────────────────────");

        // --- memory snapshot ---
        mem_log("Memory in use before processing VCFs");

        // Output subfolder (phylorust/vcf_summaries/settings/vcf/)
        let vcf_path = &name_type_location.location;
        let vcf_filename = Path::new(vcf_path).file_name().expect("Failed to extract VCF filename").to_string_lossy().to_string();
        let vcf_summary_subdir = format!("{}/{}", vcf_summary_dir, vcf_filename);

        // Skip if already written
        if Path::new(&vcf_summary_subdir).exists() {
            logger.warning(&format!("Skipping genome/region output for {} — output folder already exist.", vcf_summary_subdir));
            continue;
        }

        // Read VCF
        let entries = read_vcf::read_vcf(name_type_location, &logger, &mut global_sample_names, &args);
        logger.information(&format!("{}: {} vcf positions parsed", name_type_location.location, entries.len()));
        read_vcf::count_variants(&entries, &logger);

        // Create VCF summary folder
        fs::create_dir_all(&vcf_summary_subdir).unwrap_or_else(|error| {
            logger.error(&format!("Could not create directory '{}': {}", vcf_summary_subdir, error));
            std::process::exit(1);
        });

        // Convert 0->1 for reference, and 0->2 for snps etc. (all variants are now being written to files here. But refs kept for later)
        let (sample_genomes, base_type_map) = genome_array::fill_genome_hash_array_from_vcf(
            &logger, 
            &entries, 
            &fasta, 
            &vcf_summary_subdir, 
            args.heterozygosity_encoding.as_str());

        // Print tab files for locations of 1s (reference) etc, per-sample summaries:
        for (sample_name, contigs) in &sample_genomes {
            for (base_type, code) in &base_type_map {
                if base_type == "reference" || base_type == "ambiguous" {
                    let outfile = format!("{}/{}-{}.tab", vcf_summary_subdir, base_type, sample_name);
                    genome_array::write_regions_from_genome_array(contigs, *code, &fasta, sample_name, base_type, &outfile, &logger);
                }
            }
        }

        // drop entries
        drop(entries);
        drop(sample_genomes);
        // Encourage allocator to release pages
        #[cfg(target_os = "linux")]
        {
            unsafe { libc::malloc_trim(0); }
        }
        //mem_log("after drop(sample_genomes)");
    }
    logger.information("──────────────────────────────");

    // Dynamically gather all .tab files from subdirectories
    let (variant_paths, reference_paths) = collect_tab_paths_by_settings(&vcf_summary_dir, args.settings, &logger);

    // Load those into memory-efficient counters
    let variant_counts = read_genome_array_summary::load_contig_position_counts(&variant_paths, &logger);
    let reference_counts = read_genome_array_summary::load_contig_position_counts(&reference_paths, &logger);

    logger.information("──────────────────────────────");

    // Generate histogram
    let histogram_positions = read_genome_array_summary::load_or_generate_histogram(
        &variant_counts, 
        &reference_counts, 
        &variant_paths, 
        name_type_locations.len(), 
        &rustatools_analysis_string, 
        &args, 
        &run_dir, 
        &logger);
    let histogram_file = format!("{}/site_coverage_histogram-{}.tsv", run_dir, rustatools_analysis_string);
    utils::run_r_plotting_script(&histogram_file, &logger, args.percent_threshold, &run_dir);

    // Summarise 
    logger.information("──────────────────────────────");
    read_genome_array_summary::write_site_position_files(&histogram_positions, &run_dir, &rustatools_analysis_string, &logger);

    // Determine which percentiles to generate FASTAs for
    let percentiles_to_generate: Vec<usize> = if args.generate_fastas == "all" {
        (1..=100).collect()
    } else if !args.generate_fastas.trim().is_empty() {
        args.generate_fastas
            .split(',')
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .collect()
    } else {
        vec![args.percent_threshold]
    };

    // gather all ref and variant sites:
    let mut all_tab_paths = Vec::new();
    all_tab_paths.extend(reference_paths.clone());
    all_tab_paths.extend(variant_paths.clone());

    // Precompute and cache the sample_bases map across all percents
    logger.information("Prebuilding sample base cache from .tab files (used for all percentiles)...");
    let sample_bases_cache = read_tab::build_sample_bases_from_tabs(&all_tab_paths, &logger);
    logger.information(&format!(
        "Cached sample bases for {} isolates ({} total bases).",
        sample_bases_cache.len(),
        sample_bases_cache.values().map(|m| m.len()).sum::<usize>()
    ));

    // generate FASTA(s)
    for &percent in &percentiles_to_generate {
        // First check if histogram_positions has any sites for this percent
        if let Some(sites) = histogram_positions.get(&percent) {
            if sites.is_empty() {
                logger.warning(&format!("Skipping FASTA for {}% (no sites found in histogram).", percent));
                continue;
            }
        } else {
            logger.warning(&format!("Skipping FASTA for {}% (no entry in histogram_positions).", percent));
            continue;
        }

        fasta_from_sites::generate_fasta_for_percent_site_set_cached(
            percent,
            &histogram_positions,
            &fasta,
            &sample_bases_cache,
            &run_dir,
            &rustatools_analysis_string,
            &logger,
        );
    }

    // generate trees
    if !args.skip_fasttree {
        utils::run_fasttree_on_fastas(
            &run_dir,
            &percentiles_to_generate,
            &rustatools_analysis_string,
            &logger,
            args.fasttree_bin.as_deref(),
        );
    } else {
        logger.information("Skipping tree generation (--skip-fasttree set).");
    }
}

fn mem_log(label: &str) {
    let mut sys = System::new_all();
    sys.refresh_processes();
    if let Ok(pid) = sysinfo::get_current_pid() {
        if let Some(proc_) = sys.process(pid) {
            println!(
                "{}: {:.2} GB RSS",
                label,
                proc_.memory() as f64 / 1_073_741_824.0
            );
        } else {
            println!("{}: process info not found", label);
        }
    } else {
        println!("{}: unable to get current PID", label);
    }
}

pub fn collect_tab_paths_by_settings(
    output_dir: &str,
    settings: u8,
    logger: &Logger,
) -> (Vec<String>, Vec<String>) {
    let mut variant_paths = Vec::new();
    let mut reference_paths = Vec::new();

    logger.information(&format!("Collecting .tab files from {}", output_dir));

    for entry in fs::read_dir(output_dir).expect("Failed to read output_dir") {
        let path = entry.unwrap().path();
        if path.is_dir() {
            // scan subdirectory
            for file in fs::read_dir(&path).expect("Failed to read subdir") {
                let file_path = file.unwrap().path();
                if !file_path.is_file() {
                    continue;
                }
                let fname = file_path.file_name().unwrap().to_string_lossy();

                if fname.starts_with("reference-") {
                    reference_paths.push(file_path.to_string_lossy().to_string());
                } else {
                    let include_variant = match settings {
                        1 => fname.starts_with("snp-"),
                        2 => fname.starts_with("snp-") || fname.starts_with("heterozygous-"),
                        3 => !fname.starts_with("reference-") && !fname.starts_with("ambiguous-"),
                        _ => false,
                    };
                    if include_variant {
                        variant_paths.push(file_path.to_string_lossy().to_string());
                    }
                }
            }
        }
    }

    logger.information(&format!(
        "Collected {} reference and {} variant tab files",
        reference_paths.len(),
        variant_paths.len()
    ));

    (variant_paths, reference_paths)
}