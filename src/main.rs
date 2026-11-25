use clap::Parser;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::collections::{HashMap, HashSet};
//use sysinfo::{System};

mod args;
mod fasta_from_sites;
mod genome_array;
mod logger;
mod read_fasta;
mod read_genome_array_summary;
mod read_tab;
mod read_vcf;
mod utils;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::{Arc, Mutex};

use args::Args;
use logger::Logger;
use crate::logger::LogExt;
use crate::read_tab::NameTypeLocation;

fn main() {

    let args = Args::parse();
    let logger = Arc::new(Mutex::new(Logger));

    logger.lock().unwrap().information(&format!("Phylorust: Using up to {} thread(s)", args.threads));

    // Record settings
    let rustatools_summary_string = format!("calls-v2-m-{}-e-{}-z-{}-h-{}", args.min_read_depth, args.exclude_contig, args.restrict_contig, args.heterozygosity_encoding);
    let invariant_suffix = if args.include_invariants { "-invariants" } else { "" };
    let rustatools_analysis_string = format!("{}-s-{}{}", rustatools_summary_string, args.settings, invariant_suffix);

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
    let mut name_type_locations = read_tab::read_name_type_location_file(&args.name_type_location_filename, &logger);
    let mut input_paths = HashMap::new();
    name_type_locations.retain(|input| {
        let path = Path::new(&input.location);
        let filename = path.file_name().unwrap().to_owned();
        match input_paths.insert(filename, path.to_path_buf()) {
            None => true,
            Some(previous) if previous == path => false,
            Some(previous) => {
                logger.error(&format!("VCFs with the same filename: {} and {}", previous.display(), path.display()));
                process::exit(1);
            }
        }
    });

    // Read FASTA file to memory
    let fasta = read_fasta::read_fasta(args.fasta_filename.clone(), &logger);

    // go through each VCF converting to summary files of reference and variant positions and bases
    //let mut global_sample_names: HashSet<String> = HashSet::new();
    let global_sample_names = Arc::new(Mutex::new(HashSet::new()));
    //for name_type_location in &name_type_locations {
    logger.information("──────────────────────────────");
    logger.information("Reading VCFs...");

    // Create a custom Rayon pool limited by CLI threads arg
    let pool = ThreadPoolBuilder::new().num_threads(args.threads).build().unwrap();

    pool.install(|| {
        name_type_locations.par_iter().for_each(|name_type_location| {

            let logger = logger.clone();
            
            // Output subfolder (phylorust/vcf_summaries/settings/vcf/)
            let vcf_path = &name_type_location.location;
            let vcf_filename = Path::new(vcf_path).file_name().expect("Failed to extract VCF filename").to_string_lossy().to_string();
            let vcf_summary_subdir = format!("{}/{}", vcf_summary_dir, vcf_filename);

            // Skip if already written
            if Path::new(&vcf_summary_subdir).exists() {
                let l = logger.lock().unwrap();
                l.warning(&format!("Skipping genome/region output for {} — output folder already exists.", vcf_summary_subdir));
                return;
            }

            // Publish the cache only after every summary file has been written.
            let pending_dir = format!("{vcf_summary_subdir}.partial-{}", process::id());

            // Read VCF
            let (entries, sample_names) = read_vcf::read_vcf(
                name_type_location, 
                &mut global_sample_names.lock().unwrap(),
                &args,
                &logger);
            {
                let l = logger.lock().unwrap();
                l.information(&format!("{}: {} VCF positions parsed", name_type_location.location, entries.len()));
            }
            read_vcf::count_variants(&entries, &logger.lock().unwrap());

            // Summarize and write
            let sample_genomes = genome_array::summarise_vcf_to_tab_files_and_genome_array(
                &entries,
                &sample_names,
                &fasta,
                &pending_dir,
                &args,
                &logger.lock().unwrap());

            // Write reference/ambiguous regions
            for (sample_name, contigs) in &sample_genomes {
                for base_type in ["reference", "ambiguous"] {
                    let code = genome_array::code_for_base_type(base_type);
                    let outfile = format!("{}/{}-{}.tab", pending_dir, base_type, sample_name);
                    genome_array::write_regions_from_genome_array(
                        contigs, 
                        code,
                        &fasta, 
                        sample_name, 
                        base_type, 
                        &outfile, 
                        &logger
                    );
                }
            }

            if let Err(e) = fs::rename(&pending_dir, &vcf_summary_subdir) {
                if !Path::new(&vcf_summary_subdir).is_dir() {
                    logger.error(&format!("Could not publish {vcf_summary_subdir}: {e}"));
                    process::exit(1);
                }
                // Another process published this summary first.
                if let Err(e) = fs::remove_dir_all(&pending_dir) {
                    logger.warning(&format!("Could not remove {pending_dir}: {e}"));
                }
            }

            drop(entries);
            drop(sample_genomes);

            // Encourage allocator to release pages
            #[cfg(target_os = "linux")]
            unsafe {
                libc::malloc_trim(0);
            }
        });
    });

    logger.information("──────────────────────────────");

    let summaries = collect_sample_summary_paths(
        &vcf_summary_dir, args.settings, &name_type_locations, &logger,
    );
    let sample_calls = read_tab::build_sample_calls(&summaries, args.settings, &logger);
    logger.information(&format!("Loaded calls for {} isolates.", sample_calls.len()));

    let histogram_positions = read_genome_array_summary::load_or_generate_histogram(
        &sample_calls,
        &rustatools_analysis_string,
        &args,
        &run_dir,
        &logger,
    );

    let histogram_file = format!("{}/site_coverage_histogram-{}.tsv", run_dir, rustatools_analysis_string);
    utils::run_r_plotting_script(&histogram_file, args.percent_threshold, &run_dir, &logger);

    // Summarise 
    logger.information("──────────────────────────────");

    read_genome_array_summary::write_site_position_files(&histogram_positions, &run_dir, &rustatools_analysis_string, &logger);

    // Determine which percentiles to generate FASTAs for
    let percentiles_to_generate: Vec<usize> = if args.generate_fastas == "all" {
        (1..=100).collect()
    } else if !args.generate_fastas.trim().is_empty() {
        args.generate_fastas.split(',').filter_map(|s| s.trim().parse::<usize>().ok()).collect()
    } else {
        vec![args.percent_threshold]
    };

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
            &sample_calls,
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

/*fn mem_log(label: &str) {
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
    */

fn collect_sample_summary_paths(
    output_dir: &str,
    settings: u8,
    name_type_locations: &[NameTypeLocation],
    logger: &Arc<Mutex<Logger>>) -> Vec<(Vec<String>, Vec<String>)> {
    let mut summaries = Vec::new();

    for nt in name_type_locations {
        let vcf_filename = Path::new(&nt.location).file_name().unwrap();
        let summary_dir = Path::new(output_dir).join(vcf_filename);
        let entries = fs::read_dir(&summary_dir).unwrap_or_else(|e| {
            logger.error(&format!("Could not read summaries for '{}': {e}", nt.name));
            process::exit(1);
        });
        let mut variants = Vec::new();
        let mut references = Vec::new();
        for entry in entries {
            let path = entry.expect("Could not read summary directory entry").path();
            if path.extension().and_then(|ext| ext.to_str()) != Some("tab") {
                continue;
            }
            let filename = path.file_stem().unwrap().to_string_lossy();
            let Some((kind, _)) = filename.split_once('-') else { continue };
            match kind {
                "reference" => references.push(path.to_string_lossy().into_owned()),
                "snp" => variants.push(path.to_string_lossy().into_owned()),
                "heterozygous" if settings >= 2 => variants.push(path.to_string_lossy().into_owned()),
                "insertion" | "deletion" | "het_insertion" | "het_deletion" if settings == 3 => {
                    variants.push(path.to_string_lossy().into_owned());
                }
                _ => {}
            }
        }
        variants.sort();
        references.sort();
        summaries.push((variants, references));
    }

    summaries
}
