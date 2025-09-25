use clap::Parser;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::collections::{HashSet, HashMap};

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
use read_vcf::VCFEntry;
use logger::Logger;

fn main() {

    let args = Args::parse();
    let logger = Logger;

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

    // Make output folder
    fs::create_dir_all(&args.output_dir).unwrap_or_else(|error|{
        logger.error(&format!("Error with output directory: {}", error));
        process::exit(1);
    });

    // go through each VCF converting to summary files of reference and variant positions
    let mut global_sample_names: HashSet<String> = HashSet::new();
    let mut reference_paths: Vec<String> = Vec::new();
    let mut variant_paths: Vec<String> = Vec::new();
    let mut vcf_entries_by_sample: HashMap<String, Vec<VCFEntry>> = HashMap::new();
    for name_type_location in &name_type_locations {
        logger.information("──────────────────────────────");

        // Output files
        let vcf_path = &name_type_location.location;
        let (outfile_reference_bases, outfile_variant_bases) = generate_output_filenames(&args, &logger, vcf_path);
        reference_paths.push(outfile_reference_bases.clone());
        variant_paths.push(outfile_variant_bases.clone());

        // Read VCF
        let entries = read_vcf::read_vcf(name_type_location, &logger, &mut global_sample_names, &args);
        logger.information(&format!("{}: {} vcf positions parsed", name_type_location.location, entries.len()));
        read_vcf::count_variants(&entries, &logger);

        // Group VCF entries by sample
        for entry in &entries {
            for sample_name in entry.samples_to_base_type.keys().cloned() {
                // Check if this VCF file has a rename entry
                let final_name = if let Some(rename) = rename_map.get(&name_type_location.location) {
                    rename.clone()
                } else {
                    sample_name.clone()
                };

                //logger.information(&format!("Assigning entry to sample: {}", sample_name));
                vcf_entries_by_sample
                    .entry(final_name)
                    .or_default()
                    .push(entry.clone());
            }
        }

        // Skip if already written
        if Path::new(&outfile_reference_bases).exists() && Path::new(&outfile_variant_bases).exists() {
            logger.warning(&format!("Skipping genome/region output for {} — output files already exist.", vcf_path));
            continue;
        }

        // Save genome to hashmap of arrays
        let genome = genome_array::make_hashmap_of_arrays_for_genome(&fasta, &logger);

        // Convert 0->1 for reference, and 0->2 for variants
        let genome = genome_array::fill_genome_hash_array_from_vcf(&logger, &entries, genome, args.settings);

        // Print tab files for locations of 1s (reference)
        genome_array::write_regions_from_genome_array(&genome, 1, &outfile_reference_bases, &logger);

        // Print tab files for locations of 1s (variant)
        genome_array::write_regions_from_genome_array(&genome, 2, &outfile_variant_bases, &logger);
    }
    logger.information("──────────────────────────────");

    // Go through all reference and variant files and find 'ECA' sites
    let variant_counts = read_genome_array_summary::load_contig_position_counts(&variant_paths, &logger);
    let reference_counts = read_genome_array_summary::load_contig_position_counts(&reference_paths, &logger);
    logger.information("──────────────────────────────");

    // Generate histogram
    let histogram_positions = read_genome_array_summary::load_or_generate_histogram(&variant_counts, &reference_counts, &vcf_entries_by_sample, name_type_locations.len(), &args, &logger);
    let settings_str = format!("m-{}-s-{}-e-{}-z-{}", args.min_read_depth, args.settings, args.exclude_contig, args.restrict_contig);
    let histogram_file = format!("{}/site_coverage_histogram-{}.tsv", args.output_dir, settings_str);
    utils::run_r_plotting_script(&histogram_file, &logger, args.percent_threshold, &args.output_dir);

    // Summarise 
    logger.information("──────────────────────────────");
    read_genome_array_summary::write_site_position_files(&histogram_positions, &args.output_dir, &logger);
    logger.information(&format!("vcf_entries_by_sample has {} samples", vcf_entries_by_sample.len()));

    for (sample, entries) in &vcf_entries_by_sample {
        logger.information(&format!("Sample: {} has {} entries", sample, entries.len()));
    }

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

    // generate FASTA(s)
    for &percent in &percentiles_to_generate {
        // First check if histogram_positions has any sites for this percent
        if let Some(sites) = histogram_positions.get(&percent) {
            if sites.is_empty() {
                logger.warning(&format!(
                    "Skipping FASTA for {}% (no sites found in histogram).",
                    percent
                ));
                continue;
            }
        } else {
            logger.warning(&format!(
                "Skipping FASTA for {}% (no entry in histogram_positions).",
                percent
            ));
            continue;
        }

        // old check too (double safety)
        let site_file = format!("{}/percent_{}.tab", args.output_dir, percent);
        if let Ok(metadata) = std::fs::metadata(&site_file) {
            if metadata.len() == 0 {
                logger.warning(&format!(
                    "Skipping FASTA for {}% (tab file is empty: {}).",
                    percent, site_file
                ));
                continue;
            }
        }

        fasta_from_sites::generate_fasta_for_percent_site_set(
            percent,
            &histogram_positions,
            &fasta,
            &vcf_entries_by_sample,
            &args,
            &logger,
        );
    }

    // generate trees
    if !args.skip_fasttree {
        utils::run_fasttree_on_fastas(
            &args.output_dir,
            &percentiles_to_generate,
            &logger,
            args.fasttree_bin.as_deref(),
        );
    } else {
        logger.information("Skipping tree generation (--skip-fasttree set).");
    }
}

fn generate_output_filenames(args: &Args, logger: &Logger, vcf_path: &str) -> (String, String) {

    let vcf_filename = Path::new(vcf_path)
        .file_name()
        .expect("Failed to extract VCF filename")
        .to_string_lossy()
        .to_string();

    let rustatools_settings = format!(
        "m-{}-s-{}-e-{}-z-{}",
        args.min_read_depth,
        args.settings,
        args.exclude_contig,
        args.restrict_contig
    );

    let outfile_reference_bases = format!(
        "{}/{}-{}-reference-bases.tab",
        args.output_dir, vcf_filename, rustatools_settings
    );

    let outfile_variant_bases = format!(
        "{}/{}-{}-variant-bases.tab",
        args.output_dir, vcf_filename, rustatools_settings
    );

    logger.output(&format!("Outfile reference bases = {}", outfile_reference_bases));
    logger.output(&format!("Outfile variant bases = {}", outfile_variant_bases));

    (outfile_reference_bases, outfile_variant_bases)
}