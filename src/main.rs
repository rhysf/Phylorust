use clap::Parser;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::collections::{HashSet, HashMap};
mod logger;
use logger::Logger;

// use std::collections::HashSet;
// use std::ops::Index;
// use std::collections::HashMap;

mod args;
mod read_vcf;
mod read_fasta;
mod read_tab;
mod genome_array;
mod read_genome_array_summary;
mod fasta_from_sites;

use args::Args;
// use read_genome_array_summary::write_site_position_files;
// use read_fasta::Fasta;
// use read_vcf::VCFsamples;
use read_vcf::VCFEntry;
// use read_tab::NameTypeLocation;

fn main() {

    let args = Args::parse();
    let logger = Logger;

    // Read Name Type Location file
    let name_type_locations = read_tab::read_name_type_location_file(&args.name_type_location_filename, &logger);

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
        let entries = read_vcf::read_vcf(name_type_location, &logger, &mut global_sample_names);
        logger.information(&format!("{}: {} vcf positions parsed", name_type_location.location, entries.len()));
        read_vcf::count_variants(&entries, &logger);

        // Group VCF entries by sample
        for entry in &entries {
            for sample_name in entry.samples_to_base_type.keys().cloned() {
                //logger.information(&format!("Assigning entry to sample: {}", sample_name));
                vcf_entries_by_sample
                    .entry(sample_name)
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

    // Summarise
    let histogram_positions = read_genome_array_summary::load_or_generate_histogram(&variant_counts, &reference_counts, &vcf_entries_by_sample, name_type_locations.len(), &args, &logger);
    read_genome_array_summary::write_site_position_files(&histogram_positions, &args.output_dir, &logger);
    // let fasta_for_percent = 
    logger.information(&format!(
        "vcf_entries_by_sample has {} samples",
        vcf_entries_by_sample.len()
    ));

    for (sample, entries) in &vcf_entries_by_sample {
        logger.information(&format!("Sample: {} has {} entries", sample, entries.len()));
    }
    fasta_from_sites::generate_fasta_for_percent_site_set(args.percent_for_tree, &histogram_positions, &fasta, &vcf_entries_by_sample, &args, &logger);
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