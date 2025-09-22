use clap::Parser;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::collections::HashMap;
mod logger;
use logger::Logger;

// use std::collections::HashSet;
// use std::ops::Index;
// use std::collections::HashMap;

mod read_vcf;
mod read_fasta;
mod read_tab;
mod genome_array;

// use read_fasta::Fasta;
// use read_vcf::VCFsamples;
use read_vcf::VCFEntry;
// use read_tab::NameTypeLocation;

// setting up the command line parameters
#[derive(Parser)]
#[command(name = "RustaTools")]
#[command(about = "A simple command-line tool fopr calculating Entirely Covered in All sites", long_about = None)]

struct Args {
    /// reference FASTA
    #[arg(short='f', long="fasta")]
    fasta_filename: String, 

    /// tab delimited file containing isolate name, type (VCF) and VCF location
    #[arg(short='n', long="name_type_location")]
    name_type_location_filename: String,

    /// Output directory 
    #[arg(short='o', long="output_dir", default_value="Rustatools_output")]
    output_dir: String,

    /// Settings 
    #[arg(short='s', long="settings", default_value_t=1)]
    settings: u8,

    // Minimum read depth
    #[arg(short='m', long="min_read_depth", default_value_t=4)]
    min_read_depth: u8,

    /// Exclude variants on contig 
    #[arg(short='e', long="exclude_contig", default_value="n")]
    exclude_contig: String,

    /// Exclude variants not on contig 
    #[arg(short='z', long="restrict_contig", default_value="n")]
    restrict_contig: String,
}

fn main() {

    let args = Args::parse();
    let logger = Logger;

    // Read Name Type Location file
    let name_type_locations = read_tab::read_name_type_location_file(&args.name_type_location_filename, &logger);

    // Info about VCF's
    let number_of_vcfs = name_type_locations.len();
    println!("{} VCF files specified", number_of_vcfs);

    // Read FASTA file to memory
    let fasta = read_fasta::read_fasta(args.fasta_filename, &logger);

    // go through entries
    for entry in &fasta {
        logger.information(&format!("id and desc: {} {}", entry.id, entry.desc));
        logger.information(&format!("length of seq: {}", entry.seq.len()));
    }

    // Make output folder
    fs::create_dir_all(&args.output_dir).unwrap_or_else(|error|{
        logger.error(&format!("Error with output directory: {}", error));
        process::exit(1);
    });

    // Output files
    let name_type_location_filename_path = Path::new(&args.name_type_location_filename);
    let name_type_location_filename_wo_dir = name_type_location_filename_path.file_name().unwrap().to_str();
    logger.output(name_type_location_filename_wo_dir.unwrap());

    let rustatools_settings = &format!("m-{}-s-{}-e-{}-z-{}", 
        args.min_read_depth, 
        args.settings, 
        args.exclude_contig, 
        args.restrict_contig);
    
    let outfile_reference_bases = &format!("{}/{}-{}-reference-bases.tab", 
        args.output_dir,
        name_type_location_filename_wo_dir.unwrap(),
        rustatools_settings);

    let outfile_variant_bases = &format!("{}/{}-{}-variant-bases.tab", 
        args.output_dir,
        name_type_location_filename_wo_dir.unwrap(),
        rustatools_settings);

    logger.output(&format!("Outfile reference bases = {}", outfile_reference_bases));
    logger.output(&format!("Outfile variant bases = {}", outfile_variant_bases));

    // Save genome to hashmap of arrays
    let mut genome = genome_array::make_hashmap_of_arrays_for_genome(&fasta, &logger);

    // go through each VCF
    let mut all_entries: HashMap<String, Vec<VCFEntry>> = HashMap::new();
    for name_type_location in &name_type_locations {
        let entries = read_vcf::read_vcf(name_type_location, &logger);

        logger.information(&format!(
            "{}: {} variants parsed",
            name_type_location.location,
            entries.len()
        ));

        if let Some(first_entry) = entries.first() {
            let sample_names: Vec<String> = first_entry
                .samples
                .values()
                .cloned()
                .collect();
    
            logger.information(&format!("Samples: {:?}", sample_names));
        }

        // Variants per sample
        let mut sample_variant_counts: HashMap<String, usize> = HashMap::new();

        for entry in &entries {
            for (sample_name, genotype) in &entry.samples_to_genotype {
                if !genotype.is_empty() && genotype != "." {
                    *sample_variant_counts.entry(sample_name.clone()).or_insert(0) += 1;
                }
            }
        }

        for (sample, count) in &sample_variant_counts {
            logger.information(&format!("Sample '{}' has {} variants", sample, count));
        }

        // Show a few example variants
        logger.information("First variant:");
        for entry in entries.iter().take(1) {
            logger.information(&format!(
                "{}:{} {}>{} (GT={:?})",
                entry.contig,
                entry.position,
                entry.ref_base,
                entry.alt_base,
                entry.samples_to_genotype
            ));
    }

    logger.information("──────────────────────────────");
    
        // store all entries across files
        all_entries.insert(name_type_location.location.clone(), entries);
    }

    //# Convert ref bases to 1 and variants to 2
    // Start by doing this for a single VCF, and then put in a loop to do it for every VCF
    //$genome_array = genomearray::fill_genome_hash_array_from_vcf($genome_array, $opt_v, $opt_s, $opt_e, $opt_z, $opt_m);

    //# Print tab files for locations of 1s (reference)
    //genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 1, $outfile1);

    //# Print tab files for locations of 1s (variant)
    //genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 2, $outfile2);

    logger.output("output test");
    logger.warning("output test");
}