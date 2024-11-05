use clap::Parser;
use std::fs;
use std::process;
use std::path::Path;
use std::collections::HashMap;

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
}

fn main() {

    let args = Args::parse();
    let mut name_type_location_hash: HashMap<&str, &str> = HashMap::new();

    // Read Name Type Location file
    println!("Processing file: {}", args.name_type_location_filename);
    let content = fs::read_to_string(&args.name_type_location_filename).unwrap_or_else(|error|{
        println!("Error with file: {} {}", args.name_type_location_filename, error);
        process::exit(1);
    });

    // separate out the columns
    for(index, line) in content.lines().enumerate() {
        let line_parts:Vec<&str> = line.split('\t').collect();

        // check there are 3 columns
        if line_parts.len() != 3 {
            println!("Error with format of file: {} on line number {} = {}", args.name_type_location_filename, index, line);
            process::exit(1);
        }

        // check files exist
        if !Path::new(line_parts[2]).exists() {
            println!("Error: File {} does not exist", line_parts[2]);
            process::exit(1);
        }

        // save into hash map
        name_type_location_hash.insert(line_parts[0], line_parts[2]);
    }
    let number_of_vcfs = name_type_location_hash.len();
    println!("{} VCF files specified", number_of_vcfs);

    println!("Processing file: {}", args.fasta_filename);

    
    
    println!("Hello, ecatools world!");
}