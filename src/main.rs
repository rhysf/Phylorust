use clap::Parser;
use std::fs;
use std::process;
use std::path::Path;
use std::collections::HashMap;
mod logger;
use logger::Logger;

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

struct Fasta {
    id:String,
    desc:String,
    seq:String,
}

fn read_name_type_location_file(file : String, logger : &Logger) -> HashMap<String, String> {
    logger.information(&format!("Processing file: {}", file));

    let mut name_type_location_hash: HashMap<String, String> = HashMap::new();

    let content = fs::read_to_string(&file).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    for(index, line) in content.lines().enumerate() {
        let line_parts:Vec<&str> = line.split('\t').collect();

        // check there are 3 columns
        if line_parts.len() != 3 {
            logger.error(&format!("Error with format of file: {} on line number {} = {}", file, index, line));
            process::exit(1);
        }

        // check files exist
        if !Path::new(line_parts[2]).exists() {
            logger.error(&format!("Error: File {} does not exist", line_parts[2]));
            process::exit(1);
        }

        // save into hash map
        name_type_location_hash.insert(line_parts[0].to_string(), line_parts[2].to_string());
    }
    return name_type_location_hash;
}

fn read_fasta(file : String, logger : &Logger) -> Vec<Fasta> {
    println!("Processing file: {}", file);

    // read file
    let fasta_file = fs::read_to_string(&file).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    let mut last_id = "";
    let mut last_desc = "";
    let mut last_sequence = String::from("");
    let mut fasta: Vec<Fasta> = Vec::new();
    for line in fasta_file.lines() {

        // ID and Description
        if line.starts_with(">") {
            if last_id != "" {
                fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence });
            }
            last_sequence = String::from("");

            match line.find(" ") {
                Some(index) => {
                    last_id = &line[1..index];
                    last_desc = &line[index+1..];
                },
                None => { 
                    last_id = &line[1..];
                    last_desc = "";
                }
            };
            //println!("id and desc: {} {}", last_id, last_desc);
        }
        else {
            last_sequence.push_str(line);
        }
    }
    fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence }); 
    logger.information("Finished processing file");

    return fasta;
}

fn main() {

    let args = Args::parse();
    let logger = Logger;

    // Read Name Type Location file
    let name_type_location_hash = read_name_type_location_file(args.name_type_location_filename, &logger);

    // Info about VCF's
    let number_of_vcfs = name_type_location_hash.len();
    println!("{} VCF files specified", number_of_vcfs);

    // Read FASTA file to memory
    let fasta = read_fasta(args.fasta_filename, &logger);

    // go through entries
    for entry in fasta {
        logger.information(&format!("id and desc: {} {}", entry.id, entry.desc));
        logger.information(&format!("length of seq: {}", entry.seq.len()));
    }

    logger.output("output test");
    logger.warning("output test");
}