use crate::logger::Logger;
use std::fs::{self};
use std::process;
use std::sync::{Arc, Mutex};

pub struct Fasta {
    pub id:String,
    //pub desc:String,
    pub seq:String,
}

pub fn read_fasta(file : String, logger: &Arc<Mutex<Logger>>) -> Vec<Fasta> {
    println!("read_fasta: processing file: {}", file);

    // read file
    let fasta_file = fs::read_to_string(&file).unwrap_or_else(|error|{
        logger.lock().unwrap().error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    let mut last_id = "";
    //let mut last_desc = "";
    let mut last_sequence = String::from("");
    let mut fasta: Vec<Fasta> = Vec::new();
    for line in fasta_file.lines() {

        // ID and Description
        if line.starts_with(">") {
            if last_id != "" {
                //fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence });
                fasta.push(Fasta { id: last_id.to_string(), seq: last_sequence });
            }
            last_sequence = String::from("");

            match line.find(" ") {
                Some(index) => {
                    last_id = &line[1..index];
                    //last_desc = &line[index+1..];
                },
                None => { 
                    last_id = &line[1..];
                    //last_desc = "";
                }
            };
            //println!("id and desc: {} {}", last_id, last_desc);
        }
        else {
            last_sequence.push_str(line);
        }
    }
    //fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence }); 
    fasta.push(Fasta { id: last_id.to_string(), seq: last_sequence }); 
    logger.lock().unwrap().information("read_fasta: Finished processing file");

    return fasta;
}