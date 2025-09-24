use crate::logger::Logger;
use std::fs::{self};
use std::process;
use std::path::Path;

pub struct NameTypeLocation {
    pub name:String,
    //pub filetype:String,
    pub location:String,
}

pub fn read_name_type_location_file(file : &str, logger : &Logger) -> Vec<NameTypeLocation> {
    logger.information(&format!("read_name_type_location_file: Processing file: {}", file));

    let mut name_type_locations = Vec::new();

    let content = fs::read_to_string(file).unwrap_or_else(|error|{
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

        // check that 2nd column says vcf (case-insensitive)
        let filetype = line_parts[1].trim();
        if filetype.to_lowercase() != "vcf" {
            logger.warning(&format!("Ignoring line {}: '{}'. Unexpected filetype '{}', expected 'VCF'", index + 1, line, filetype));
            continue;
        }

        // save into hash map
        name_type_locations.push(NameTypeLocation{
            name:line_parts[0].to_string(),
            //filetype:line_parts[1].to_string(),
            location:line_parts[2].to_string(),
        });
    }
    return name_type_locations;
}