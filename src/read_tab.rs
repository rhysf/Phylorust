use crate::logger::Logger;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::fs::File;
use std::collections::{HashMap};
use std::io::{BufRead, BufReader};

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

/// Builds sample_bases from the .tab summary files: sample_name → (contig, pos) → base
pub fn build_sample_bases_from_tabs(
    variant_tab_paths: &[String],
    logger: &Logger,
) -> HashMap<String, HashMap<(String, usize), char>> {

    logger.information("build_sample_bases_from_tabs: rebuilding sample bases from .tab files");
    let mut sample_bases: HashMap<String, HashMap<(String, usize), char>> = HashMap::new();

    for path in variant_tab_paths {
        // Derive sample name from filename (everything after the last '-')
        let filename = Path::new(path)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
        let sample_name = filename
            .replace("reference-", "")
            .replace("snp-", "")
            .replace("deletion-", "")
            .replace("insertion-", "")
            .replace("ambiguous-", "")
            .replace("het_insertion-", "")
            .replace("het_deletion-", "")
            .replace("het-", "");

        logger.information(&format!("build_sample_bases_from_tabs: reading sample '{}'", sample_name));

        let file = File::open(path).unwrap_or_else(|e| {
            logger.error(&format!("Could not open {}: {}", path, e));
            std::process::exit(1);
        });

        let reader = BufReader::new(file);
        let mut contig: Option<String> = None;
        
        // Get (or create) the map for this sample — don't overwrite!
        let sample_map = sample_bases
            .entry(sample_name.clone())
            .or_insert_with(HashMap::new);

        for line in reader.lines() {
            let line = line.unwrap();

            if line.starts_with("##") {
                contig = Some(line.trim_start_matches("##").to_string());
                continue;
            }
            let contig = match &contig {
                Some(c) => c,
                None => continue,
            };

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 3 {
                continue;
            }

            let start: usize = parts[0].parse().unwrap_or(0);
            //let stop: usize = parts[1].parse().unwrap_or(start);
            let bases = parts[2];

            // Write every position within the range
            for (i, base_char) in bases.chars().enumerate() {
                let pos = start + i;
                sample_map.insert((contig.clone(), pos), base_char);
            }
        }
    }

    logger.information(&format!(
        "build_sample_bases_from_tabs: collected {} samples total",
        sample_bases.len()
    ));

    sample_bases
}