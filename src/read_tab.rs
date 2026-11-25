use crate::logger::Logger;
use std::fs::{self};
use std::process;
use std::path::Path;
use std::fs::File;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};

pub struct NameTypeLocation {
    pub name:String,
    //pub filetype:String,
    pub location:String,
}

pub fn read_name_type_location_file(file : &str, logger: &Arc<Mutex<Logger>>) -> Vec<NameTypeLocation> {
    logger.lock().unwrap().information(&format!("read_name_type_location_file: Processing file: {}", file));

    let mut name_type_locations = Vec::new();

    let content = fs::read_to_string(file).unwrap_or_else(|error|{
        logger.lock().unwrap().error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    for(index, line) in content.lines().enumerate() {
        let line_parts:Vec<&str> = line.split('\t').collect();

        // check there are 3 columns
        if line_parts.len() != 3 {
            logger.lock().unwrap().error(&format!("Error with format of file: {} on line number {} = {}", file, index, line));
            process::exit(1);
        }

        // check files exist
        if !Path::new(line_parts[2]).exists() {
            logger.lock().unwrap().error(&format!("Error: File {} does not exist", line_parts[2]));
            process::exit(1);
        }

        // check that 2nd column says vcf (case-insensitive)
        let filetype = line_parts[1].trim();
        if filetype.to_lowercase() != "vcf" {
            logger.lock().unwrap().warning(&format!("Ignoring line {}: '{}'. Unexpected filetype '{}', expected 'VCF'", index + 1, line, filetype));
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

/// One encoded call per sample, contig, and position.
pub type SampleCalls = BTreeMap<String, HashMap<String, HashMap<usize, u8>>>;

pub fn build_sample_calls(
    summaries: &[(Vec<String>, Vec<String>)],
    settings: u8,
    logger: &Arc<Mutex<Logger>>,
) -> SampleCalls {
    let mut calls = SampleCalls::new();
    let mut candidates: HashMap<String, BTreeSet<usize>> = HashMap::new();
    let mut reference_paths = Vec::new();

    for (variants, references) in summaries {
        let mut names = HashMap::new();
        for path in references.iter().chain(variants) {
            let filename = Path::new(path).file_stem().unwrap().to_string_lossy();
            let (_, raw_name) = filename.split_once('-').unwrap();
            names.entry(raw_name.to_owned()).or_insert_with(|| {
                let mut name = raw_name.to_owned();
                let mut suffix = 1;
                while calls.contains_key(&name) {
                    name = format!("{raw_name}-{suffix}");
                    suffix += 1;
                }
                calls.insert(name.clone(), HashMap::new());
                name
            });
        }

        for path in variants {
            let filename = Path::new(path).file_stem().unwrap().to_string_lossy();
            let (_, raw_name) = filename.split_once('-').unwrap();
            let sample = calls.get_mut(&names[raw_name]).unwrap();
            read_summary_rows(path, logger, |contig, start, end, bases| {
                if settings < 3 && (start != end || bases.len() != 1) {
                    return;
                }
                let positions = candidates.entry(contig.to_owned()).or_default();
                let sample_contig = sample.entry(contig.to_owned()).or_default();
                // Summary rows use reference coordinates. Inserted bases have no
                // corresponding position beyond the end of the row.
                for (position, base) in (start..=end).zip(bases.bytes()) {
                    if let Some(base) = normalize_call(base) {
                        positions.insert(position);
                        sample_contig.insert(position, base);
                    }
                }
            });
        }

        for path in references {
            let filename = Path::new(path).file_stem().unwrap().to_string_lossy();
            let (_, raw_name) = filename.split_once('-').unwrap();
            reference_paths.push((path, names[raw_name].clone()));
        }
    }

    for (path, name) in reference_paths {
        let sample = calls.get_mut(&name).unwrap();
        read_summary_rows(path, logger, |contig, start, end, bases| {
            let Some(positions) = candidates.get(contig) else { return };
            let sample_contig = sample.entry(contig.to_owned()).or_default();
            for &position in positions.range(start..=end) {
                if let Some(&base) = bases.as_bytes().get(position - start) {
                    let base = base.to_ascii_uppercase();
                    if b"ACGT".contains(&base) {
                        sample_contig.entry(position).or_insert(base);
                    }
                }
            }
        });
    }

    calls
}

fn normalize_call(base: u8) -> Option<u8> {
    let base = base.to_ascii_uppercase();
    // Retain IUPAC heterozygous encodings and deletion gaps.
    b"ACGTRYSWKMBDHV-".contains(&base).then_some(base)
}

fn read_summary_rows(
    path: &str,
    logger: &Arc<Mutex<Logger>>,
    mut visit: impl FnMut(&str, usize, usize, &str),
) {
    let file = File::open(path).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Could not open {path}: {e}"));
        process::exit(1);
    });
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut contig = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line).unwrap_or_else(|e| {
            logger.lock().unwrap().error(&format!("Could not read {path}: {e}"));
            process::exit(1);
        }) == 0 {
            break;
        }
        let line = line.trim_end_matches(['\r', '\n']);
        if let Some(name) = line.strip_prefix("##") {
            contig.clear();
            contig.push_str(name);
            continue;
        }

        let mut fields = line.splitn(3, '\t');
        let start = fields.next().and_then(|value| value.parse::<usize>().ok());
        let end = fields.next().and_then(|value| value.parse::<usize>().ok());
        if let (Some(start), Some(end), Some(bases)) = (start, end, fields.next()) {
            if !contig.is_empty() && start > 0 && start <= end {
                visit(&contig, start, end, bases);
            }
        }
    }
}
