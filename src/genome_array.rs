use crate::logger::Logger;
use crate::read_fasta::Fasta;
use crate::read_vcf::VCFEntry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn make_hashmap_of_arrays_for_genome(fasta : &Vec<Fasta>, logger : &Logger) -> HashMap<String, Vec<i32>> {
    logger.information("make_hashmap_of_arrays_for_genome: initialising genome array...");

    let mut genome = HashMap::new();

    // go through fasta and make contig array of all zeros
    for entry in fasta {
        let contig_array = vec![0;entry.seq.len()];
        genome.insert(entry.id.to_string(), contig_array);
    }

    return genome;
}

pub fn fill_genome_hash_array_from_vcf(logger: &Logger, entries: &Vec<VCFEntry>, mut genome: HashMap<String, Vec<i32>>, settings: u8) -> HashMap<String, Vec<i32>> {
    logger.information("make_hashmap_of_arrays_for_genome: filling genome array with vcf entries (1=ref, 2=variant)...");

    let mut refs_in_array = 0;
    let mut snps_in_array = 0;
    let mut indels_in_array = 0;

    for entry in entries {
        let contig = &entry.contig;
        let position = entry.position;

        // wrongly formatted position in vcf
        if position == 0 {
            logger.error(&format!("fill_genome_hash_array_from_vcf: Encountered VCF position 0 (invalid in 1-based VCF format) on contig {}", contig));
            std::process::exit(1);
        }

        // convert vcf (1-based coords, to the rust 0-based array
        let pos = &entry.position - 1; 

        // Check if contig exists in genome
        let contig_array = match genome.get_mut(contig) {
            Some(array) => array,
            None => {
                logger.warning(&format!("fill_genome_hash_array_from_vcf: Contig '{}' not found in genome array — skipping entry at VCF position {}.", contig, position));
                continue;
            }
        };

        // Check if position is within the length of the contig
        if pos >= contig_array.len() {
            logger.error(&format!("fill_genome_hash_array_from_vcf: Position {} (0-based: {}) out of bounds for contig '{}' (length = {}). Skipping.", position, pos, contig, contig_array.len()));
            std::process::exit(1);
        }

        for (_sample, base_type) in &entry.samples_to_base_type {
            let should_write = match settings {
                1 => base_type == "snp" || base_type == "reference",
                2 => base_type == "snp" || base_type == "reference" || base_type == "heterozygous",
                3 => base_type != "ambiguous", // accept everything else
                _ => {
                    logger.error(&format!("Invalid settings value: {}", settings));
                    std::process::exit(1);
                }
            };

            if !should_write {
                continue;
            }

            let val = if base_type == "reference" {
                1
            } else {
                2
            };

            // update the genome array
            if let Some(array) = genome.get_mut(contig) {
                if pos < array.len() {
                    array[pos] = val;

                    match base_type.as_str() {
                        "reference" => refs_in_array += 1,
                        "snp" | "heterozygous" => snps_in_array += 1,
                        "insertion" | "het_insertion" | "deletion" | "het_deletion" => indels_in_array += 1,
                        _ => {}, // do nothing for ambiguous or others
                    }
                } else {
                    logger.warning(&format!("Position {} out of bounds for contig '{}' (length = {})", pos, contig, array.len()));
                }
            }
        }
    }
    logger.information(&format!("make_hashmap_of_arrays_for_genome: array changes = reference bases: {}, SNPs: {}, indels: {}", refs_in_array, snps_in_array, indels_in_array));
    genome
}

pub fn write_regions_from_genome_array(genome: &HashMap<String, Vec<i32>>, find_value: i32, outfile_path: &str, logger: &Logger) {
    logger.information(&format!("write_regions_from_genome_array: finding value {}...", find_value));

    let file = File::create(outfile_path).unwrap_or_else(|e| {
        logger.error(&format!("Could not create file '{}': {}", outfile_path, e));
        std::process::exit(1);
    });

    let mut writer = BufWriter::new(file);
    let mut count_found = 0;

    for (contig, values) in genome {
        let mut start: Option<usize> = None;
        let mut printed_header = false;
        let length = values.len();

        for i in 0..=length {
            let current = if i < length { values[i] } else { -1 }; // Sentinel: fake value to close run

            if current != find_value {
                if let Some(start_pos) = start {
                    let end = i;
                    count_found += end - start_pos;

                    if !printed_header {
                        writeln!(writer, "##{}", contig).expect("Failed to write header");
                        printed_header = true;
                    }

                    writeln!(writer, "{}\t{}", start_pos, end).expect("Failed to write region");
                    start = None;
                }
            } else {
                if start.is_none() {
                    start = Some(i);
                }
            }
        }
    }
    logger.information(&format!("write_regions_from_genome_array: found {} positions with value {}", count_found, find_value));
}