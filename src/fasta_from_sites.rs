use std::collections::{HashMap};
use std::fs::{File};
use std::io::{BufWriter, Write};
use std::path::Path;
use crate::args::Args;
use crate::logger::Logger;
use crate::read_vcf::VCFEntry;
use crate::read_fasta::Fasta;

pub fn generate_fasta_for_percent_site_set(
    percent: usize,
    histogram_positions: &HashMap<usize, Vec<(String, usize)>>,
    fasta: &Vec<Fasta>,
    vcf_entries_by_sample: &HashMap<String, Vec<VCFEntry>>,
    args: &Args,
    logger: &Logger,
) {
    let out_fasta_path = format!("{}/percent_{}.fasta", args.output_dir, percent);
    println!("generate_fasta_for_percent_site_set: output: {}", out_fasta_path);

    // Skip if already exists
    if Path::new(&out_fasta_path).exists() {
        logger.warning(&format!("FASTA for {}% already exists — skipping generation.", percent));
        return;
    }

    // Print/log a few values from histogram_positions.get(&90) right before processing
    //if let Some(sites) = histogram_positions.get(&args.percent_for_tree) {
    //    logger.information(&format!("Found {} sites for {}% threshold", sites.len(), args.percent_for_tree));
    //    for (contig, pos) in sites.iter().take(5) {
    //        logger.information(&format!("Example site: {}:{}", contig, pos));
    //    }
    //}

    logger.information(&format!("Generating FASTA for {}% coverage...", percent));

    // Build reference map
    let mut ref_map: HashMap<(String, usize), char> = HashMap::new();
    for fasta_entry in fasta {
        let chars: Vec<char> = fasta_entry.seq.chars().collect();
        for (i, base) in chars.iter().enumerate() {
            ref_map.insert((fasta_entry.id.clone(), i + 1), *base);
        }
    }

    // Build a lookup: sample → (contig,pos) → base
    let mut sample_bases: HashMap<String, HashMap<(String, usize), char>> = HashMap::new();
    for (sample, entries) in vcf_entries_by_sample {
        //logger.information(&format!("Processing sample: {}", sample));

        let mut pos_map = HashMap::new();
        for entry in entries {
            if let Some(base) = entry.samples_to_base1.get(sample) {
                let base = base.chars().next().unwrap_or('N');
                pos_map.insert((entry.contig.clone(), entry.position), base);
            }
            else {
                logger.warning(&format!(
                    "Sample '{}' not found in entry.samples_to_base1 for contig {}, pos {}",
                    sample, entry.contig, entry.position
                ));
            }
        }
        sample_bases.insert(sample.clone(), pos_map);
    }

    // Prepare output: sample_name → sequence string
    let mut sample_to_sequence: HashMap<String, String> = HashMap::new();

    let positions = match histogram_positions.get(&percent) {
        Some(p) => p,
        None => {
            logger.warning(&format!("No positions found for {}% — skipping FASTA.", percent));
            return;
        }
    };

    // Get set of all samples
    let all_samples: Vec<String> = sample_bases.keys().cloned().collect();

    // For each sample, build sequence
    for sample in &all_samples {
        let mut seq = String::with_capacity(positions.len());

        let base_map = match sample_bases.get(sample) {
            Some(map) => map,
            None => {
                logger.warning(&format!("Sample '{}' not found in sample_bases", sample));
                continue;
            }
        };

        for (contig, pos) in positions {
            if let Some(&base) = base_map.get(&(contig.clone(), *pos)) {
                seq.push(base);
            } else if let Some(&ref_base) = ref_map.get(&(contig.clone(), *pos)) {
                seq.push(ref_base);  // fallback to reference
            } else {
                seq.push('N');  // unknown
            }
        }

        sample_to_sequence.insert(sample.clone(), seq);
    }

    if sample_bases.is_empty() {
        logger.warning("No sample bases were collected. Skipping FASTA writing.");
        return;
    }

    // Write FASTA file
    let file = File::create(&out_fasta_path).unwrap_or_else(|e| {
        logger.error(&format!("Could not create FASTA '{}': {}", out_fasta_path, e));
        std::process::exit(1);
    });

    let mut writer = BufWriter::new(file);

    for (sample, seq) in sample_to_sequence {
        writeln!(writer, ">{}", sample).unwrap();
        writeln!(writer, "{}", seq).unwrap();
    }

    logger.information(&format!("Saved FASTA to {}", out_fasta_path));
}