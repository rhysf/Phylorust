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

    logger.information(&format!("Generating FASTA for {}% coverage...", percent));

    // Lookup positions for this threshold
    let positions = match histogram_positions.get(&percent) {
        Some(p) if !p.is_empty() => p,
        _ => {
            logger.warning(&format!("No positions found for {}% — skipping FASTA.", percent));
            return;
        }
    };

    // Build a minimal reference map only for required contig/pos
    let mut ref_map: HashMap<(String, usize), char> = HashMap::new();
    for (contig, pos) in positions {
        if let Some(fasta_entry) = fasta.iter().find(|f| f.id == *contig) {
            if let Some(base) = fasta_entry.seq.chars().nth(*pos - 1) {
                ref_map.insert((contig.clone(), *pos), base);
            }
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

    if sample_bases.is_empty() {
        logger.warning("No sample bases were collected. Skipping FASTA writing.");
        return;
    }

    // Write directly to file without holding all sequences in memory
    let file = File::create(&out_fasta_path).unwrap_or_else(|e| {
        logger.error(&format!("Could not create FASTA '{}': {}", out_fasta_path, e));
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    // Stream each sample’s sequence
    for sample in sample_bases.keys() {
        writeln!(writer, ">{}", sample).unwrap();
        let base_map = &sample_bases[sample];

        for (contig, pos) in positions {
            if let Some(&base) = base_map.get(&(contig.clone(), *pos)) {
                write!(writer, "{}", base).unwrap();
            } else if let Some(&ref_base) = ref_map.get(&(contig.clone(), *pos)) {
                write!(writer, "{}", ref_base).unwrap();
            } else {
                write!(writer, "N").unwrap();
            }
        }
        writeln!(writer).unwrap();
    }

    logger.information(&format!("Saved FASTA to {}", out_fasta_path));
}