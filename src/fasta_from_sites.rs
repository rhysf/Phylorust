use crate::logger::Logger;
use crate::read_fasta::Fasta;
use std::collections::{HashMap};
use std::fs::{File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

pub fn generate_fasta_for_percent_site_set_cached(
    percent: usize,
    histogram_positions: &HashMap<usize, Vec<(String, usize)>>,
    fasta: &Vec<Fasta>,
    sample_bases_cache: &HashMap<String, HashMap<(String, usize), char>>,
    target_dir: &str,
    settings_str: &str,
    logger: &Arc<Mutex<Logger>>) {

    let out_fasta_path = format!("{}/percent_{}-{}.fasta", target_dir, percent, settings_str);

    // Skip if already exists
    if Path::new(&out_fasta_path).exists() {
        logger.lock().unwrap().warning(&format!("generate_fasta_for_percent_site_set_cached: FASTA for {}% already exists — skipping generation.", percent));
        return;
    }

    logger.lock().unwrap().information(&format!("generate_fasta_for_percent_site_set_cached: Generating FASTA for {}% coverage...", percent));

    // Lookup positions for this threshold
    let positions = match histogram_positions.get(&percent) {
        Some(p) if !p.is_empty() => p,
        _ => {
            logger.lock().unwrap().warning(&format!("No positions found for {}% — skipping FASTA.", percent));
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

    // Write directly to file without holding all sequences in memory
    let file = File::create(&out_fasta_path).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Could not create FASTA '{}': {}", out_fasta_path, e));
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    for (sample, base_map) in sample_bases_cache {
        writeln!(writer, ">{}", sample).unwrap();

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

    logger.lock().unwrap().information(&format!("generate_fasta_for_percent_site_set: Saved FASTA to {}", out_fasta_path));
}