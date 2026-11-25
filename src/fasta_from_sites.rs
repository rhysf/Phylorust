use crate::logger::Logger;
use crate::read_tab::SampleCalls;
use std::collections::{HashMap};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

pub fn generate_fasta_for_percent_site_set_cached(
    percent: usize,
    histogram_positions: &HashMap<usize, Vec<(String, usize)>>,
    sample_calls: &SampleCalls,
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

    // Write directly to file without holding all sequences in memory
    let pending = Path::new(&out_fasta_path).with_extension(format!("partial-{}", std::process::id()));
    let file = File::create(&pending).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Could not create FASTA '{}': {}", out_fasta_path, e));
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    for (sample, calls) in sample_calls {
        writeln!(writer, ">{sample}").unwrap();

        for (contig, position) in positions {
            let base = calls.get(contig).and_then(|sites| sites.get(position)).copied();
            writer.write_all(&[base.unwrap_or(b'N')]).unwrap();
        }
        writeln!(writer).unwrap();
    }
    writer.flush().unwrap();
    fs::rename(pending, &out_fasta_path).expect("Cannot publish FASTA");

    logger.lock().unwrap().information(&format!("generate_fasta_for_percent_site_set: Saved FASTA to {}", out_fasta_path));
}