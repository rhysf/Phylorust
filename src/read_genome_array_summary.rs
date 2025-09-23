use crate::logger::Logger;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use crate::args::Args;
use crate::read_vcf::VCFEntry;

pub fn load_contig_position_counts(file_paths: &[String], logger: &Logger) -> HashMap<String, HashMap<usize, usize>> {
    let mut contig_position_counts: HashMap<String, HashMap<usize, usize>> = HashMap::new();

    for path in file_paths {
        logger.information(&format!("load_contig_position_counts: Loading position counts from file: {}", path));

        let file = File::open(path).unwrap_or_else(|e| {
            logger.error(&format!("Could not open {}: {}", path, e));
            std::process::exit(1);
        });

        let reader = BufReader::new(file);
        let mut current_contig: Option<String> = None;
        let mut positions_added_in_file = 0;

        for line in reader.lines() {
            let line = line.expect("Error reading line");

            if line.starts_with("##") {
                current_contig = Some(line.trim_start_matches("##").to_string());
                continue;
            }

            let contig = match &current_contig {
                Some(c) => c,
                None => {
                    logger.error("No contig header found before position line.");
                    continue;
                }
            };

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 2 {
                logger.warning(&format!("Invalid line format: {}", line));
                continue;
            }

            let start: usize = parts[0].parse().unwrap_or(0);
            let stop: usize = parts[1].parse().unwrap_or(start + 1);

            let contig_counts = contig_position_counts.entry(contig.clone()).or_default();
            for pos in start..stop {
                *contig_counts.entry(pos).or_insert(0) += 1;
                positions_added_in_file += 1;
            }
        }
        logger.information(&format!("load_contig_position_counts: Positions added from {}: {}", path, positions_added_in_file));
    }
    contig_position_counts
}

/// Builds sample_bases: sample_name → (contig, pos) → base
pub fn build_sample_bases(vcf_entries_by_sample: &HashMap<String, Vec<VCFEntry>>, logger: &Logger) -> HashMap<String, HashMap<(String, usize), char>> {
    let mut sample_bases: HashMap<String, HashMap<(String, usize), char>> = HashMap::new();

    for (sample, entries) in vcf_entries_by_sample {
        let mut pos_map = HashMap::new();

        for entry in entries {
            if let Some(base) = entry.samples_to_base1.get(sample) {
                let base_char = base.chars().next().unwrap_or('N');
                pos_map.insert((entry.contig.clone(), entry.position), base_char);
            }
        }

        sample_bases.insert(sample.clone(), pos_map);
    }

    sample_bases
}

/// Either loads the histogram + site position data from disk,
/// or regenerates it if the file doesn't exist.
pub fn load_or_generate_histogram(
    variant_counts: &HashMap<String, HashMap<usize, usize>>, 
    reference_counts: &HashMap<String, HashMap<usize, usize>>, 
    vcf_entries_by_sample: &HashMap<String, Vec<VCFEntry>>,
    num_vcfs: usize, 
    args: &Args, 
    logger: &Logger,
) -> HashMap<usize, Vec<(String, usize)>> {
    let settings_str = format!(
        "m-{}-s-{}-e-{}-z-{}",
        args.min_read_depth, args.settings, args.exclude_contig, args.restrict_contig
    );
    let outfile_path = format!("{}/site_coverage_histogram-{}.tsv", args.output_dir, settings_str);

    // Check if summary already exists
    if Path::new(&outfile_path).exists() {
        logger.information(&format!("Found existing histogram file '{}'. Loading instead of recalculating.", outfile_path));
        let histogram_positions = load_histogram_positions_from_disk(&args.output_dir, logger);
        let histogram = read_histogram_from_file(&outfile_path, logger);
        visualize_variant_site_coverage(&histogram, args, logger);
        return histogram_positions;
    }

    let sample_bases = build_sample_bases(vcf_entries_by_sample, logger);

    // Else: fallback to real calculation
    summarize_variant_site_coverage(
        variant_counts,
        reference_counts,
        &sample_bases,
        num_vcfs,
        args,
        logger,
    )
}

pub fn summarize_variant_site_coverage(
    variant_counts: &HashMap<String, HashMap<usize, usize>>, 
    reference_counts: &HashMap<String, HashMap<usize, usize>>, 
    sample_bases: &HashMap<String, HashMap<(String, usize), char>>,
    num_vcfs: usize, 
    args: &Args, 
    logger: &Logger) -> HashMap<usize, Vec<(String, usize)>> {
    logger.information(&format!("summarize_variant_site_coverage: Flatten contig/pos to global position count across {} VCFs...", num_vcfs));

    // Flatten contig/pos to global position count
    let mut position_to_count: HashMap<(String, usize), usize> = HashMap::new();
    for (contig, positions) in reference_counts {
        for (&pos, &count) in positions {
            *position_to_count.entry((contig.clone(), pos)).or_insert(0) += count;
        }
    }
    for (contig, positions) in variant_counts {
        for (&pos, &count) in positions {
            *position_to_count.entry((contig.clone(), pos)).or_insert(0) += count;
        }
    }

    logger.information(&format!("summarize_variant_site_coverage: Compute histogram..."));

    // Compute histogram
    let mut skipped_invariant_sites = 0;
    let mut histogram: Vec<(usize, usize)> = Vec::new(); // (percent, num_sites)
    let mut histogram_positions: HashMap<usize, Vec<(String, usize)>> = HashMap::new();

    for percent in (1..=100).rev() {
        let threshold = ((percent as f64 / 100.0) * num_vcfs as f64).ceil() as usize;
        let mut count = 0;

        for (&(ref contig, pos), &total_count) in &position_to_count {
            if total_count < threshold {
                continue;
            }

            // Collect bases at this site across all samples
            let mut base_set = HashSet::new();
            for sample in sample_bases.keys() {
                if let Some(base_map) = sample_bases.get(sample) {
                    if let Some(base) = base_map.get(&(contig.clone(), pos)) {
                        base_set.insert(*base);
                    }
                }
            }

            // Only include if there's variation across samples
            if base_set.len() <= 1 {
                if percent == 100 {
                    skipped_invariant_sites += 1;
                }
                continue;
            }

            count += 1;
            histogram_positions
                .entry(percent)
                .or_default()
                .push((contig.clone(), pos));
        }

        histogram.push((percent, count));
    }

    // Visualise
    visualize_variant_site_coverage(&histogram, args, logger);

    // Output skipped invariant count
    if !args.include_invariants && skipped_invariant_sites > 0 {
        logger.output(&format!(
            "\nNote: {} invariant-only positions were skipped (use --include_invariants to include them)",
            skipped_invariant_sites
        ));
    }

    histogram_positions
}

fn visualize_variant_site_coverage(histogram: &Vec<(usize, usize)>, args: &Args, logger: &Logger) {

    //logger.output("\nPhylogenetically informative site coverage:");
    //logger.output("Percent_VCFs\tNum_Sites");
    //for (percent, count) in histogram {
    //    logger.output(&format!("{}\t{}", percent, count));
    //}

    logger.output("visualize_variant_site_coverage: ASCII Plot:");
    let max_val = histogram.first().map(|(_, v)| *v).unwrap_or(1) as f64;
    for (percent, count) in histogram {
        let bar_len = (*count as f64 / max_val * 50.0).round() as usize;
        let bar = "▇".repeat(bar_len);
        logger.output(&format!("{:>3}% | {:>6} | {}", percent, count, bar));
    }

    // outfile
    let settings_str = format!("m-{}-s-{}-e-{}-z-{}", args.min_read_depth, args.settings, args.exclude_contig, args.restrict_contig);
    let outfile_path = format!("{}/site_coverage_histogram-{}.tsv", args.output_dir, settings_str);

    // Skip writing if file exists
    if Path::new(&outfile_path).exists() {
        logger.warning(&format!("Histogram file '{}' already exists — skipping write.", outfile_path));
        return;
    }

    let file = File::create(&outfile_path).unwrap_or_else(|e| {
        logger.error(&format!("Could not create histogram file '{}': {}", outfile_path, e));
        std::process::exit(1);
    });

    let mut writer = BufWriter::new(file);
    writeln!(writer, "Percent_VCFs\tNum_Sites").unwrap();
    for (percent, count) in histogram {
        writeln!(writer, "{}\t{}", percent, count).unwrap();
    }
    logger.output(&format!("Saved histogram to: {}\n", outfile_path));
}

/// Writes position files for each % threshold (1-100) into a subfolder.
pub fn write_site_position_files(histogram_positions: &HashMap<usize, Vec<(String, usize)>>, output_dir: &str, logger: &Logger) {
    let summary_dir = Path::new(output_dir).join("site_position_files");

    // Avoid re-generating if already exists
    if summary_dir.exists() {
        logger.warning("Site position subfolder already exists. Skipping generation of tab files.");
        return;
    } else {
        fs::create_dir_all(&summary_dir).expect("Failed to create site_position_files folder");
    }

    for (percent, positions) in histogram_positions {
        let file_path = summary_dir.join(format!("percent_{}.tab", percent));

        // Skip if file already exists
        if file_path.exists() {
            logger.warning(&format!("File already exists for percent {} — skipping.", percent));
            continue;
        }

        let mut writer = std::fs::File::create(&file_path).expect("Cannot write site file");
        let mut contig_map: HashMap<String, Vec<usize>> = HashMap::new();

        // Group by contig
        for (contig, pos) in positions {
            contig_map.entry(contig.clone()).or_default().push(*pos);
        }

        use std::io::Write;
        for (contig, mut positions) in contig_map {
            positions.sort();
            writeln!(writer, "##{}", contig).unwrap();
            for pos in positions {
                writeln!(writer, "{}", pos).unwrap();
            }
        }
    }
    logger.output("All site position files written to site_position_files/");
}

fn load_histogram_positions_from_disk(
    output_dir: &str,
    logger: &Logger
) -> HashMap<usize, Vec<(String, usize)>> {
    let mut histogram_positions: HashMap<usize, Vec<(String, usize)>> = HashMap::new();
    let site_dir = format!("{}/site_position_files", output_dir);

    for percent in 1..=100 {
        let file_path = format!("{}/percent_{}.tab", site_dir, percent);
        if !Path::new(&file_path).exists() {
            logger.warning(&format!("Missing file '{}'; skipping.", file_path));
            continue;
        }

        let file = File::open(&file_path).expect("Failed to open tab file");
        let reader = BufReader::new(file);
        let mut contig = String::new();

        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with("##") {
                contig = line[2..].to_string();
            } else if let Ok(pos) = line.trim().parse::<usize>() {
                histogram_positions.entry(percent).or_default().push((contig.clone(), pos));
            }
        }
    }

    logger.information("Histogram positions loaded from disk.");
    histogram_positions
}

/// Read the histogram TSV (Percent_VCFs\tNum_Sites) from file into a Vec
pub fn read_histogram_from_file(path: &str, logger: &Logger) -> Vec<(usize, usize)> {
    let file = File::open(path).unwrap_or_else(|e| {
        logger.error(&format!("Failed to open histogram file '{}': {}", path, e));
        std::process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut histogram: Vec<(usize, usize)> = Vec::new();

    for line in reader.lines().skip(1) { // Skip header
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                logger.warning(&format!("Error reading line in '{}': {}", path, e));
                continue;
            }
        };

        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() != 2 {
            logger.warning(&format!("Invalid histogram line: '{}'", line));
            continue;
        }

        let percent = parts[0].parse::<usize>().unwrap_or_else(|_| {
            logger.warning(&format!("Invalid percent in line: '{}'", line));
            0
        });

        let count = parts[1].parse::<usize>().unwrap_or_else(|_| {
            logger.warning(&format!("Invalid count in line: '{}'", line));
            0
        });

        if percent > 0 {
            histogram.push((percent, count));
        }
    }

    histogram
}