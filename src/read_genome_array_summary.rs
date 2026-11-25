use crate::logger::Logger;
use std::collections::{BTreeMap, HashMap};
use crate::read_tab::SampleCalls;
use std::fs;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use crate::args::Args;
use std::sync::{Arc, Mutex};

/// Either loads the histogram + site position data from disk,
/// or regenerates it if the file doesn't exist.
pub fn load_or_generate_histogram(
    sample_calls: &SampleCalls,
    rustatools_settings_string: &str,
    args: &Args, 
    rundir: &str,
    logger: &Arc<Mutex<Logger>>,
) -> HashMap<usize, Vec<(String, usize)>> {

    let outfile_path = format!("{}/site_coverage_histogram-{}.tsv", rundir, rustatools_settings_string);
    logger.lock().unwrap().information(&format!("load_or_generate_histogram: Output file = '{}'", outfile_path));

    // Check if summary already exists
    let site_dir = Path::new(rundir).join(format!("site_positions-{rustatools_settings_string}"));
    if Path::new(&outfile_path).exists()
        && (1..=100).all(|percent| site_dir.join(format!("percent_{percent}.tab")).is_file())
    {
        logger.lock().unwrap().information(&format!("load_or_generate_histogram: Found existing histogram file '{}'. Loading instead of recalculating.", outfile_path));
        let histogram_positions = load_histogram_positions_from_disk(&rundir, &rustatools_settings_string, &logger);
        let histogram = read_histogram_from_file(&outfile_path, &logger);
        visualize_variant_site_coverage(&histogram, &rustatools_settings_string, &rundir, &logger);
        return histogram_positions;
    }

    // Else: fallback to real calculation
    summarise_variant_site_coverage(
        sample_calls,
        &rustatools_settings_string,
        args,
        &rundir,
        &logger,
    )
}

struct SiteCalls {
    count: usize,
    first_base: u8,
    variable: bool,
}

fn summarise_variant_site_coverage(
    sample_calls: &SampleCalls,
    settings: &str,
    args: &Args,
    rundir: &str,
    logger: &Arc<Mutex<Logger>>,
) -> HashMap<usize, Vec<(String, usize)>> {
    let mut called_sites: BTreeMap<&str, BTreeMap<usize, SiteCalls>> = BTreeMap::new();
    for contigs in sample_calls.values() {
        for (contig, calls) in contigs {
            let sites = called_sites.entry(contig).or_default();
            for (&position, &base) in calls {
                let site = sites.entry(position).or_insert(SiteCalls {
                    count: 0,
                    first_base: base,
                    variable: false,
                });
                site.count += 1;
                site.variable |= site.first_base != base;
            }
        }
    }

    let mut histogram = Vec::new();
    let mut histogram_positions = HashMap::new();
    let mut skipped_invariant_sites = 0;
    for percent in 1..=100 {
        let threshold = (percent * sample_calls.len()).div_ceil(100);
        let mut positions = Vec::new();
        for (&contig, sites) in &called_sites {
            for (&position, site) in sites {
                if site.count < threshold {
                    continue;
                }
                if site.variable || args.include_invariants {
                    positions.push((contig.to_owned(), position));
                } else if percent == 100 {
                    skipped_invariant_sites += 1;
                }
            }
        }
        histogram.push((percent, positions.len()));
        histogram_positions.insert(percent, positions);
    }

    visualize_variant_site_coverage(&histogram, settings, rundir, logger);
    if skipped_invariant_sites > 0 {
        logger.lock().unwrap().output(&format!(
            "\n{skipped_invariant_sites} fully called candidate positions were invariant."
        ));
    }

    histogram_positions
}

fn visualize_variant_site_coverage(
    histogram: &Vec<(usize, usize)>, 
    rustatools_settings_string: &str,
    rundir: &str, 
    logger: &Arc<Mutex<Logger>>) {
    logger.lock().unwrap().output("visualize_variant_site_coverage: ASCII Plot...");

    // Find max value safely
    let max_val = histogram.first().map(|(_, v)| *v).unwrap_or(0).max(1) as f64;

    for (percent, count) in histogram {
        // Show 100%, 1%, and every 10% step
        if *percent == 100 || *percent == 1 || percent % 10 == 0 {
            let mut bar_len = (*count as f64 / max_val * 50.0).round() as usize;

            // Cap at 200 characters to prevent runaway allocation
            if bar_len > 200 {
                logger.lock().unwrap().warning(&format!(
                    "Bar length too large ({}), capping at 200 (percent={}, count={}, max_val={})",
                    bar_len, percent, count, max_val
                ));
                bar_len = 200;
            }

            let bar = "▇".repeat(bar_len);
            logger.lock().unwrap().output(&format!("{:>3}% | {:>6} | {}", percent, count, bar));
        }
    }

    // outfile
    let outfile_path = format!("{}/site_coverage_histogram-{}.tsv", rundir, rustatools_settings_string);

    // Skip writing if file exists
    if Path::new(&outfile_path).exists() {
        logger.lock().unwrap().warning(&format!("Histogram file '{}' already exists — skipping write.", outfile_path));
        return;
    }

    let pending = Path::new(&outfile_path).with_extension(format!("partial-{}", std::process::id()));
    let file = File::create(&pending).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Could not create histogram file '{}': {}", outfile_path, e));
        std::process::exit(1);
    });

    let mut writer = BufWriter::new(file);
    writeln!(writer, "Percent_VCFs\tNum_Sites").unwrap();
    for (percent, count) in histogram {
        writeln!(writer, "{}\t{}", percent, count).unwrap();
    }
    writer.flush().unwrap();
    fs::rename(pending, &outfile_path).expect("Cannot publish histogram");
    logger.lock().unwrap().output(&format!("Saved histogram to: {}\n", outfile_path));
}

/// Writes position files for each % threshold (1-100) into a subfolder.
pub fn write_site_position_files(
    histogram_positions: &HashMap<usize, Vec<(String, usize)>>, 
    output_dir: &str, 
    settings_str: &str,
    logger: &Arc<Mutex<Logger>>) {

    let summary_dir = Path::new(output_dir).join(format!("site_positions-{}", settings_str));
    logger.lock().unwrap().information("write_site_position_files: summary_dir");

    fs::create_dir_all(&summary_dir).expect("Failed to create site position folder");

    for (percent, positions) in histogram_positions {
        let file_path = summary_dir.join(format!("percent_{}.tab", percent));

        // Skip if file already exists
        if file_path.exists() {
            continue;
        }

        let pending = file_path.with_extension(format!("partial-{}", std::process::id()));
        let mut writer = BufWriter::new(File::create(&pending).expect("Cannot write site file"));
        let mut contig_map: BTreeMap<&str, Vec<usize>> = BTreeMap::new();

        // Group by contig
        for (contig, pos) in positions {
            contig_map.entry(contig).or_default().push(*pos);
        }

        for (contig, mut positions) in contig_map {
            positions.sort();
            writeln!(writer, "##{}", contig).unwrap();
            for pos in positions {
                writeln!(writer, "{}", pos).unwrap();
            }
        }
        writer.flush().unwrap();
        fs::rename(pending, file_path).expect("Cannot publish site file");
    }
    logger.lock().unwrap().output("All site position files written to site_position_files/");
}

fn load_histogram_positions_from_disk(
    output_dir: &str,
    settings_str: &str,
    logger: &Arc<Mutex<Logger>>
) -> HashMap<usize, Vec<(String, usize)>> {
    let mut histogram_positions: HashMap<usize, Vec<(String, usize)>> = HashMap::new();
    let site_dir = format!("{}/site_positions-{}", output_dir, settings_str);

    for percent in 1..=100 {
        let file_path = format!("{}/percent_{}.tab", site_dir, percent);
        if !Path::new(&file_path).exists() {
            logger.lock().unwrap().warning(&format!("Missing file '{}'; skipping.", file_path));
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

    logger.lock().unwrap().information("Histogram positions loaded from disk.");
    histogram_positions
}

/// Read the histogram TSV (Percent_VCFs\tNum_Sites) from file into a Vec
pub fn read_histogram_from_file(path: &str, logger: &Arc<Mutex<Logger>>) -> Vec<(usize, usize)> {
    let file = File::open(path).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Failed to open histogram file '{}': {}", path, e));
        std::process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut histogram: Vec<(usize, usize)> = Vec::new();

    for line in reader.lines().skip(1) { // Skip header
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                logger.lock().unwrap().warning(&format!("Error reading line in '{}': {}", path, e));
                continue;
            }
        };

        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() != 2 {
            logger.lock().unwrap().warning(&format!("Invalid histogram line: '{}'", line));
            continue;
        }

        let percent = parts[0].parse::<usize>().unwrap_or_else(|_| {
            logger.lock().unwrap().warning(&format!("Invalid percent in line: '{}'", line));
            0
        });

        let count = parts[1].parse::<usize>().unwrap_or_else(|_| {
            logger.lock().unwrap().warning(&format!("Invalid count in line: '{}'", line));
            0
        });

        if percent > 0 {
            histogram.push((percent, count));
        }
    }

    histogram
}