use crate::logger::Logger;
use std::fs::write;
use std::fs::File;
use std::path::Path;
use std::process::{Command, Stdio};
use which::which;
use std::collections::{HashSet};
use std::io::{Write};
use std::sync::{Arc, Mutex};

pub fn check_r_installed() -> bool {
    Command::new("Rscript")
        .arg("--version")
        .output()
        .is_ok()
}

/// Run FastTree on a set of FASTA files (one per percent threshold).
///
/// - `output_dir`: where FASTAs and trees live
/// - `percents`: list of percent thresholds (e.g., [90, 95, 100])
/// - `logger`: logging utility
/// - `fasttree_bin`: optional path to FastTree binary. If `None`, will search PATH.
pub fn run_fasttree_on_fastas(
    output_dir: &str,
    percents: &[usize],
    settings_str: &str,
    logger: &Arc<Mutex<Logger>>,
    fasttree_bin: Option<&str>,
) {
    // Decide which binary to use
    let fasttree_exe = if let Some(path) = fasttree_bin {
        Path::new(path).to_path_buf()
    } else if let Ok(path) = which("FastTree") {
        path
    } else {
        logger.lock().unwrap().warning(
            "FastTree not found (neither in PATH nor provided with --fasttree-bin). Skipping tree generation.",
        );
        return;
    };

    logger.lock().unwrap().information(&format!("Using FastTree binary: {}", fasttree_exe.display()));

    for percent in percents {
        let fasta_file = format!("{}/percent_{}-{}.fasta", output_dir, percent, settings_str);
        let tree_file = format!("{}/percent_{}-{}-FastTree.tree", output_dir, percent, settings_str);

        if !Path::new(&fasta_file).exists() {
            logger.lock().unwrap().warning(&format!(
                "FASTA file not found for {}%: {}",
                percent, fasta_file
            ));
            continue;
        }

        logger.lock().unwrap().information(&format!("Running FastTree on {}", fasta_file));

        let tree_out = match File::create(&tree_file) {
            Ok(f) => f,
            Err(e) => {
                logger.lock().unwrap().error(&format!("Failed to create tree output file: {}", e));
                continue;
            }
        };

        match Command::new(&fasttree_exe)
            .arg("-nt")
            .arg(&fasta_file)
            .stdout(Stdio::from(tree_out))
            .status()
        {
            Ok(s) if s.success() => {
                logger.lock().unwrap().information(&format!("Tree written to {}", tree_file));

                // Read Newick and pretty-print ASCII tree
                match std::fs::read_to_string(&tree_file) {
                    Ok(newick) => {
                        logger.lock().unwrap().information(&format!("ASCII tree for {}%:", percent));
                        print_newick_tree(&newick, logger);
                    }
                    Err(e) => {
                        logger.lock().unwrap().warning(&format!("Could not read {}: {}", tree_file, e));
                    }
                }
            }
            Ok(s) => {
                logger.lock().unwrap().warning(&format!("FastTree failed with exit code: {}", s));
            }
            Err(e) => {
                logger.lock().unwrap().error(&format!("Failed to run FastTree: {}", e));
            }
        }
    }
}

/// Print a phylogenetic tree from a Newick string as an indented hierarchy.
///
/// # Arguments
/// * `newick` - A string slice containing a Newick tree.
/// * `logger` - Your logger instance for output.
///
/// Example Newick:
/// (A:0.1,(B:0.2,C:0.3):0.4,D:0.5);
pub fn print_newick_tree(newick: &str, logger: &Arc<Mutex<Logger>>) {
    // Strip trailing semicolon if present
    let newick = newick.trim_end_matches(';');

    fn recurse(subtree: &str, depth: usize, logger: &Arc<Mutex<Logger>>) {
        // Split by top-level commas
        let mut balance = 0;
        let mut parts = Vec::new();
        let mut start = 0;

        for (i, ch) in subtree.char_indices() {
            match ch {
                '(' => balance += 1,
                ')' => balance -= 1,
                ',' if balance == 0 => {
                    parts.push(&subtree[start..i]);
                    start = i + 1;
                }
                _ => {}
            }
        }
        parts.push(&subtree[start..]);

        // Print each part
        for part in parts {
            let part = part.trim();

            if part.starts_with('(') {
                // Nested clade
                let inner = part.trim_start_matches('(').trim_end_matches(|c| c == ')' || c == ':');
                logger.lock().unwrap().output(&format!("{}[clade]", "  ".repeat(depth)));
                recurse(inner, depth + 1, logger);
            } else {
                // Leaf: may contain a branch length (e.g. "A:0.1")
                let name = part.split(':').next().unwrap_or("").trim();
                if !name.is_empty() {
                    logger.lock().unwrap().output(&format!("{}{}", "  ".repeat(depth), name));
                }
            }
        }
    }

    recurse(newick, 0, logger);
}

pub fn run_r_plotting_script(histogram_file: &str, percent_for_tree: usize, output_dir: &str, logger: &Arc<Mutex<Logger>>) {
    if !check_r_installed() {
        logger.lock().unwrap().warning("run_r_plotting_script: R not found. Skipping graphical histogram generation.");
        return;
    }

    logger.lock().unwrap().information("run_r_plotting_script: Writing Rscript...");

    let r_script = r#"
#!/usr/bin/env Rscript

suppressMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_histogram.R <input.tsv> <percent_for_tree>")
}

input_file <- args[1]
selected_percent <- as.numeric(args[2])
output_prefix <- tools::file_path_sans_ext(input_file)

# Read histogram data
data <- read_tsv(input_file, col_types = cols(
  Percent_VCFs = col_integer(),
  Num_Sites = col_integer()
))

# Filter out invalid or out-of-bounds rows
data <- data %>%
  filter(
    !is.na(Percent_VCFs),
    !is.na(Num_Sites),
    Percent_VCFs >= 0,
    Percent_VCFs <= 100
  )

# Aggregate into 10% bins with numeric midpoints
bin_breaks <- seq(0, 100, 10)

data <- data %>%
  mutate(
    BinIndex = cut(Percent_VCFs,
                   breaks = bin_breaks,
                   include.lowest = TRUE,
                   labels = FALSE)
  ) %>%
  group_by(BinIndex) %>%
  summarise(
    Num_Sites = mean(Num_Sites),
    BinMid = mean(Percent_VCFs),   # midpoint from actual values
    .groups = "drop"
  )

# Compute Y-axis padding
y_max <- max(data$Num_Sites, na.rm = TRUE)
y_pad <- ceiling(y_max * 0.1)
y_limit <- y_max + y_pad

# Get y value for red dot, if available
bin_for_dot <- as.character(cut(
  selected_percent,
  breaks = seq(0, 100, 5),
  include.lowest = TRUE
))

# Main plot (, color = "white", width = 0.9 to geom_col
p <- ggplot(data, aes(x = BinMid, y = Num_Sites)) +
  geom_col(fill = "steelblue", color = "white", width = 8) +
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 100),
    expand = c(0, 5)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 0)
  )

# Switch to “thousands” view
if (y_max > 10000) {
  # scale values to thousands
  p <- p + scale_y_continuous(
    limits = c(0, y_limit),
    expand = c(0, 0),
    labels = function(x) x / 1000
  )
  y_label <- "Number of sites (thousands)"
} else {
  p <- p + scale_y_continuous(
    limits = c(0, y_limit),
    expand = c(0, 0)
  )
  y_label <- "Number of sites"
}

# Add labels (after y_label is set)
p <- p + labs(
  title = "Phylogenetically informative sites",
  x = "Percent of VCFs (10% bins)",
  y = y_label
)

# Add the red dot above the corresponding bar (if found)
dot_y <- data %>%
  filter(BinMid <= selected_percent & selected_percent < BinMid + 10) %>%
  pull(Num_Sites)

if (length(dot_y) == 1 && !is.na(dot_y)) {
  # Add padding so the dot is slightly above the bar
  dot_y_pos <- dot_y + y_max * 0.02

  p <- p + annotate("point", x = selected_percent, y = dot_y_pos, color = "red", size = 3)
}

# Save plots
ggsave(paste0(output_prefix, ".png"), p, width = 8, height = 8, dpi = 300)
ggsave(paste0(output_prefix, ".pdf"), p, width = 8, height = 8)

message("Saved plots to: ", output_prefix, ".png and .pdf")
"#;

    let script_path = format!("{}/plot_histogram.R", output_dir);
    if let Err(e) = write(&script_path, r_script) {
        logger.lock().unwrap().error(&format!("Failed to write R script: {}", e));
        return;
    }

    logger.lock().unwrap().information("Rscript detected. Generating graphical histogram...");

    let status = Command::new("Rscript")
        .arg(script_path)
        .arg(histogram_file)
        .arg(percent_for_tree.to_string())
        .status();

    match status {
        Ok(s) if s.success() => logger.lock().unwrap().information("Histogram plot generated successfully."),
        Ok(s) => logger.lock().unwrap().warning(&format!("Rscript failed with exit code: {}", s)),
        Err(e) => logger.lock().unwrap().error(&format!("Failed to run Rscript: {}", e)),
    }

    // Optional: clean up
    // let _ = std::fs::remove_file(script_path);
}

/// Ensures that a header line "##{contig}" is written once per contig.
///
/// Returns `true` if a new header was written, `false` if it was already present.
pub fn ensure_contig_header(
    writer: &mut std::io::BufWriter<std::fs::File>,
    written_contigs: &mut HashSet<(String, String, String)>,
    sample: &str,
    base_type: &str,
    contig: &str,
) {
    let key = (sample.to_string(), base_type.to_string(), contig.to_string());
    if !written_contigs.contains(&key) {
        writeln!(writer, "##{}", contig).ok();
        written_contigs.insert(key);
    }
}

/// Converts a heterozygous pair (ref, alt) into an IUPAC ambiguity code.
pub fn iupac_code(ref_base: char, alt_base: char) -> char {
    match (ref_base, alt_base) {
        ('A', 'G') | ('G', 'A') => 'R',
        ('C', 'T') | ('T', 'C') => 'Y',
        ('A', 'C') | ('C', 'A') => 'M',
        ('G', 'T') | ('T', 'G') => 'K',
        ('A', 'T') | ('T', 'A') => 'W',
        ('C', 'G') | ('G', 'C') => 'S',
        // If same base or invalid, return one of them
        (a, b) if a == b => a,
        _ => 'N',
    }
}