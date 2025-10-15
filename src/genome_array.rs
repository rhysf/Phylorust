use crate::logger::Logger;
use crate::read_fasta::Fasta;
use crate::read_vcf::VCFEntry;
use crate::utils::ensure_contig_header;
use crate::utils::iupac_code;
use crate::args::Args;

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};

// Use 0 only as "unset"
pub const CODE_UNSET:        u8 = 0;
pub const CODE_REFERENCE:    u8 = 1;
pub const CODE_SNP:          u8 = 2;
pub const CODE_INSERTION:    u8 = 3;
pub const CODE_DELETION:     u8 = 4;
pub const CODE_AMBIGUOUS:    u8 = 5;
pub const CODE_HET:          u8 = 6;
pub const CODE_HET_INS:      u8 = 7;
pub const CODE_HET_DEL:      u8 = 8;

#[inline]
pub fn code_for_base_type(bt: &str) -> u8 {
    match bt {
        "reference"      => CODE_REFERENCE,
        "snp"            => CODE_SNP,
        "insertion"      => CODE_INSERTION,
        "deletion"       => CODE_DELETION,
        "ambiguous"      => CODE_AMBIGUOUS,
        "heterozygous"   => CODE_HET,
        "het_insertion"  => CODE_HET_INS,
        "het_deletion"   => CODE_HET_DEL,
        _                => CODE_UNSET,
    }
}

pub fn summarise_vcf_to_tab_files_and_genome_array(
    entries: &[VCFEntry],
    fasta: &[Fasta],
    vcf_summary_subdir: &str,
    args: &Args,
    logger: &Logger) -> HashMap<String, HashMap<String, Vec<u8>>> {

    logger.information(&format!(
        "summarise_vcf_to_tab_files_and_genome_array: {} ({} entries)",
        vcf_summary_subdir,
        entries.len(),
    ));

    // 1. Setup output directory
    fs::create_dir_all(&vcf_summary_subdir).unwrap_or_else(|error| {
        logger.error(&format!("Could not create directory '{}': {}", vcf_summary_subdir, error));
        std::process::exit(1);
    });

    // 2. Setup writer structures
    let mut base_type_writers: HashMap<(String, String), BufWriter<File>> = HashMap::new();
    let mut written_contigs: HashSet<(String, String, String)> = HashSet::new();

    // 3. Pre-build sample genomes
    let mut sample_genomes: HashMap<String, HashMap<String, Vec<u8>>> = HashMap::new();
    for entry in entries {
        for sample_name in entry.samples_to_base_type.keys() {
            sample_genomes.entry(sample_name.clone()).or_insert_with(|| {
                fasta.iter()
                    .map(|f| (f.id.clone(), vec![0u8; f.seq.len()]))
                    .collect()
            });
        }
    }

    // 4. Process all variants (batched I/O, fewer lookups)
    //let mut counters: HashMap<String, usize> = HashMap::new();

    for entry in entries {
        let contig = &entry.contig;
        let pos = entry.position - 1; // convert from 1-based to 0-based

        for (sample_name, base_type) in &entry.samples_to_base_type {
            // --- Handle reference & ambiguous bases ---
            if base_type == "reference" || base_type == "ambiguous" {
                // Assign numeric code for ref/ambiguous (0 = uninitialised, 1 = ref, 4 = amb)
                let code = code_for_base_type(base_type);
                if let Some(contigs) = sample_genomes.get_mut(sample_name) {
                    if let Some(array) = contigs.get_mut(contig) {
                        if pos < array.len() {
                            array[pos] = code;
                        }
                    }
                }
                continue; // don't write .tab file entry for these
            }

            let code = code_for_base_type(base_type);

            // --- open writer lazily ---
            let key = (sample_name.clone(), base_type.clone());
            let writer = base_type_writers.entry(key.clone()).or_insert_with(|| {
                let path = format!("{}/{}-{}.tab", vcf_summary_subdir, base_type, sample_name);
                let file = File::create(&path).unwrap_or_else(|e| {
                    logger.error(&format!("Could not create '{}': {}", path, e));
                    std::process::exit(1);
                });
                BufWriter::new(file)
            });

            // --- header check ---
            ensure_contig_header(writer, &mut written_contigs, sample_name, base_type, contig);

            // --- write the variant row ---
            let (start, stop, bases_string) = make_variant_tab_row(
                entry,
                sample_name,
                base_type,
                args.heterozygosity_encoding.as_str(),
            );
            writeln!(writer, "{}\t{}\t{}", start, stop, bases_string).ok();

            // --- write numeric code to genome array ---
            if let Some(contigs) = sample_genomes.get_mut(sample_name) {
                if let Some(array) = contigs.get_mut(contig) {
                    if pos < array.len() {
                        array[pos] = code;
                    }
                }
            }

            //*counters.entry(base_type.clone()).or_insert(0) += 1;
        }
    }

    // Log summary
    // --- 6. Summarize ---
    //logger.information("Variant writing complete. Summary:");
    //for (bt, count) in &counters {
    //    logger.information(&format!("  {}: {}", bt, count));
    //}

    sample_genomes
}

fn make_variant_tab_row(
    entry: &VCFEntry,
    sample_name: &str,
    base_type: &str,
    heterozygosity_encoding: &str,
) -> (usize, usize, String) {
    // Compute start/stop range (1-based inclusive)
    let ref_len = entry.ref_base.len().max(1);
    let start = entry.position;
    let stop = entry.position + ref_len - 1;

    // Build sequence string depending on variant type
    let bases_string = match base_type {
        // --- Homozygous SNP ---
        "snp" => {
            entry.samples_to_base1
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // --- Heterozygous SNP ---
        "heterozygous" => match heterozygosity_encoding {
            "iupac" => {
                let ref_base = entry
                    .samples_to_base1
                    .get(sample_name)
                    .and_then(|s| s.chars().next())
                    .unwrap_or('N');
                let alt_base = entry
                    .samples_to_base2
                    .get(sample_name)
                    .and_then(|s| s.chars().next())
                    .unwrap_or('N');
                iupac_code(ref_base, alt_base).to_string()
            }
            _ => {
                // Default "alt" mode
                entry.samples_to_base2
                    .get(sample_name)
                    .cloned()
                    .unwrap_or_else(|| entry.alt_base.clone())
            }
        },

        // --- Homozygous insertion ---
        "insertion" => {
            entry.samples_to_base1
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // --- Heterozygous insertion ---
        "het_insertion" => {
            entry.samples_to_base2
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // --- Homozygous deletion ---
        "deletion" => {
            let anchor = entry
                .samples_to_base1
                .get(sample_name)
                .and_then(|s| s.chars().next())
                .unwrap_or('N');
            format!("{}{}", anchor, "-".repeat(ref_len.saturating_sub(1)))
        }

        // --- Heterozygous deletion ---
        "het_deletion" => {
            let anchor = entry
                .samples_to_base2
                .get(sample_name)
                .and_then(|s| s.chars().next())
                .unwrap_or('N');
            format!("{}{}", anchor, "-".repeat(ref_len.saturating_sub(1)))
        }

        _ => String::new(),
    };

    (start, stop, bases_string)
}

/// Writes contiguous genomic regions of a given base-type value, with bases as the 3rd column.
///
/// Output format:
/// 
//
//##contig
//start<TAB>end<TAB>bases
pub fn write_regions_from_genome_array(
    genome: &HashMap<String, Vec<u8>>,
    find_value: u8,
    fasta: &Vec<Fasta>,
    sample_name: &str,
    base_type: &str,
    outfile_path: &str,
    logger: &Arc<Mutex<Logger>>
) {

    logger.lock().unwrap().information(&format!("write_regions_from_genome_array: finding contiguous regions for value {}...", find_value));

    let file = File::create(outfile_path).unwrap_or_else(|e| {
        logger.lock().unwrap().error(&format!("Could not create file '{}': {}", outfile_path, e));
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    let mut count_regions = 0;
    let mut count_bases = 0;

    // Track whether we've already written the header for a contig
    let mut written_contigs: HashSet<(String, String, String)> = HashSet::new();

    for (contig, values) in genome {
        let mut start: Option<usize> = None;
        let length = values.len();

        // Find contiguous runs of `find_value`
        for i in 0..=length {
            let current = if i < length { values[i] } else { 0 }; // Sentinel: fake zero to close final run

            if current != find_value {
                if let Some(start_pos) = start {
                    let end = i;
                    count_regions += 1;
                    count_bases += end - start_pos;

                    // Write header only once per contig
                    ensure_contig_header(&mut writer, &mut written_contigs, sample_name, base_type, contig);

                    // Extract the reference sequence
                    let seq = fasta
                        .iter()
                        .find(|f| f.id == *contig)
                        .map(|f| {
                            f.seq
                                .chars()
                                .skip(start_pos)
                                .take(end - start_pos)
                                .collect::<String>()
                        })
                        .unwrap_or_else(|| "N".repeat(end - start_pos));

                    // Write region to file
                    writeln!(writer, "{}\t{}\t{}", start_pos + 1, end, seq).ok();

                    start = None;
                }
            } else if start.is_none() {
                start = Some(i);
            }
        }
    }

    // "write_regions_from_genome_array: wrote {} regions covering {} bases to {}", count_regions, count_bases, outfile_path
    logger.lock().unwrap().information(&format!("write_regions_from_genome_array: wrote {} regions covering {} bases", count_regions, count_bases));
}