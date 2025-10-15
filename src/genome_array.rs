use crate::logger::Logger;
use crate::read_fasta::Fasta;
use crate::read_vcf::VCFEntry;
use crate::utils::ensure_contig_header;
use crate::utils::iupac_code;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn summarise_vcf_to_tab_files_and_genome_array(
    logger: &Logger, 
    entries: &Vec<VCFEntry>, 
    fasta: &Vec<Fasta>, 
    vcf_summary_subdir: &str,
    heterozygosity_encoding: &str) -> (
        HashMap<String, HashMap<String, Vec<u8>>>, // sample → contig → array
        HashMap<String, u8>,                        // base_type → code{
    ) {

    logger.information("summarise_vcf_to_tab_files_and_genome_array: filling genome arrays for each sample...");

    // sample_name → (contig → Vec<u8>)
    let mut sample_genomes: HashMap<String, HashMap<String, Vec<u8>>> = HashMap::new();

    // base_type → numeric code
    let mut type_to_code: HashMap<String, u8> = HashMap::new();
    let mut next_code: u8 = 1;

    // --- NEW: prepare per-base-type writers ---
    let mut base_type_writers: HashMap<(String, String), BufWriter<File>> = HashMap::new();
    let mut last_contig_written: HashMap<(String, String), String> = HashMap::new();
    let mut written_contigs: HashSet<(String, String, String)> = HashSet::new();

    // Pre-initialize arrays for all samples & contigs
    for entry in entries {
        for sample_name in entry.samples_to_base_type.keys() {
            sample_genomes.entry(sample_name.clone()).or_insert_with(|| {
                let mut contigs = HashMap::new();
                for f in fasta {
                    contigs.insert(f.id.clone(), vec![0u8; f.seq.len()]);
                }
                contigs
            });
        }
    }

    // Fill arrays with variant codes
    for entry in entries {
        let contig = &entry.contig;
        let pos = entry.position - 1; // convert from 1-based to 0-based

        for (sample_name, base_type) in &entry.samples_to_base_type {
            // Dynamically assign numeric codes to new base types
            let code = *type_to_code.entry(base_type.clone()).or_insert_with(|| {
                let c = next_code;
                next_code += 1;
                c
            });

            // --- NEW: Handle variants with true ALT bases ---
            match base_type.as_str() {
                "reference" | "ambiguous" => {
                    // We'll handle these later via contiguous-region detection
                }
                _ => {
                    write_variant_entry(
                        &mut base_type_writers,
                        &mut last_contig_written,
                        &mut written_contigs,
                        logger,
                        sample_name,
                        base_type,
                        &entry.contig,
                        entry,
                        heterozygosity_encoding,
                        vcf_summary_subdir,
                    );
                }
            }

            // --- Continue marking the genome array as before ---

            // Write code into the appropriate sample + contig + position
            if let Some(contigs) = sample_genomes.get_mut(sample_name) {
                if let Some(array) = contigs.get_mut(contig) {
                    if pos < array.len() {
                        array[pos] = code;
                    } else {
                        logger.warning(&format!(
                            "Position {} out of bounds for contig '{}' (len={})",
                            pos, contig, array.len()
                        ));
                    }
                }
            }
        }
    }

    // --- Flush all open variant writers ---
    for ((_sample, _btype), mut writer) in base_type_writers {
        writer.flush().ok();
    }

    // Log summary
    logger.information(&format!(
        "summarise_vcf_to_tab_files_and_genome_array: processed {} entries, {} base types, {} samples.",
        entries.len(),
        type_to_code.len(),
        sample_genomes.len()
    ));

    //for (t, c) in &type_to_code {
    //    logger.information(&format!("  code {:>2} → {}", c, t));
    //}

    (sample_genomes, type_to_code)
}

fn write_variant_entry(
    base_type_writers: &mut HashMap<(String, String), BufWriter<File>>,
    last_contig_written: &mut HashMap<(String, String), String>,
    written_contigs: &mut HashSet<(String, String, String)>,
    logger: &Logger,
    sample_name: &str,
    base_type: &str,
    contig: &str,
    entry: &VCFEntry,
    heterozygosity_encoding: &str,
    vcf_summary_subdir: &str,) {

    // 1) Get (or open) writer
    let key = (sample_name.to_string(), base_type.to_string());
    let writer = base_type_writers.entry(key.clone()).or_insert_with(|| {
        let path = format!("{}/{}-{}.tab", vcf_summary_subdir, base_type, sample_name);
        let file = File::create(&path).unwrap_or_else(|e| {
            logger.error(&format!("Could not create '{}': {}", path, e));
            std::process::exit(1);
        });
        BufWriter::new(file)
    });

    // 2) contig header once per (sample,base_type,contig)
    if last_contig_written.get(&key) != Some(&contig.to_string()) {
        ensure_contig_header(writer, written_contigs, sample_name, base_type, contig);
        last_contig_written.insert(key.clone(), contig.to_string());
    }

    // 3) compute span (1-based inclusive)
    let ref_len = entry.ref_base.len().max(1);
    let start = entry.position;
    let stop = entry.position + ref_len - 1;

    // 4) bases string per type
    let bases_string = match base_type {
        // Homozygous SNP
        "snp" => {
            entry.samples_to_base1
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // Heterozygous SNP — use base2 (or optionally compute IUPAC later)
        "heterozygous" => {
            match heterozygosity_encoding {
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
                    // "alt" mode — simply return the alt base
                    entry.samples_to_base2
                        .get(sample_name)
                        .cloned()
                        .unwrap_or_else(|| entry.alt_base.clone())
                }
            }
        }

        // Homozygous insertion — inserted bases after anchor
        "insertion" => {
            entry.samples_to_base1
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // Heterozygous insertion — the other allele’s inserted bases
        "het_insertion" => {
            entry.samples_to_base2
                .get(sample_name)
                .cloned()
                .unwrap_or_else(|| entry.alt_base.clone())
        }

        // Homozygous deletion — anchor + '-' padding
        "deletion" => {
            let ref_len = entry.ref_base.len();
            let anchor = entry
                .samples_to_base1
                .get(sample_name)
                .and_then(|s| s.chars().next())
                .unwrap_or('N');
            format!("{}{}", anchor, "-".repeat(ref_len.saturating_sub(1)))
        }

        // Heterozygous deletion — base2 + '-' padding
        "het_deletion" => {
            let ref_len = entry.ref_base.len();
            let anchor = entry
                .samples_to_base2
                .get(sample_name)
                .and_then(|s| s.chars().next())
                .unwrap_or('N');
            format!("{}{}", anchor, "-".repeat(ref_len.saturating_sub(1)))
        }

        _ => String::new(),
    };

    // 5) write line
    writeln!(writer, "{}\t{}\t{}", start, stop, bases_string).ok();
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
    logger: &Logger,
) {

    logger.information(&format!("write_regions_from_genome_array: finding contiguous regions for value {}...", find_value));

    let file = File::create(outfile_path).unwrap_or_else(|e| {
        logger.error(&format!("Could not create file '{}': {}", outfile_path, e));
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
    logger.information(&format!("write_regions_from_genome_array: wrote {} regions covering {} bases", count_regions, count_bases));
}