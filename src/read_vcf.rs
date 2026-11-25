use crate::logger::Logger;
use crate::read_tab::NameTypeLocation;
use crate::args::Args;
use std::fs::{self};
use std::process;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

pub struct VCFsamples {
    //header:String,
    samples:HashMap<usize, String>,
}

#[derive(Debug, Clone)]
pub struct VCFEntry {
    pub contig:String,
    pub position:usize,
    //pub id:String,
    pub ref_base:String,
    pub alt_base:String,
    //pub cons_qual:String,
    //pub filter:String,
    //pub info:String,
    //pub format:String,
    //pub samples:HashMap<usize, String>,
    //pub samples_to_genotype:HashMap<String, String>,
    pub samples_to_base1:HashMap<String, String>,
    pub samples_to_base2:HashMap<String, String>,
    pub samples_to_base_type:HashMap<String, String>,
}

pub fn read_vcf(
    name_type_location : &NameTypeLocation, 
    global_sample_names: &mut HashSet<String>, 
    args: &Args,
    logger: &Arc<Mutex<Logger>>, ) -> (Vec<VCFEntry>, Vec<String>) {

    logger.lock().unwrap().information(&format!("read VCF: {}", name_type_location.location));

    // read file
    let vcf_file = fs::read_to_string(&name_type_location.location).unwrap_or_else(|error|{
        logger.lock().unwrap().error(&format!("Error with file: {} {}", name_type_location.location, error));
        process::exit(1);
    });

    let mut vcf_entries: Vec<VCFEntry> = Vec::new();
    let mut vcf_samples:Option<VCFsamples> = None;

    // go through each line of the VCF
    for line in vcf_file.lines() {

        // Metadata
        if line.starts_with("##") { 
            continue; 
        }
        // Header
        if line.starts_with("#CHROM") { 

            let line_parts : Vec<&str> = line.split('\t').collect();
            if line_parts.len() < 10 {
                logger.lock().unwrap().error("VCF header line does not contain sample columns.");
                process::exit(1);
            }
            let mut sample_names = HashMap::new();

            //logger.information(&format!("Header line: {}", line));
            //logger.information(&format!("Parsed sample names: {:?}", &line_parts[9..]));
            for(index, sample_name_raw) in line_parts[9..].iter().enumerate() {
                let mut sample_name = sample_name_raw.to_string();

                // Make name unique if it's already been used
                if global_sample_names.contains(&sample_name) {
                    let mut counter = 1;
                    let original = sample_name.clone();
                    while global_sample_names.contains(&sample_name) {
                        sample_name = format!("{}-{}", original, counter);
                        counter += 1;
                    }

                    logger.lock().unwrap().warning(&format!(
                        "Duplicate sample name detected: Sample '{}' renamed to '{}'.",
                        original, sample_name
                    ));
                }

                global_sample_names.insert(sample_name.clone());
                sample_names.insert(index, sample_name);
            }

            vcf_samples = Some(VCFsamples{
                //header:line.to_string(),
                samples:sample_names,
            });
            continue; 
        }

        // Entries
        let vcf_samples_ref = match vcf_samples.as_ref() {
            Some(samples) => samples,
            None => {
                logger.lock().unwrap().error("read_vcf: VCF header line missing.");
                process::exit(1);
            }
        };
        
        if let Some(entry) = read_vcf_line(line, &args, logger, vcf_samples_ref) {
            vcf_entries.push(entry);
        }
    }
    let samples = vcf_samples.unwrap_or_else(|| {
        logger.lock().unwrap().error("read_vcf: VCF header line missing.");
        process::exit(1);
    });
    (vcf_entries, samples.samples.into_values().collect())
}

fn read_vcf_line(line : &str, args: &Args, logger: &Arc<Mutex<Logger>>, samples: &VCFsamples) -> Option<VCFEntry> {

    let line_parts : Vec<&str> = line.split('\t').collect();
    if line_parts.len() < 10 {
        logger.lock().unwrap().error(&format!("Error with vcf line: {line}"));
        process::exit(1);
    }
    let contig = line_parts[0].to_string();
    let position: usize = line_parts[1].parse::<usize>().unwrap_or_else(|_| {
        logger.lock().unwrap().warning(&format!("Invalid position value in VCF: '{}'", line_parts[1]));
        process::exit(1);
    });
    //let id = line_parts[2].to_string();
    let ref_base = line_parts[3].to_ascii_uppercase();
    let alt_base = line_parts[4].to_ascii_uppercase();
    //let cons_qual = line_parts[5].to_string();
    //let filter = line_parts[6].to_string();
    //let info = line_parts[7].to_string();
    //let format = line_parts[8].to_string();

    // Ignore contig
    if args.exclude_contig != "n" && contig == args.exclude_contig {
        return None;
    }
    // Restrict to contig
    if args.restrict_contig != "n" && contig != args.restrict_contig {
        return None;
    }

    // get sample info
    //let sample_names = samples.samples.clone();

    let (gt_index, depth_indices) = read_genotype_and_depth_positions(line_parts[8], logger);
    let samples_to_genotype = create_sample_to_genotype(
        &line_parts[9..], gt_index, depth_indices, &samples.samples, args.min_read_depth,
    );
    if samples_to_genotype.is_empty() {
        return None;
    }

    let (samples_to_base1, samples_to_base2, samples_to_base_type) = determine_bases_and_base_type(
        &samples_to_genotype,
        &ref_base,
        &alt_base,
        logger,
    );

    // return various info for vcf entry
    Some(VCFEntry{
        contig,
        position,
        //id,
        ref_base,
        alt_base,
        //cons_qual,
        //filter,
        //info,
        //format,
        //samples:sample_names,
        //samples_to_genotype,
        samples_to_base1,
        samples_to_base2,
        samples_to_base_type,
    })
}

fn read_genotype_and_depth_positions(
    format: &str,
    logger: &Arc<Mutex<Logger>>,
) -> (usize, [Option<usize>; 2]) {
    let mut gt = None;
    let mut dp = None;
    let mut min_dp = None;
    for (index, key) in format.split(':').enumerate() {
        match key.to_ascii_uppercase().as_str() {
            "GT" => gt = Some(index),
            "DP" => dp = Some(index),
            "MIN_DP" => min_dp = Some(index),
            _ => {}
        }
    }
    let gt = gt.unwrap_or_else(|| {
        logger.lock().unwrap().error(&format!("Error with vcf line (no genotype found): {format}"));
        process::exit(1);
    });
    (gt, [dp, min_dp])
}

fn create_sample_to_genotype(
    samples: &[&str],
    gt_index: usize,
    depth_indices: [Option<usize>; 2],
    sample_names: &HashMap<usize, String>,
    min_depth: u8,
) -> HashMap<String, String> {
    let mut genotypes = HashMap::with_capacity(samples.len());
    let filter_depth = min_depth > 0 && depth_indices.iter().any(Option::is_some);

    for (index, sample_field) in samples.iter().enumerate() {
        let fields: Vec<&str> = sample_field.split(':').collect();
        let Some(&genotype) = fields.get(gt_index) else { continue };

        if filter_depth {
            // Prefer DP; expanded gVCF reference rows can provide only MIN_DP.
            let depth = depth_indices.iter().flatten().find_map(|&field| {
                fields.get(field).and_then(|value| value.parse::<usize>().ok())
            });
            if !depth.is_some_and(|depth| depth >= usize::from(min_depth)) {
                continue;
            }
        }

        let genotype = match genotype.split_once(['/', '|']) {
            Some((first, second)) if first == second => first,
            _ => genotype,
        };
        genotypes.insert(sample_names[&index].clone(), genotype.to_owned());
    }
    genotypes
}

fn determine_bases_and_base_type(
    samples_to_genotype: &HashMap<String, String>, 
    ref_base : &String, 
    alt_base : &String, 
    logger: &Arc<Mutex<Logger>>
) -> (HashMap<String, String>, HashMap<String, String>, HashMap<String, String>) {
    let mut sample_to_base1 = HashMap::new();
    let mut sample_to_base2 = HashMap::new();
    let mut sample_to_base_type = HashMap::new();

    // loop over each sample and fill hash maps
    for (sample, genotype) in samples_to_genotype.iter() {
        let (base1, base2, base_type) = determine_bases_and_base_type_single_sample(&genotype, &ref_base, &alt_base, &logger);
        sample_to_base1.insert(sample.to_string(), base1);
        sample_to_base2.insert(sample.to_string(), base2);
        sample_to_base_type.insert(sample.to_string(), base_type);
    }

    (sample_to_base1, sample_to_base2, sample_to_base_type)
}

fn determine_bases_and_base_type_single_sample(genotype: &String, ref_base : &String, alt_base : &String, logger: &Arc<Mutex<Logger>>) -> (String, String, String) {

    // Save bases and GT parts
    let consensus_bases : Vec<&str> = alt_base.split(',').collect();
    let mut all_bases = vec![ref_base.as_str()];
    all_bases.extend(consensus_bases.iter().cloned());

    //logger.lock().unwrap().warning("determine_bases_and_base_type_single_sample");

    let gt_parts : Vec<&str> = genotype.split(|c| c == '/' || c == '|').collect();
    let mut gt_bases : Vec<&str> = Vec::new();
    for gt_part in gt_parts.iter() {

        // Try to parse genotype index
        let n = match gt_part.parse::<usize>() {
            Ok(n) => n,
            Err(_) => {
                return ("N".to_string(), String::new(), "ambiguous".to_string());
            }
        };

        // Validate index range
        if all_bases.len() <= n {
            logger.lock().unwrap().error(&format!("vcf parsing genotype error: No {} defined from ref {} and cons {}", n, ref_base, alt_base));
            process::exit(1);
        }

        let gt_base = all_bases[n];

        // Ambiguous bases
        if matches!(gt_base, "N" | "." | "*") || gt_base.starts_with('<') {
            return ("N".to_string(), String::new(), "ambiguous".to_string());
        }

        // save the bases
        gt_bases.push(gt_base);
    }

	// Homozygous ref-calls
    if genotype == "0" {
        return (ref_base.to_string(), String::new(), "reference".to_string());
    }

	// Homozygous non-reference
    if genotype != "0" && gt_bases.len() == 1 {
        let get_base = gt_bases[0];

        // Homozygous SNP
        if ref_base.len() == get_base.len() {

            // Simple Homozygous SNP
            if get_base.len() == 1 {
                return (get_base.to_string(), String::new(), "snp".to_string());
            }

            // SNP(s) disguised as an indel
            let mut snp_count = 0;
            let mut ref_base_saved = String::new();
            let mut cons_base_saved = String::new();
            for(index, consensus_base_char) in get_base.chars().enumerate() {
                let ref_base_char = ref_base.chars().nth(index).unwrap();
                if ref_base_char != consensus_base_char {
                    ref_base_saved = ref_base_char.to_string();
					cons_base_saved = consensus_base_char.to_string();
					snp_count+=1;
                }
            }

            return match snp_count{
                0 => (
                    ref_base.to_string(),
                    String::new(),
                    "reference".to_string(),
                ),
                1 => (
                    ref_base_saved,
                    cons_base_saved,
                    "snp".to_string(),
                ),
                _ => (
                    get_base.to_string(),
                    String::new(),
                    format!("snp_multi_{}", snp_count),
                ),
            };
        }

        // Homozygous Deletion
        if ref_base.len() > get_base.len() {
	 		return (
                get_base.to_string(),
                String::new(),
                "deletion".to_string(),
            );
        }
        // Homozygous Insertion
	 	return (
            get_base.to_string(),
            String::new(),
            "insertion".to_string(),
        );
    }

	// Bi-allelic heterozygous positions & indels
    if gt_bases.len() == 2 {
        let base1 = gt_bases[0].to_string();
        let base2 = gt_bases[1].to_string();
        
        let base_type = if base1.len() < base2.len() {
            "het_insertion"
        } else if base1.len() > base2.len() {
            "het_deletion"
        } else {
            "heterozygous"
        };

        return (base1, base2, base_type.to_string());
    }

    // polyploid?
    if gt_bases.len() > 2 {
        logger.lock().unwrap().error(&format!("polyploid not implemented currently: {}", genotype));
        process::exit(1);
    }

    logger.lock().unwrap().error(&format!("Unexpected genotype: {}", genotype));
    process::exit(1);
}

pub fn count_variants(entries: &[VCFEntry], logger: &Logger) {
    let mut sample_type_counts: HashMap<String, HashMap<String, usize>> = HashMap::new();

    for entry in entries {
        for (sample, base_type) in &entry.samples_to_base_type {
            // Skip non-variants or ambiguous calls
            if base_type == "reference" || base_type == "ambiguous" {
                continue;
            }

            let per_sample = sample_type_counts.entry(sample.clone()).or_default();
            *per_sample.entry(base_type.clone()).or_insert(0) += 1;
        }
    }

    for (sample, buckets) in &sample_type_counts {
        let snps = buckets.get("snp").copied().unwrap_or(0) + buckets.get("heterozygous").copied().unwrap_or(0);
        let ins = buckets.get("insertion").copied().unwrap_or(0) + buckets.get("het_insertion").copied().unwrap_or(0);
        let del = buckets.get("deletion").copied().unwrap_or(0) + buckets.get("het_deletion").copied().unwrap_or(0);

        logger.information(&format!(
            "Sample '{}' — SNPs: {}, insertions: {}, deletions: {}",
            sample, snps, ins, del
        ));
    }
}