use crate::logger::Logger;
use crate::read_tab::NameTypeLocation;
use std::fs::{self};
use std::process;
use std::collections::{HashMap, HashSet};

pub struct VCFsamples {
    header:String,
    samples:HashMap<usize, String>,
}

#[derive(Debug, Clone)]
pub struct VCFEntry {
    pub contig:String,
    pub position:usize,
    pub id:String,
    pub ref_base:String,
    pub alt_base:String,
    pub cons_qual:String,
    pub filter:String,
    pub info:String,
    pub format:String,
    pub samples:HashMap<usize, String>,
    pub samples_to_genotype:HashMap<String, String>,
    pub samples_to_base1:HashMap<String, String>,
    pub samples_to_base2:HashMap<String, String>,
    pub samples_to_base_type:HashMap<String, String>,
}

pub fn read_vcf(name_type_location : &NameTypeLocation, logger : &Logger, global_sample_names: &mut HashSet<String>) -> Vec<VCFEntry> {
    logger.information(&format!("read VCF: {}", name_type_location.location));

    // read file
    let vcf_file = fs::read_to_string(&name_type_location.location).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", name_type_location.location, error));
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
                logger.error("VCF header line does not contain sample columns.");
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

                    logger.warning(&format!(
                        "Duplicate sample name detected: Sample '{}' renamed to '{}'.",
                        original, sample_name
                    ));
                }

                global_sample_names.insert(sample_name.clone());
                sample_names.insert(index, sample_name);
            }

            vcf_samples = Some(VCFsamples{
                header:line.to_string(),
                samples:sample_names,
            });
            continue; 
        }

        // Entries
        let vcf_samples_ref = match vcf_samples.as_ref() {
            Some(samples) => samples,
            None => {
                logger.error("read_vcf: VCF header line missing.");
                process::exit(1);
            }
        };
        
        let vcfentry = read_vcf_line(line, logger, vcf_samples_ref);
        vcf_entries.push(vcfentry);
    }
    vcf_entries
}

fn read_vcf_line(line : &str, logger : &Logger, samples: &VCFsamples) -> VCFEntry {

    let line_parts : Vec<&str> = line.split('\t').collect();
    let contig = line_parts[0].to_string();
    let position: usize = line_parts[1].parse::<usize>().unwrap_or_else(|_| {
        logger.warning(&format!("Invalid position value in VCF: '{}'", line_parts[1]));
        process::exit(1);
    });
    //let position = line_parts[1];
    let id = line_parts[2].to_string();
    let ref_base = line_parts[3].to_string();
    let alt_base = line_parts[4].to_string();
    let cons_qual = line_parts[5].to_string();
    let filter = line_parts[6].to_string();
    let info = line_parts[7].to_string();
    let format = line_parts[8].to_string();

    // Initial quality check
    if line_parts.len() < 9 {
        logger.error(&format!("Error with vcf line: {}", line));
        process::exit(1);
    }

    // get sample info
    let sample_names = samples.samples.clone();

    // get genotype
    let gt_position = read_genotype_position(line_parts[8], logger);
    let samples_to_genotype = create_sample_to_genotype(&line_parts[9..], gt_position, &samples.samples);

    // Determine what base it is
    let (samples_to_base1, samples_to_base2, samples_to_base_type) = determine_bases_and_base_type(&samples_to_genotype, &ref_base, &alt_base, &logger);

    // return various info for vcf entry
    VCFEntry{
        contig,
        position,
        id,
        ref_base,
        alt_base,
        cons_qual,
        filter,
        info,
        format,
        samples:sample_names,
        samples_to_genotype,
        samples_to_base1,
        samples_to_base2,
        samples_to_base_type,
    }
}

fn read_genotype_position(format : &str, logger : &Logger) -> usize {
    if format.starts_with("GT") {
        return 0;
    }
    
    let format_string = format.to_string();
    let format_parts : Vec<&str> = format_string.split(':').collect();

    for(index, format_part) in format_parts.iter().enumerate() {
        if format_part.to_string().to_uppercase() == "GT" {
            return index;
        }
    }

    logger.error(&format!("Error with vcf line (no genotype found): {}", format));
    process::exit(1);
}

fn create_sample_to_genotype(samples : &[&str], gt_pos : usize, sample_names : &HashMap<usize, String>) -> HashMap<String, String> {

    let mut sample_to_genotype = HashMap::new();

    for(index, sample) in samples.iter().enumerate() {
        let sample_parts : Vec<&str> = sample.split(':').collect();
        let sample_name = &sample_names[&index];
        let genotype = match sample_parts[gt_pos] {
            "0/0" => "0",
            "1/1" => "1",
            "0|0" => "0",
            "1|1" => "1",
            _ => sample_parts[gt_pos]
        };

        sample_to_genotype.insert(sample_name.to_string(), genotype.to_string());
    }

    sample_to_genotype
}

fn determine_bases_and_base_type(samples_to_genotype: &HashMap<String, String>, ref_base : &String, alt_base : &String, logger : &Logger) -> (HashMap<String, String>, HashMap<String, String>, HashMap<String, String>) {
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

fn determine_bases_and_base_type_single_sample(genotype: &String, ref_base : &String, alt_base : &String, logger : &Logger) -> (String, String, String) {

    let mut base1 = "N".to_string();
    let mut base2 = "None".to_string();
    let mut base_type = "ambiguous".to_string();

    // Save bases and GT parts
    let consensus_bases : Vec<&str> = alt_base.split(',').collect();
    let mut all_bases = vec![ref_base.as_str()];
    all_bases.extend(consensus_bases.iter().cloned());

    //logger.warning("determine_bases_and_base_type_single_sample");

    let gt_parts : Vec<&str> = genotype.split(|c| c == '/' || c == '|').collect();
    let mut gt_bases : Vec<&str> = Vec::new();
    for gt_part in gt_parts.iter() {
        let mut index = 0;

        // Ambiguous GT (all things non numerical such as . or *)
        match gt_part.parse::<usize>() {
            Ok(n) => {
                index = n;
            },
            Err(_e) => {
                base1 = "N".to_string();
                base_type = "ambiguous".to_string();
                return (base1, base2, base_type);
            }
        }

        // No base defined
        if(all_bases.len()-1) < index {
            logger.error(&format!("vcf parsing genotype error: No {} defined from ref {} and cons {}", index, ref_base, alt_base));
            process::exit(1);
        }

        let gt_base = all_bases[index];

        // Ambiguous bases
        if(gt_base == "N") || (gt_base == ".") || (gt_base == "*") {
            base1 = "N".to_string();
            base_type = "ambiguous".to_string();
            return (base1, base2, base_type);
        }

        // save the bases
        gt_bases.push(gt_base);
    }

	// Homozygous ref-calls
    if genotype == "0" {
        base1 = ref_base.to_string();
	    base_type = "reference".to_string();
	 	return (base1, base2, base_type);
    }

	// Homozygous non-reference
    if genotype != "0" && gt_bases.len() == 1 {
        let get_base = gt_bases[0];

        // Homozygous SNP
        if ref_base.len() == get_base.len() {

            // Simple Homozygous SNP
            if get_base.len() == 1 {
                base1 = get_base.to_string();
	 			base_type = "snp".to_string();
	 			return (base1, base2, base_type);
            }

            // SNP(s) disguised as an indel
            let mut snp_count = 0;
            let mut ref_base_saved = "".to_string();
            let mut cons_base_saved = "".to_string();
            for(index, consensus_base_char) in get_base.chars().enumerate() {
                let ref_base_char = ref_base.chars().nth(index).unwrap();
                if ref_base_char != consensus_base_char {
                    ref_base_saved = ref_base_char.to_string();
					cons_base_saved = consensus_base_char.to_string();
					snp_count+=1;
                }
            }

            return match snp_count{
                0 => {
                    base1 = ref_base.to_string();
	 			    base_type = "reference".to_string();
	 			    (base1, base2, base_type)
                }
                1 => {
                    base1 = ref_base_saved;
				    base2 = cons_base_saved;
				    base_type = "snp".to_string();
				    (base1, base2, base_type)
                }
                _ => {
                    base1 = get_base.to_string();
	 			    base_type = format!("snp_multi_{}", snp_count);
	 			    (base1, base2, base_type)
                }
            };
        }

        // Homozygous indel
        base1 = get_base.to_string();
        // Deletion
        if ref_base.len() > get_base.len() {
	 		base_type = "deletion".to_string();
	 		return (base1, base2, base_type);
        }
        // Insertion
	 	base_type = "insertion".to_string();
	 	return (base1, base2, base_type);
    }

	// Bi-allelic heterozygous positions & indels
    if gt_bases.len() == 2 {
        let base1 = gt_bases[0].to_string();
        let base2 = gt_bases[1].to_string();
        base_type = "heterozygous".to_string();
        if base1.len() < base2.len() { base_type = "het_insertion".to_string(); }
        if base1.len() > base2.len() { base_type = "het_deletion".to_string(); }
        return (base1, base2, base_type);
    }

    // polyploid?
    if gt_bases.len() > 2 {
        logger.error(&format!("polyploid not implemented currently: {}", genotype));
        process::exit(1);
    }

    logger.error(&format!("Unexpected genotype: {}", genotype));
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