use clap::Parser;
use std::collections::HashSet;
use std::fs::{self, exists};
use std::ops::Index;
use std::process;
use std::path::Path;
use std::collections::HashMap;
mod logger;
use logger::Logger;

// setting up the command line parameters
#[derive(Parser)]
#[command(name = "RustaTools")]
#[command(about = "A simple command-line tool fopr calculating Entirely Covered in All sites", long_about = None)]

struct Args {
    /// reference FASTA
    #[arg(short='f', long="fasta")]
    fasta_filename: String, 

    /// tab delimited file containing isolate name, type (VCF) and VCF location
    #[arg(short='n', long="name_type_location")]
    name_type_location_filename: String,

    /// Output directory 
    #[arg(short='o', long="output_dir", default_value="Rustatools_output")]
    output_dir: String,

    /// Settings 
    #[arg(short='s', long="settings", default_value_t=1)]
    settings: u8,

    // Minimum read depth
    #[arg(short='m', long="min_read_depth", default_value_t=4)]
    min_read_depth: u8,

    /// Exclude variants on contig 
    #[arg(short='e', long="exclude_contig", default_value="n")]
    exclude_contig: String,

    /// Exclude variants not on contig 
    #[arg(short='z', long="restrict_contig", default_value="n")]
    restrict_contig: String,
}

struct Fasta {
    id:String,
    desc:String,
    seq:String,
}

struct NameTypeLocation {
    name:String,
    filetype:String,
    location:String,
}

fn read_name_type_location_file(file : &str, logger : &Logger) -> Vec<NameTypeLocation> {
    logger.information(&format!("Processing file: {}", file));

    let mut name_type_locations = Vec::new();

    let content = fs::read_to_string(file).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    for(index, line) in content.lines().enumerate() {
        let line_parts:Vec<&str> = line.split('\t').collect();

        // check there are 3 columns
        if line_parts.len() != 3 {
            logger.error(&format!("Error with format of file: {} on line number {} = {}", file, index, line));
            process::exit(1);
        }

        // check files exist
        if !Path::new(line_parts[2]).exists() {
            logger.error(&format!("Error: File {} does not exist", line_parts[2]));
            process::exit(1);
        }

        // check that 2nd column says vcf in lc or uc

        // save into hash map
        name_type_locations.push(NameTypeLocation{
            name:line_parts[0].to_string(),
            filetype:line_parts[1].to_string(),
            location:line_parts[2].to_string(),
        });
    }
    return name_type_locations;
}

fn read_fasta(file : String, logger : &Logger) -> Vec<Fasta> {
    println!("Processing file: {}", file);

    // read file
    let fasta_file = fs::read_to_string(&file).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", file, error));
        process::exit(1);
    });

    // separate out the columns
    let mut last_id = "";
    let mut last_desc = "";
    let mut last_sequence = String::from("");
    let mut fasta: Vec<Fasta> = Vec::new();
    for line in fasta_file.lines() {

        // ID and Description
        if line.starts_with(">") {
            if last_id != "" {
                fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence });
            }
            last_sequence = String::from("");

            match line.find(" ") {
                Some(index) => {
                    last_id = &line[1..index];
                    last_desc = &line[index+1..];
                },
                None => { 
                    last_id = &line[1..];
                    last_desc = "";
                }
            };
            //println!("id and desc: {} {}", last_id, last_desc);
        }
        else {
            last_sequence.push_str(line);
        }
    }
    fasta.push(Fasta { id: last_id.to_string(), desc: last_desc.to_string(), seq: last_sequence }); 
    logger.information("Finished processing file");

    return fasta;
}

fn make_hashmap_of_arrays_for_genome(fasta : &Vec<Fasta>, logger : &Logger) -> HashMap<String, Vec<i32>> {
    logger.information("make_hashmap_of_arrays_for_genome: filling genome array...");

    let mut genome = HashMap::new();

    // go through fasta and make contig array of all zeros
    for entry in fasta {
        let contig_array = vec![0;entry.seq.len()];
        genome.insert(entry.id.to_string(), contig_array);
    }

    return genome;
}

//sub fill_genome_hash_array_from_vcf {
//	my ($genome_hash, $vcf, $settings, $exclude, $include, $min_depth) = @_;
//
//	# sometimes 2 vcf lines for the same position are given (Pilon). Only consider the first one
//	my $last_contig_and_position;
//	my %variants_found;
//
//	my ($ref_count, $variant_count, $min_depth_filtered, $other_filtered) = (0, 0, 0, 0);
//	warn "fill_genome_hash_array_from_vcf: Saving reference positions (1) and variants (2) over genome array from $vcf (setting=$settings, exclude=$exclude, include=$include, min depth=$min_depth)...\n";
//	open my $fh, '<', $vcf or die "Cannot open $vcf: $!\n";
//	VCF: while(my $line=<$fh>) {
//		chomp $line;
//
//		my $VCF_line = vcflines::read_VCF_lines($line);
//		my $supercontig = $$VCF_line{'supercontig'};
//		my $position = $$VCF_line{'position'};
//		my $base_type = $$VCF_line{'base_type0'};
//
//		# Ignore ambigious, or contigs of uninterest
//		next VCF if($$VCF_line{'next'} eq 1);
//		next VCF if(($exclude !~ m/n/i) && ($supercontig eq $exclude));
//		next VCF if(($include !~ m/n/i) && ($supercontig ne $include));
//
//		# Ignore lines given twice
//		if(!defined $last_contig_and_position) { $last_contig_and_position = "$supercontig\t$position"; }
//		else {
//			my $current_contig_and_position = "$supercontig\t$position";
//			if($current_contig_and_position eq $last_contig_and_position) {
//				$last_contig_and_position = $current_contig_and_position;
//				next VCF;
//			}
//			$last_contig_and_position = $current_contig_and_position;
//		}
//
//		# Min depth filtering
//		if((defined $$VCF_line{'DP0'}) && ($$VCF_line{'DP0'} ne '?'))  {
//			if($$VCF_line{'DP0'} < $min_depth) {
//				$min_depth_filtered++;
//				next VCF;
//			}
//		}
//
//		# heterozygous sites
//		if(($settings eq 1) && ($base_type eq 'heterozygous')) {
//			$other_filtered++;
//			next VCF;
//		}
//		# indels (snp_multi = E.g. ref: ACTCGTCCTGACT consensus: AACTCGTACTGAC)
//		if(($settings =~ m/[12]/) && ($base_type =~ m/insertion|deletion|snp_multi/)) {
//			$other_filtered++;
//			next VCF;
//		}
//
//		# Everything i do want to save
//		if($base_type !~ m/heterozygous|snp|insertion|deletion|reference/) {
//			$other_filtered++;
//			next VCF;
//		}
//
//		# save reference bases
//		if($base_type =~ m/reference/) {
//			$ref_count++;
//			$$genome_hash{$supercontig}[$position] = 1;
//		} else {
//			$variant_count++;
//			$variants_found{$base_type}++;
//			$$genome_hash{$supercontig}[$position] = 2;
//		}
//	}
//
//	# summary
//	warn "Reference bases:\t$ref_count\n";
//	warn "Variant bases:\t$variant_count\n";
//	foreach my $type(keys %variants_found) {
//		my $count = $variants_found{$type};
//		warn "Variant bases (type=$type)\t$count\n";
//	}
//	warn "Excluded for < min depth:\t$min_depth_filtered\n";
//	warn "Excluded for base type:\t$other_filtered\n";
//	return $genome_hash;
//}

struct VCFsamples {
    header:String,
    samples:HashMap<usize, String>,
}

struct VCFEntry {
    contig:String,
    position:String,
    id:String,
    ref_base:String,
    alt_base:String,
    cons_qual:String,
    filter:String,
    info:String,
    format:String,
    samples:HashMap<usize, String>,
    samples_to_genotype:HashMap<String, String>,
    samples_to_base1:HashMap<String, String>,
    samples_to_base2:HashMap<String, String>,
    samples_to_base_type:HashMap<String, String>,
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
            Err(e) => {
                base1 = "N".to_string();
                base_type = "ambiguous".to_string();
                return (base1, base2, base_type);
            }
        }

        // No base defined
        if((all_bases.len()-1) < index) {
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

	// // polyploid?
	// die "VCF_struct_determine_bases_and_base_type: polyploid not implemented yet: $GT\n" if(scalar(@GT_parts) > 2);

	// // Homozygous ref-calls
	// if($GT eq 0) {
	// 	$base1 = $ref_base;
	// 	$base_type = 'reference';
	// 	return ($base1, $base2, $base_type);
	// }

	// # Homozygous non-reference
	// elsif(($GT ne 0) && (scalar(@GT_bases) eq 1)) {
	// 	my $consensus = $GT_bases[0];

	// 	# Homozygous SNP
	// 	if(length($ref_base) eq length($consensus)) { 

	// 		# SNP
	// 		if((length($ref_base) eq 1) && (length($consensus) eq 1)) { 
	// 			$base1 = $consensus;
	// 			$base_type = 'snp';
	// 			return ($base1, $base2, $base_type);
	// 		}

	// 		// SNP(s) disguised as an indel
	// 		if(length($ref_base) eq length($consensus)) {
	// 			my @bases_reference = split //, $ref_base;
	// 			my @bases_consensus = split //, $consensus;
	// 			my $snp_count = 0;
	// 			my ($ref_base_saved, $cons_base_saved);
	// 			for(my $i=0; $i<scalar(@bases_reference); $i++) {
	// 				my $ref_base = $bases_reference[$i];
	// 				my $cons_base = $bases_consensus[$i];
	// 				if($ref_base ne $cons_base) { 
	// 					$ref_base_saved = $ref_base;
	// 					$cons_base_saved = $cons_base;
	// 					$snp_count++; 
	// 				}
	// 			}
	// 			if($snp_count eq 0) {
	// 				$base1 = $ref_base;
	// 				$base_type = 'reference';
	// 				return ($base1, $base2, $base_type);
	// 			}
	// 			elsif($snp_count eq 1) {
	// 				$base1 = $ref_base_saved;
	// 				$base2 = $cons_base_saved;
	// 				$base_type = 'snp';
	// 				return ($base1, $base2, $base_type);
	// 			}
	// 			if($snp_count > 1) {
	// 				$base1 = $consensus;
	// 				$base_type = ('snp_multi' . $snp_count);
	// 				return ($base1, $base2, $base_type);
	// 			}
	// 		}
	
	// 		# Ambiguous?
	// 		warn "Nothing found for this apparant homozygous snp:\n";
	// 		warn Dumper($VCF_struct);
	// 		$base1 = 'N';
	// 		$base_type = 'ambiguous';
	// 		return ($base1, $base2, $base_type);
	// 	}

	// 	# Homozygous indel
	// 	if(length($ref_base) ne length($consensus)) {

	// 		# Deletion (maybe with snps in there too!)
	// 		if(length($ref_base) > length($consensus)) { 
	// 			$base1 = $consensus;
	// 			$base_type = 'deletion';
	// 			return ($base1, $base2, $base_type);
	// 		}

	// 		# Insertion (maybe with snps in there too!)
	// 		if(length($ref_base) < length($consensus)) { 
	// 			$base1 = $consensus;
	// 			$base_type = 'insertion';
	// 			return ($base1, $base2, $base_type);
	// 		}

	// 		# Ambiguous?
	// 		warn "Nothing found for this apparent homozygous indel:\n";
	// 		warn Dumper($VCF_struct);
	// 		$base1 = 'N';
	// 		$base_type = 'ambiguous';
	// 		return ($base1, $base2, $base_type);
	// 	}
	// }

	// # Bi-allelic heterozygous positions & indels
	// elsif(scalar(@GT_bases) eq 2) {
	// 	$base1 = $GT_bases[0];
	// 	$base2 = $GT_bases[1];
	// 	$base_type = 'heterozygous';
	// 	if(length($base1) < length($base2)) { $base_type = 'het_insertion'; }
	// 	if(length($base1) > length($base2)) { $base_type = 'het_deletion'; }
	// 	return ($base1, $base2, $base_type);
	// }

	// # Ambiguous?
	// else {
	// 	die "VCF_struct_determine_bases_and_base_type: Error. Undefined variant: ref $ref_base and cons $cons_base and GT ($GT) = @GT_bases\n";
	// }

    (base1, base2, base_type)
}



fn read_vcf_line(line : &str, logger : &Logger) -> VCFEntry {

    let line_parts : Vec<&str> = line.split('\t').collect();
    let contig = line_parts[0].to_string();
    let position = line_parts[1].to_string();
    let id = line_parts[2].to_string();
    let ref_base = line_parts[3].to_string();
    let alt_base = line_parts[4].to_string();

    // Initial quality check
    if line_parts.len() < 9 {
        logger.error(&format!("Error with vcf line: {}", line));
        process::exit(1);
    }

    // get sample info
    let mut sample_names = HashMap::new();
    for(index, sample_name) in line_parts[9..].iter().enumerate() {
        sample_names.insert(index, sample_name.to_string());
    }

    // get genotype
    let gt_position = read_genotype_position(line_parts[8], logger);
    let samples_to_genotype = create_sample_to_genotype(&line_parts[9..], gt_position, &sample_names);

    // Determine what base it is
    let (samples_to_base1, samples_to_base2, samples_to_base_type) = determine_bases_and_base_type(&samples_to_genotype, &ref_base, &alt_base, &logger);

    // return various info for vcf entry
    VCFEntry{
        contig,
        position,
        id,
        ref_base,
        alt_base,
        cons_qual:line_parts[5].to_string(),
        filter:line_parts[6].to_string(),
        info:line_parts[7].to_string(),
        format:line_parts[8].to_string(),
        samples:sample_names,
        samples_to_genotype,
        samples_to_base1,
        samples_to_base2,
        samples_to_base_type,
    }
}

fn read_vcf(name_type_location : &NameTypeLocation, logger : &Logger) {
    logger.information(&format!("read VCF: {}", name_type_location.location));

    // read file
    let vcf_file = fs::read_to_string(&name_type_location.location).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", name_type_location.location, error));
        process::exit(1);
    });

    let mut vcf_samples:Option<VCFsamples> = None;

    // go through each line of the VCF
    for line in vcf_file.lines() {

        // Metadata
        if line.starts_with("##") { 
            continue; 
        }
        // Header
        if line.starts_with("#CHROM\tPOS\tID\tREF") { 

            let line_parts : Vec<&str> = line.split('\t').collect();
            let mut sample_names = HashMap::new();
            for(index, sample_name) in line_parts[9..].iter().enumerate() {
                logger.warning(&format!("index = {}", index));
                sample_names.insert(index, sample_name.to_string());
            }

            vcf_samples = Some(VCFsamples{
                header:line.to_string(),
                samples:sample_names,
            });

            logger.warning(line);

            continue; 
        }

        // Entries
        let vcfentry = read_vcf_line(line, logger);

    }
}

fn main() {

    let args = Args::parse();
    let logger = Logger;

    // Read Name Type Location file
    let name_type_locations = read_name_type_location_file(&args.name_type_location_filename, &logger);

    // Info about VCF's
    let number_of_vcfs = name_type_locations.len();
    println!("{} VCF files specified", number_of_vcfs);

    // Read FASTA file to memory
    let fasta = read_fasta(args.fasta_filename, &logger);

    // go through entries
    for entry in &fasta {
        logger.information(&format!("id and desc: {} {}", entry.id, entry.desc));
        logger.information(&format!("length of seq: {}", entry.seq.len()));
    }

    // Make output folder
    fs::create_dir_all(&args.output_dir).unwrap_or_else(|error|{
        logger.error(&format!("Error with output directory: {}", error));
        process::exit(1);
    });

    // Output files
    let name_type_location_filename_path = Path::new(&args.name_type_location_filename);
    let name_type_location_filename_wo_dir = name_type_location_filename_path.file_name().unwrap().to_str();
    logger.output(name_type_location_filename_wo_dir.unwrap());

    let rustatools_settings = &format!("m-{}-s-{}-e-{}-z-{}", 
        args.min_read_depth, 
        args.settings, 
        args.exclude_contig, 
        args.restrict_contig);
    
    let outfile_reference_bases = &format!("{}/{}-{}-reference-bases.tab", 
        args.output_dir,
        name_type_location_filename_wo_dir.unwrap(),
        rustatools_settings);

    let outfile_variant_bases = &format!("{}/{}-{}-variant-bases.tab", 
        args.output_dir,
        name_type_location_filename_wo_dir.unwrap(),
        rustatools_settings);

    logger.output(&format!("Outfile reference bases = {}", outfile_reference_bases));
    logger.output(&format!("Outfile variant bases = {}", outfile_variant_bases));

    // Save genome to hashmap of arrays
    let mut genome = make_hashmap_of_arrays_for_genome(&fasta, &logger);

    // go through each VCF
    for name_type_location in &name_type_locations {
        read_vcf(name_type_location, &logger);
    }

    //# Convert ref bases to 1 and variants to 2
    // Start by doing this for a single VCF, and then put in a loop to do it for every VCF
    //$genome_array = genomearray::fill_genome_hash_array_from_vcf($genome_array, $opt_v, $opt_s, $opt_e, $opt_z, $opt_m);

    //# Print tab files for locations of 1s (reference)
    //genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 1, $outfile1);

    //# Print tab files for locations of 1s (variant)
    //genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 2, $outfile2);

    logger.output("output test");
    logger.warning("output test");
}