use clap::Parser;
use std::fs;
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

fn read_VCF(name_type_location : &NameTypeLocation, logger : &Logger) {
    logger.information(&format!("read VCF: {}", name_type_location.location));

    // read file
    let vcf_file = fs::read_to_string(&name_type_location.location).unwrap_or_else(|error|{
        logger.error(&format!("Error with file: {} {}", name_type_location.location, error));
        process::exit(1);
    });

    // go through each line
    for line in vcf_file.lines() {

        // Header
        if line.starts_with("#") { continue; }

        // Entries
        let line_parts : Vec<&str> = line.split('\t').collect();

        // Initial quality check
        if line_parts.len() < 9 {
            logger.error(&format!("Error with vcf line: {}", line));
            process::exit(1);
        }

        let contig = line_parts[0];
        let position = line_parts[1];
        let id = line_parts[2];
        let ref_base = line_parts[3];
        let alt_base = line_parts[4];
        let cons_qual = line_parts[5];
        let filter = line_parts[6];
        let info = line_parts[7];
        let format = line_parts[8];
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
        read_VCF(name_type_location, &logger);
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