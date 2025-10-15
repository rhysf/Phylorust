use clap::Parser;

// setting up the command line parameters
#[derive(Parser, Debug)]
#[command(name = "Phylorust")]
#[command(about = "Generate phylogenetically informative SNP site sets and FASTAs from multi-sample VCFs", long_about = None)]

pub struct Args {
    /// reference FASTA
    #[arg(short='f', long="fasta")]
    pub fasta_filename: String, 

    /// tab delimited file containing isolate name, type (VCF) and VCF location
    #[arg(short='n', long="name_type_location")]
    pub name_type_location_filename: String,

    /// Output directory 
    #[arg(short='o', long="output_dir", default_value="phylorust_output")]
    pub output_dir: String,

    /// Variant filtering mode: 
    /// 1 = homozygous SNPs only,
    /// 2 = homozygous SNPs + heterozygous positions,
    /// 3 = all variants (including indels)
    #[arg(short='s', long="settings", default_value_t=1, value_parser = clap::value_parser!(u8).range(1..=3))]
    pub settings: u8,

    /// Include invariant sites in downstream analysis (optional flag)
    /// Default = false (invariant sites excluded)
    #[arg(long = "include_invariants", default_value_t = false, action)]
    pub include_invariants: bool,

    // Minimum read depth
    #[arg(short='m', long="min_read_depth", default_value_t=4)]
    pub min_read_depth: u8,

    /// Exclude variants on contig 
    #[arg(short='e', long="exclude_contig", default_value="n")]
    pub exclude_contig: String,

    /// Exclude variants not on contig 
    #[arg(short='z', long="restrict_contig", default_value="n")]
    pub restrict_contig: String,

    /// How to represent heterozygous variants: "alt" (default) or "iupac"
    #[arg(short='h', long="heterozygosity_encoding", default_value="alt")]
    pub heterozygosity_encoding: String,

    /// Percent threshold to extract sites for tree construction
    #[arg(short='p', long="percent_for_tree", default_value_t=90)]
    pub percent_threshold: usize,

    /// Generate FASTAs for a comma-separated list of percentiles (e.g., 80,90,95) or "all"
    #[arg(long = "generate_fastas", default_value = "")]
    pub generate_fastas: String,

    /// Optional path to the FastTree binary.
    /// 
    /// If not provided, the program will search for `FastTree` in the system PATH.
    /// Example: `--fasttree-bin ./bin/FastTree`
    ///
    /// By default, nucleotide mode (`-nt`) is used.
    #[arg(long)]
    pub fasttree_bin: Option<String>,

    /// Skip FastTree tree generation, even if FastTree is available.
    #[arg(long)]
    pub skip_fasttree: bool,
}