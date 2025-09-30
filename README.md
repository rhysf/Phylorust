<p align="center">
  <img src="images/logo.png" alt="Phylorust" width="250"/>
</p>

**Phylorust** is a Rust-based command-line tool to generate phylogenetically informative SNP site sets (as FASTA files) and associated Tree's from VCF file(s) and FASTA reference genome files.  

---

## Features
- Reads a reference FASTA and single sample or multi-sample VCF(s).
- Generates phylogenetically informative SNP site sets at configurable coverage thresholds.
- Produces per-sample FASTA alignments.
- Runs FastTree automatically (if installed) to generate trees.
- ASCII tree rendering directly in the terminal.
- Simple tab-delimited input file (`Name_Type_Location.tab`) for managing multiple samples.

---

## Prerequisites
To build and run Phylorust, you’ll need:

- [Rust](https://www.rust-lang.org/tools/install) (stable, installed via `rustup`).  
  - Verify install with:  
    ```bash
    rustc --version
    cargo --version
    ```
- [R](https://www.r-project.org/) (for plotting histograms).  
  - Packages: `ggplot2`, `readr`, `dplyr` (install inside R with):  
    ```R
    install.packages(c("ggplot2", "readr", "dplyr"))
    ```
- [FastTree](http://www.microbesonline.org/fasttree/) (optional, for tree generation).  
  - Must be in your system `PATH`.  
  - Verify with:  
    ```bash
    FastTree -help
    ```

---

## Installation
Clone the repo and install with Cargo:

```bash
git clone https://github.com/<your-username>/Phylorust.git
cd Phylorust
cargo install --path .
```

If you installed Rust with rustup, ~/.cargo/bin is normally already in your $PATH.
If not, you can add it:

```bash
echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

You can now run:

```bash
phylorust --help
```

Alternative (manual install): If you prefer to place the binary in ~/.local/bin

```bash
cargo build --release
cp target/release/phylorust ~/.local/bin/
```

## Input file formats

Phylorust uses a tab-delimited file (Name_Type_Location.tab) with three columns:

```SampleName VCF /path/to/sample.vcf```

- SampleName = Your preferred sample label (used in output FASTAs/trees).
- Filetype (must be VCF or vcf).
- /path/to/sample.vcf = Path to the sample’s VCF file.

## Example pipeline

```bash
phylorust \
  --fasta ./examples/Cryp_gatt_R265.genome.fa-scaffold3.14.fasta \
  --name_type_location ./examples/Name_Type_location.tab
```

This will:
  1.  Parse the reference FASTA.
  2.  Parse VCFs listed in Name_Type_location.tab.
  3.  Generate SNP site sets and coverage histograms.
  4.  Produce FASTA alignments for each coverage threshold.
  5.  Run FastTree (if available) and print ASCII trees in the terminal.

## Building with Docker

```bash
git clone https://github.com/rhysf/Phylorust.git
cd Phylorust
docker build -t phylorust .

docker run --rm -v $(pwd)/examples:/examples phylorust \
  --fasta /examples/Cryp_gatt_R265.genome.fa-scaffold3.14.fasta \
  --name_type_location /examples/Name_Type_location.tab
```

## On HPC systems without Docker, you can convert the Docker image into a Singularity (Apptainer) image

```bash
apptainer build phylorust.sif docker-daemon://phylorust:latest

apptainer run phylorust.sif \
  --fasta examples/Cryp_gatt_R265.genome.fa-scaffold3.14.fasta \
  --name_type_location examples/Name_Type_location.tab
```

## Command-line arguments

Key options (full list available with --help):
  --fasta <FILE> → Reference FASTA file.
  --name_type_location <FILE> → Tab-delimited file of sample names, file type, and VCF paths.
  --output_dir <DIR> → Directory for results (default: Phylorust_output).
  --generate_fastas <MODE> → FASTA generation mode (all or specific thresholds).
  --skip-fasttree → Skip tree generation.
  --fasttree-bin <PATH> → Path to FastTree binary (if not in PATH).

## Plotting

Histograms are generated with R and saved to both .png and .pdf.
You can also run the plotting script directly:

```bash
Rscript plot_histogram.R site_coverage_histogram.tsv 90
```

## License

This project is licensed under the MIT License.
