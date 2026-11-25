use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

struct TestDirectory(PathBuf);

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn read_alignment(path: &Path) -> BTreeMap<String, String> {
    let fasta = fs::read_to_string(path).unwrap();
    let mut lines = fasta.lines();
    let mut alignment = BTreeMap::new();
    while let Some(header) = lines.next() {
        let sample = header.strip_prefix('>').unwrap();
        alignment.insert(sample.to_owned(), lines.next().unwrap().to_owned());
    }
    alignment
}

#[test]
fn counts_reference_and_variant_calls_once_per_isolate() {
    let dir = TestDirectory(
        std::env::temp_dir().join(format!("phylorust-called-sites-{}", std::process::id())),
    );
    fs::create_dir(&dir.0).unwrap();
    fs::write(dir.0.join("ref.fa"), ">chr\nccccc\n").unwrap();
    fs::write(dir.0.join("inputs.tab"), "cohort\tVCF\tcohort.vcf\n").unwrap();
    fs::write(
        dir.0.join("cohort.vcf"),
        concat!(
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\tthree\tfour\tfive\tsix\n",
            "chr\t1\t.\tc\tA\t50\tPASS\t.\tGT:DP:MIN_DP\t1/1:7:.\t1/1:7:.\t1/1:7:.\t0/0:.:7\t0/0:.:7\t0/0:.:7\n",
            "chr\t2\t.\tc\tG\t50\tPASS\t.\tGT:MIN_DP\t1/1:7\t0/0:7\t0/0:7\t0/0:7\t0/0:7\t0/0:0\n",
            "chr\t3\t.\tc\tT\t50\tPASS\t.\tGT:DP:MIN_DP\t1|0:7:.\t0/0:7:.\t0/0:7:.\t0/0:7:.\t0/0:7:.\t0/0:7:.\n",
            "chr\t4\t.\tc\tA\t50\tPASS\t.\tGT:DP:MIN_DP\t1/1:7:.\t0/0:7:.\t0/0:7:.\t0/0:7:.\t0/0:7:.\t0/0:0:7\n",
            "chr\t5\t.\tc\tA\t50\tPASS\t.\tGT:DP\t1/1:7\t0/0:7\t0/0:7\t0/0:7\t0/0:7\t./.:7\n",
        ),
    ).unwrap();

    for (mode, full, partial) in [("1", "A", "AGAA"), ("2", "AT", "AGTAA")] {
        let output = Command::new(env!("CARGO_BIN_EXE_phylorust"))
            .current_dir(&dir.0)
            .args([
                "-f",
                "ref.fa",
                "-n",
                "inputs.tab",
                "-s",
                mode,
                "-t",
                "1",
                "--generate_fastas",
                "80,90",
                "--skip-fasttree",
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );

        let output_dir = dir.0.join("phylorust_output/inputs.tab");
        let settings = format!("calls-v2-m-4-e-n-z-n-h-alt-s-{mode}");
        let full_alignment =
            read_alignment(&output_dir.join(format!("percent_90-{settings}.fasta")));
        let partial_alignment =
            read_alignment(&output_dir.join(format!("percent_80-{settings}.fasta")));
        assert_eq!(full_alignment.len(), 6);
        assert_eq!(partial_alignment.len(), 6);
        assert_eq!(full_alignment["one"], full);
        assert_eq!(partial_alignment["one"], partial);
        assert_eq!(full_alignment["six"], "C".repeat(full.len()));
        assert_eq!(
            partial_alignment["six"],
            if mode == "1" { "CNNN" } else { "CNCNN" }
        );
    }
}
