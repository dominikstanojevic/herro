// A corrupted gzip reads file must stop HERRO with exit code 1, the "corrupted" message and no
// output, on the rapidgzip path; legacy must fail too (it panics). Runs the built binary.
use std::{
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};

use flate2::{write::GzEncoder, Compression};

// ~2.5 MB of FASTQ, ~1 MB gzipped.
fn fastq() -> Vec<u8> {
    let mut x: u64 = 0x2545_F491_4F6C_DD1D;
    let mut out = Vec::new();
    for i in 0..1000 {
        write!(out, "@read{i} d\n").unwrap();
        let mut seq = Vec::new();
        for _ in 0..1200 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            seq.push(b"ACGT"[(x % 4) as usize]);
        }
        out.extend_from_slice(&seq);
        out.extend_from_slice(b"\n+\n");
        out.extend(seq.iter().map(|b| b'!' + b % 40));
        out.push(b'\n');
    }
    out
}

fn run(reads: &Path, out: &Path, decoder: &str) -> (Option<i32>, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_herro"))
        .args(["features", "--gzip-decoder", decoder])
        .arg(reads)
        .arg(out)
        .output()
        .unwrap();
    (o.status.code(), String::from_utf8_lossy(&o.stderr).into_owned())
}

#[test]
fn corrupted_gzip_exits_without_output() {
    let dir: PathBuf = std::env::temp_dir().join(format!("herro-c25-it-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    let mut e = GzEncoder::new(Vec::new(), Compression::new(6));
    e.write_all(&fastq()).unwrap();
    let good = e.finish().unwrap();
    let n = good.len();

    let mut cases: Vec<(&str, Vec<u8>)> = Vec::new();
    cases.push(("truncated", good[..n * 2 / 3].to_vec()));
    let mut flipped = good.clone();
    flipped[n / 2] ^= 0xff;
    cases.push(("flipped", flipped));
    let mut bad_crc = good.clone();
    bad_crc[n - 8] ^= 1;
    cases.push(("bad_crc", bad_crc));
    let mut bad_isize = good.clone();
    bad_isize[n - 4] ^= 1;
    cases.push(("bad_isize", bad_isize));
    cases.push(("trailing_garbage", [&good[..], b"trailing garbage\n"].concat()));

    for (name, data) in cases {
        let reads = dir.join(format!("{name}.fq.gz"));
        std::fs::write(&reads, data).unwrap();
        let out = dir.join(format!("{name}.out"));

        let (code, stderr) = run(&reads, &out, "rapidgzip");
        assert_eq!(code, Some(1), "{name} rapidgzip: {stderr}");
        let corrupted = format!(
            "Reads file {} is corrupted or not valid FASTQ: ",
            reads.display()
        );
        assert!(stderr.contains(&corrupted), "{name} rapidgzip: {stderr}");
        assert!(stderr.contains("\nrapidgzip-core 0.3.1 error: "), "{name} rapidgzip: {stderr}");
        assert!(!out.exists(), "{name} rapidgzip: output created");
        eprintln!("{name} rapidgzip: {}", stderr.trim_end());

        let (code, stderr) = run(&reads, &out, "legacy");
        assert_ne!(code, Some(0), "{name} legacy: {stderr}");
        assert!(stderr.contains("Error parsing fastx file."), "{name} legacy: {stderr}");
        assert!(!out.exists(), "{name} legacy: output created");
    }
    std::fs::remove_dir_all(dir).unwrap();
}
