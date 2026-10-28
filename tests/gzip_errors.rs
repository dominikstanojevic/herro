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
    let dir: PathBuf = std::env::temp_dir().join(format!("herro-test-it-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    let good = gz(&fastq());
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
        // Second line: the decoder's error, or (flipped: bad bytes before the CRC) the parse error.
        let second = match name {
            "flipped" => "\nParse error on the rapidgzip path: ",
            _ => "\nrapidgzip-core 0.3.1 error: ",
        };
        assert!(stderr.contains(second), "{name} rapidgzip: {stderr}");
        assert!(!out.exists(), "{name} rapidgzip: output created");
        eprintln!("{name} rapidgzip: {}", stderr.trim_end());

        let (code, stderr) = run(&reads, &out, "legacy");
        assert_ne!(code, Some(0), "{name} legacy: {stderr}");
        assert!(stderr.contains("Error parsing fastx file."), "{name} legacy: {stderr}");
        assert!(!out.exists(), "{name} legacy: output created");
    }
    std::fs::remove_dir_all(dir).unwrap();
}

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = GzEncoder::new(Vec::new(), Compression::new(6));
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

// bzip2 -9 of "@r1\nACGT\n+\nIIII\n".
const BZ2_FASTQ: [u8; 59] = [
    0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26, 0x53, 0x59, 0xaf, 0x85, 0x72, 0x8b, 0x00, 0x00,
    0x03, 0xde, 0x80, 0x40, 0x10, 0x00, 0x08, 0x20, 0x00, 0x68, 0xa0, 0x04, 0x00, 0x10, 0x00, 0x20,
    0x00, 0x22, 0x01, 0xa3, 0x4d, 0x08, 0x06, 0x9a, 0x68, 0x3d, 0x20, 0x05, 0x0c, 0x78, 0xbd, 0x25,
    0xe2, 0xee, 0x48, 0xa7, 0x0a, 0x12, 0x15, 0xf0, 0xae, 0x51, 0x60,
];

// Like parse_fastx_file, the rapidgzip path does not look for a second compression layer.
#[test]
fn nested_compression_is_rejected_by_both() {
    let dir: PathBuf = std::env::temp_dir().join(format!("herro-test-nest-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    for (name, data) in [("gz_gz", gz(&gz(&fastq()))), ("gz_bz2", gz(&BZ2_FASTQ))] {
        let reads = dir.join(format!("{name}.fq.gz"));
        std::fs::write(&reads, data).unwrap();
        let out = dir.join(format!("{name}.out"));

        let (code, stderr) = run(&reads, &out, "rapidgzip");
        assert_eq!(code, Some(1), "{name} rapidgzip: {stderr}");
        let corrupted = format!(
            "Reads file {} is corrupted or not valid FASTQ: Expected '@' or '>' at the start",
            reads.display()
        );
        assert!(stderr.contains(&corrupted), "{name} rapidgzip: {stderr}");
        assert!(stderr.contains("\nParse error on the rapidgzip path: "), "{name}: {stderr}");
        assert!(!out.exists(), "{name} rapidgzip: output created");
        eprintln!("{name} rapidgzip: {}", stderr.trim_end());

        let (code, stderr) = run(&reads, &out, "legacy");
        assert_ne!(code, Some(0), "{name} legacy: {stderr}");
        assert!(stderr.contains("Cannot open file containing reads."), "{name} legacy: {stderr}");
        assert!(stderr.contains("UnknownFormat"), "{name} legacy: {stderr}");
        assert!(!out.exists(), "{name} legacy: output created");
    }
    std::fs::remove_dir_all(dir).unwrap();
}
