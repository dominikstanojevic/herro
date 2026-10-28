use core::panic;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, Read},
    ops::RangeBounds,
    path::Path,
    sync::{
        atomic::{AtomicBool, Ordering},
        Mutex, Once,
    },
    thread,
};

use crossbeam_channel::{bounded, unbounded};
use flate2::read::MultiGzDecoder;
use needletail::{
    errors::ParseError,
    parse_fastx_file,
    parser::{FastaReader, FastqReader},
    FastxReader,
};
use rapidgzip_core::Decoder;
use rustc_hash::FxHashSet;

// get_reads: worker threads (id split, 2-bit encode, core/neighbour filter) and records per chunk.
const PARSE_WORKERS: usize = 4;
const PARSE_CHUNK_SIZE: usize = 64;

// rapidgzip-core decoder threads (a ceiling; the library adapts below it).
const RAPIDGZIP_THREADS: usize = 16;
// Consecutive Interrupted errors retried before one is passed on as a hard error (see TrackedRead).
const MAX_INTERRUPTED_RETRIES: u32 = 10;

const BASE_ENCODING: [u64; 256] = [
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255,
];

const BASE_DECODING: [u8; 4] = [b'A', b'C', b'G', b'T'];

pub struct HAECRecord {
    pub id: Vec<u8>,
    pub description: Option<Vec<u8>>,
    pub seq: HAECSeq,
    pub qual: Vec<u8>,
}

impl HAECRecord {
    fn new(id: Vec<u8>, description: Option<Vec<u8>>, seq: HAECSeq, qual: Vec<u8>) -> Self {
        HAECRecord {
            id,
            description,
            seq,
            qual,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum GzipDecoder {
    /// needletail's own decoder (flate2, single-threaded)
    Legacy,
    /// rapidgzip-core 0.3.1 (parallel) for gzip files; other inputs as legacy
    Rapidgzip,
}

pub fn get_reads<P: AsRef<Path>>(
    path: P,
    min_length: u32,
    core: &Option<FxHashSet<String>>,
    neighbour: &Option<FxHashSet<String>>,
    gzip_decoder: GzipDecoder,
) -> Vec<HAECRecord> {
    let path = path.as_ref();
    // Known difference: rapidgzip-core accepts a dynamic block that declares
    // 31 or 32 distance codes (HDIST 30/31, allowed by RFC 1951, the extra codes unused), which
    // zlib-rs rejects. On such a file legacy aborts and rapidgzip decodes the CRC-checked bytes.
    // Test: tests::rapidgzip_accepts_hdist31.
    if gzip_decoder == GzipDecoder::Rapidgzip && is_gzip_file(path) {
        return get_reads_rapidgzip(path, min_length, core, neighbour);
    }

    let reader = parse_fastx_file(path).expect("Cannot open file containing reads.");
    read_records(reader, min_length, core, neighbour, false).unwrap()
}

// A regular file starting with the gzip magic. Anything else stays on needletail's own path.
fn is_gzip_file(path: &Path) -> bool {
    let mut magic = [0; 2];
    match File::open(path) {
        Ok(mut file) => {
            file.metadata().is_ok_and(|m| m.is_file())
                && file.read_exact(&mut magic).is_ok()
                && magic == [0x1f, 0x8b]
        }
        Err(_) => false,
    }
}

fn get_reads_rapidgzip(
    path: &Path,
    min_length: u32,
    core: &Option<FxHashSet<String>>,
    neighbour: &Option<FxHashSet<String>>,
) -> Vec<HAECRecord> {
    install_rapidgzip_panic_hook();
    let error = Mutex::new(None);

    let decoder = Decoder::builder()
        .decoder_threads(RAPIDGZIP_THREADS)
        .build()
        .expect("Valid rapidgzip-core configuration.");
    let result = match decoder.open(path) {
        Ok(inner) => {
            let tracked = TrackedRead {
                inner,
                delivered: 0,
                error: &error,
            };
            // The reader (and with it all decoding) is dropped before this returns.
            decoded_fastx_reader(tracked)
                .and_then(|reader| read_records(reader, min_length, core, neighbour, true))
        }
        Err(e) => {
            *error.lock().unwrap() = Some((io::Error::other(e.clone()), 0));
            Err(ParseError::from(io::Error::other(e)))
        }
    };

    // Nothing is kept on error: classify, then exit before any output exists.
    result.unwrap_or_else(|parse_error| {
        let decoder_error = error.into_inner().unwrap();
        eprintln!(
            "{}",
            rapidgzip_failure_message(path, decoder_error, parse_error)
        );
        std::process::exit(1);
    })
}

// Release has panic = "abort": a panic in a rapidgzip-core thread ends HERRO at once, so no
// WorkerPanicked error reaches get_reads. The hook names the culprit, then runs the previous hook
// (default message, then abort). Installed once per process.
fn install_rapidgzip_panic_hook() {
    static HOOK: Once = Once::new();
    HOOK.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if let Some(msg) = rapidgzip_panic_message(info.location().map(|l| l.file())) {
                eprintln!("{msg}");
            }
            previous(info);
        }));
    });
}

fn rapidgzip_panic_message(file: Option<&str>) -> Option<&'static str> {
    file.filter(|f| f.contains("rapidgzip-core")).map(|_| {
        "rapidgzip-core 0.3.1 crashed while decoding the reads file. This is a rapidgzip-core problem; rerun with --gzip-decoder legacy."
    })
}

// Message for a failed rapidgzip read. Re-reads the file with the legacy decoder (this thread,
// records discarded) to tell a bad file from a rapidgzip-core problem.
fn rapidgzip_failure_message(
    path: &Path,
    decoder_error: Option<(io::Error, u64)>,
    parse_error: ParseError,
) -> String {
    let rapidgzip_error = match &decoder_error {
        Some((e, _)) => e.to_string(),
        None => parse_error.to_string(),
    };

    match legacy_check(path) {
        Ok(()) => format!(
            "rapidgzip-core 0.3.1 failed on {} ({}) while the standard decoder reads the file fine. This is a rapidgzip-core problem; rerun with `--gzip-decoder legacy`.",
            path.display(),
            rapidgzip_error
        ),
        Err(legacy_error) => format!(
            "Reads file {} is corrupted or not valid FASTQ: {}\n{}: {}",
            path.display(),
            legacy_error,
            match decoder_error {
                Some(_) => "rapidgzip-core 0.3.1 error",
                None => "Parse error on the rapidgzip path",
            },
            rapidgzip_error
        ),
    }
}

// Parses the gzip file to the first error or the end with needletail on flate2's MultiGzDecoder
// (the decoder parse_fastx_file uses for gzip), keeping nothing. A decoder error carries the
// decompressed byte offset at which it occurred.
fn legacy_check(path: &Path) -> Result<(), String> {
    let error = Mutex::new(None);
    let file = File::open(path).map_err(|e| e.to_string())?;
    let tracked = TrackedRead {
        inner: MultiGzDecoder::new(file),
        delivered: 0,
        error: &error,
    };
    let result = decoded_fastx_reader(tracked).and_then(|mut reader| {
        while let Some(record) = reader.next() {
            record?;
        }
        Ok(())
    });
    result.map_err(|parse_error| match error.into_inner().unwrap() {
        Some((e, offset)) => format!("{} (at decompressed byte {})", e, offset),
        None => parse_error.to_string(),
    })
}

// needletail 0.5.1's gzip branch of parse_fastx_reader, on an already decoded stream: the format
// comes from the first decoded byte, with no sniffing for another compression layer (so nested
// gz/bz2/xz is rejected as by parse_fastx_file). No first byte: the io error, as there.
fn decoded_fastx_reader<'a, R: Read + Send + 'a>(
    mut reader: R,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    let mut first = [0; 1];
    reader.read_exact(&mut first)?;
    let reader = io::Cursor::new(first).chain(reader);
    match first[0] {
        b'>' => Ok(Box::new(FastaReader::new(reader))),
        b'@' => Ok(Box::new(FastqReader::new(reader))),
        b => Err(ParseError::new_unknown_format(b)),
    }
}

// Counts the bytes it delivers and keeps the first read error with the count at that point,
// before passing it on (needletail keeps only the message, or none at all for the
// first two bytes). rapidgzip-core 0.3.1 returns a failed decode's error on every later read, so an
// Interrupted error that repeats is final (e.g. EINTR from pread): it is passed on as a
// hard error, because needletail and read_exact would retry it forever.
struct TrackedRead<'a, R> {
    inner: R,
    delivered: u64,
    error: &'a Mutex<Option<(io::Error, u64)>>,
}

impl<R: Read> Read for TrackedRead<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut retries = 0;
        loop {
            match self.inner.read(buf) {
                Ok(n) => {
                    self.delivered += n as u64;
                    return Ok(n);
                }
                Err(e)
                    if e.kind() == io::ErrorKind::Interrupted
                        && retries < MAX_INTERRUPTED_RETRIES =>
                {
                    retries += 1
                }
                Err(e) => {
                    let kind = match e.kind() {
                        io::ErrorKind::Interrupted => io::ErrorKind::Other,
                        kind => kind,
                    };
                    let passed_on = io::Error::new(kind, e.to_string());
                    self.error.lock().unwrap().get_or_insert((e, self.delivered));
                    return Err(passed_on);
                }
            }
        }
    }
}

// Reads all records into HAECRecords, in file order. With `stop_on_error`, a parse error ends the
// read and is returned (workers skip the chunks still queued); without it, it panics in place.
fn read_records<'a>(
    mut reader: Box<dyn FastxReader + 'a>,
    min_length: u32,
    core: &Option<FxHashSet<String>>,
    neighbour: &Option<FxHashSet<String>>,
    stop_on_error: bool,
) -> Result<Vec<HAECRecord>, ParseError> {
    // This thread inflates and splits records, appending each record's id, seq and qual (as returned
    // by needletail) to a chunk buffer, with their end offsets; workers build HAECRecords from it.
    // Chunks carry a sequence number and are appended in that order, so `reads` keeps file order.
    let (chunk_sender, chunk_receiver) =
        bounded::<(usize, Vec<u8>, Vec<[usize; 3]>)>(2 * PARSE_WORKERS);
    let (recs_sender, recs_receiver) = unbounded::<(usize, Vec<HAECRecord>)>();
    // Workers return emptied chunk buffers, so this thread refills warm memory instead of allocating.
    let (pool_sender, pool_receiver) = unbounded::<(Vec<u8>, Vec<[usize; 3]>)>();
    // Set only on a parse error with `stop_on_error`: records near a decoder error may be garbage
    // and must not reach invalid_base_exit before the error is reported.
    let stopped = AtomicBool::new(false);

    thread::scope(|s| {
        for _ in 0..PARSE_WORKERS {
            let chunk_receiver = chunk_receiver.clone();
            let recs_sender = recs_sender.clone();
            let pool_sender = pool_sender.clone();
            let stopped = &stopped;

            s.spawn(move || {
                for (i, mut buf, mut ends) in chunk_receiver {
                    if stopped.load(Ordering::Relaxed) {
                        continue;
                    }
                    let mut recs = Vec::with_capacity(ends.len());
                    let mut start = 0;
                    for &[id_end, seq_end, qual_end] in &ends {
                        let line = &buf[start..id_end];
                        let seq = &buf[id_end..seq_end];
                        let qual = buf[seq_end..qual_end].to_vec();
                        start = qual_end;

                        let mut split = line.splitn(2, |c| *c == b' ' || *c == b'\t');
                        let id = split.next().expect("Invalid read id.").to_owned();
                        let description = split.next().map(|d| d.to_owned());

                        if let (Some(core), Some(neighbour)) = (&core, &neighbour) {
                            if !neighbour.contains(std::str::from_utf8(&id).unwrap())
                                && !core.contains(std::str::from_utf8(&id).unwrap())
                            {
                                continue;
                            }
                        }

                        let seq = match encode(seq) {
                            Ok((data, length)) => HAECSeq::new(data, length),
                            Err(i) => invalid_base_exit(seq[i], i, &id),
                        };

                        recs.push(HAECRecord::new(id, description, seq, qual));
                    }
                    recs_sender.send((i, recs)).unwrap();

                    buf.clear();
                    ends.clear();
                    pool_sender.send((buf, ends)).unwrap();
                }
            });
        }
        drop(recs_sender);

        let collector = s.spawn(move || {
            let mut reads = Vec::new();
            let mut pending = BTreeMap::new();
            let mut next = 0;
            for (i, recs) in recs_receiver {
                pending.insert(i, recs);
                while let Some(recs) = pending.remove(&next) {
                    reads.extend(recs);
                    next += 1;
                }
            }

            reads.shrink_to_fit();
            reads
        });

        let mut n_chunks = 0;
        let mut parse_error = None;
        let (mut buf, mut ends) = pool_receiver.try_recv().unwrap_or_default();
        while let Some(record) = reader.next() {
            let record = match record {
                Ok(record) => record,
                Err(e) if stop_on_error => {
                    stopped.store(true, Ordering::Relaxed);
                    parse_error = Some(e);
                    break;
                }
                Err(e) => panic!("Error parsing fastx file.: {e:?}"),
            };
            if record.num_bases() < min_length as usize {
                continue;
            }

            let qual = record.qual().expect("Qualities should be present.");
            buf.extend_from_slice(record.id());
            let id_end = buf.len();
            buf.extend_from_slice(&record.seq());
            let seq_end = buf.len();
            buf.extend_from_slice(qual);
            ends.push([id_end, seq_end, buf.len()]);

            if ends.len() == PARSE_CHUNK_SIZE {
                chunk_sender.send((n_chunks, buf, ends)).unwrap();
                n_chunks += 1;
                (buf, ends) = pool_receiver.try_recv().unwrap_or_default();
            }
        }
        if !ends.is_empty() && parse_error.is_none() {
            chunk_sender.send((n_chunks, buf, ends)).unwrap();
        }
        drop(chunk_sender);

        let reads = collector.join().unwrap();
        match parse_error {
            Some(e) => Err(e),
            None => Ok(reads),
        }
    })
}

#[derive(PartialEq, Debug)]
pub struct HAECSeq {
    data: Vec<u64>,
    length: usize,
}

impl HAECSeq {
    pub fn new(data: Vec<u64>, length: usize) -> Self {
        HAECSeq { data, length }
    }

    pub fn len(&self) -> usize {
        return self.length;
    }

    pub fn get_sequence(&self, buffer: &mut [u8]) {
        decode(&self.data, self.length, .., false, buffer)
    }

    pub fn get_subseq<R: RangeBounds<usize>>(&self, range: R, buffer: &mut [u8]) {
        decode(&self.data, self.length, range, false, buffer)
    }

    pub fn get_rc_subseq<R: RangeBounds<usize>>(&self, range: R, buffer: &mut [u8]) {
        decode(&self.data, self.length, range, true, buffer)
    }
}

impl From<&[u8]> for HAECSeq {
    fn from(value: &[u8]) -> Self {
        let (data, length) =
            encode(value).unwrap_or_else(|i| panic!("Invalid base at position {i}."));
        HAECSeq::new(data, length)
    }
}

impl From<&HAECSeq> for Vec<u8> {
    fn from(value: &HAECSeq) -> Self {
        let mut seq = vec![0; value.length];
        decode(&value.data, value.length, .., false, &mut seq);

        seq
    }
}

// Err holds the position of the first base that is not A, C, G or T (either case).
fn encode(sequence: &[u8]) -> Result<(Vec<u64>, usize), usize> {
    let mut data = Vec::with_capacity((sequence.len() + 31) / 32);
    let mut block = 0;
    // ORs all codes; any invalid base sets bits above the low two, checked once after the loop.
    let mut acc = 0;

    for (i, b) in sequence.iter().enumerate() {
        let c = BASE_ENCODING[*b as usize];
        acc |= c;

        block |= c << ((i << 1) & 63);
        if (i + 1) & 31 == 0 || i == sequence.len() - 1 {
            data.push(block);
            block = 0;
        }
    }

    if acc > 3 {
        return Err(sequence
            .iter()
            .position(|b| BASE_ENCODING[*b as usize] > 3)
            .unwrap());
    }

    Ok((data, sequence.len()))
}

// Called from get_reads' worker threads, before any output exists. Exits with a clear message and
// exit code 1 instead of a panic (release has panic = "abort": SIGABRT and a panic dump).
fn invalid_base_exit(base: u8, pos: usize, id: &[u8]) -> ! {
    eprintln!(
        "Invalid base '{}' at position {} (0-based) in read {}. HERRO accepts only A, C, G, T (upper or lower case).",
        base.escape_ascii(),
        pos,
        String::from_utf8_lossy(id)
    );
    std::process::exit(1);
}

fn decode<R: RangeBounds<usize>>(
    sequence: &[u64],
    length: usize,
    range: R,
    is_reversed: bool,
    buffer: &mut [u8],
) {
    let start = match range.start_bound() {
        std::ops::Bound::Unbounded => 0,
        std::ops::Bound::Included(s) => *s,
        std::ops::Bound::Excluded(s) => *s + 1,
    };

    let end = match range.end_bound() {
        std::ops::Bound::Unbounded => length,
        std::ops::Bound::Included(e) => *e + 1,
        std::ops::Bound::Excluded(e) => *e,
    };

    assert!(end <= length, "Out of bounds for 2-bit sequence decoding.");

    if start >= end {
        return;
    }

    let out = &mut buffer[..end - start];
    let mut idx = 0;
    if !is_reversed {
        // Walk words forward; base i is bits 2*(i % 32).. of word i / 32.
        let mut pos = start;
        while pos < end {
            let n = (32 - (pos & 31)).min(end - pos);
            let mut word = sequence[pos >> 5] >> ((pos & 31) << 1);
            for b in &mut out[idx..idx + n] {
                *b = BASE_DECODING[(word & 3) as usize];
                word >>= 2;
            }
            idx += n;
            pos += n;
        }
    } else {
        // Walk words backward from end - 1, emitting complemented bases.
        // The current base is kept in the top 2 bits of `word`.
        let mut pos = end;
        while pos > start {
            let last = pos - 1;
            let n = ((last & 31) + 1).min(pos - start);
            let mut word = sequence[last >> 5] << ((31 - (last & 31)) << 1);
            for b in &mut out[idx..idx + n] {
                *b = BASE_DECODING[((word >> 62) ^ 3) as usize];
                word <<= 2;
            }
            idx += n;
            pos -= n;
        }
    }
}

pub(crate) fn bytes_to_u32(bytes: &[u8]) -> u32 {
    bytes.iter().fold(0, |acc, &d| {
        if d.is_ascii_digit() {
            acc * 10 + (d - b'0') as u32
        } else {
            panic!("Character is not a valid digit");
        }
    })
}

#[cfg(test)]
mod tests {
    use std::{
        io::{self, Read, Write},
        path::{Path, PathBuf},
        sync::Mutex,
    };

    use flate2::{read::MultiGzDecoder, write::DeflateEncoder, write::GzEncoder, Compression, Crc};
    use needletail::errors::ParseError;
    use rapidgzip_core::{DecodeError, Decoder};

    use crate::haec_io::HAECSeq;

    use super::{
        decode, encode, get_reads, rapidgzip_failure_message, rapidgzip_panic_message, GzipDecoder,
        HAECRecord, TrackedRead, RAPIDGZIP_THREADS,
    };

    // Deterministic FASTQ: ids with and without a description, mixed case, lengths 1..~3000.
    fn fastq(n: usize) -> Vec<u8> {
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut out = Vec::new();
        for i in 0..n {
            let len = 1 + (next() % 3000) as usize;
            if i % 3 == 0 {
                write!(out, "@read{i}\n").unwrap();
            } else {
                write!(out, "@read{i} desc\tx={i}\n").unwrap();
            }
            for _ in 0..len {
                out.push(b"ACGTacgt"[(next() % 8) as usize]);
            }
            out.extend_from_slice(b"\n+\n");
            for _ in 0..len {
                out.push(b'!' + (next() % 40) as u8);
            }
            out.push(b'\n');
        }
        out
    }

    fn gz(data: &[u8]) -> Vec<u8> {
        let mut e = GzEncoder::new(Vec::new(), Compression::new(6));
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    // BGZF: gzip members of <= 60000 input bytes with the BC extra field, then the EOF block.
    fn bgzf(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for block in data.chunks(60000) {
            let mut e = DeflateEncoder::new(Vec::new(), Compression::new(6));
            e.write_all(block).unwrap();
            let cdata = e.finish().unwrap();
            let mut crc = Crc::new();
            crc.update(block);
            let bsize = (18 + cdata.len() + 8 - 1) as u16;
            out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
            out.extend_from_slice(&bsize.to_le_bytes());
            out.extend_from_slice(&cdata);
            out.extend_from_slice(&crc.sum().to_le_bytes());
            out.extend_from_slice(&(block.len() as u32).to_le_bytes());
        }
        out.extend_from_slice(&[
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
        ]);
        out
    }

    fn test_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("herro-test-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write(dir: &Path, name: &str, data: &[u8]) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, data).unwrap();
        path
    }

    fn same(a: &[HAECRecord], b: &[HAECRecord]) -> bool {
        a.len() == b.len()
            && a.iter().zip(b).all(|(a, b)| {
                a.id == b.id && a.description == b.description && a.seq == b.seq && a.qual == b.qual
            })
    }

    #[test]
    fn both_decoders_give_identical_records() {
        let dir = test_dir("records");
        let fq = fastq(3000);
        // Three members (one empty); the cut falls mid-record, the decoded stream is the same.
        let (first, second) = fq.split_at(fq.len() / 2);
        let multi = [gz(first), gz(b""), gz(second)].concat();

        let plain = write(&dir, "reads.fq", &fq);
        let expected = get_reads(&plain, 500, &None, &None, GzipDecoder::Legacy);
        assert!(expected.len() > 1000);
        for (name, data) in [
            ("reads.fq", fq.clone()),
            ("gzip6.fq.gz", gz(&fq)),
            ("multi.fq.gz", multi),
            ("bgzf.fq.gz", bgzf(&fq)),
        ] {
            let path = write(&dir, name, &data);
            for decoder in [GzipDecoder::Legacy, GzipDecoder::Rapidgzip] {
                let reads = get_reads(&path, 500, &None, &None, decoder);
                assert!(same(&reads, &expected), "{name} {decoder:?}");
            }
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn classification_of_valid_file_is_library_problem() {
        let dir = test_dir("classify");
        let path = write(&dir, "reads.fq.gz", &gz(&fastq(200)));
        let injected = DecodeError::ChecksumMismatch {
            member: 0,
            expected: 1,
            actual: 2,
        };
        let msg = rapidgzip_failure_message(
            &path,
            Some((io::Error::new(io::ErrorKind::InvalidData, injected), 123)),
            ParseError::from(io::Error::other("parse")),
        );
        assert!(msg.contains("while the standard decoder reads the file fine"), "{msg}");
        assert!(msg.contains("rapidgzip-core problem"), "{msg}");
        assert!(msg.contains("--gzip-decoder legacy"), "{msg}");

        // A needletail parse error with no decoder error (bad bytes before the CRC check).
        let msg = rapidgzip_failure_message(&path, None, ParseError::new_empty_file());
        assert!(msg.contains("while the standard decoder reads the file fine"), "{msg}");

        // A corrupted file is reported as such, with legacy's error and offset first.
        let mut bad = gz(&fastq(200));
        bad.truncate(bad.len() / 2);
        let bad = write(&dir, "bad.fq.gz", &bad);
        let msg = rapidgzip_failure_message(
            &bad,
            Some((io::Error::other("injected"), 0)),
            ParseError::new_empty_file(),
        );
        assert!(msg.starts_with(&format!("Reads file {} is corrupted", bad.display())), "{msg}");
        assert!(msg.contains("(at decompressed byte "), "{msg}");
        assert!(msg.ends_with("\nrapidgzip-core 0.3.1 error: injected"), "{msg}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn panic_message_only_for_rapidgzip_core() {
        let lib = "/home/u/.cargo/registry/src/index.crates.io-x/rapidgzip-core-0.3.1/src/backend.rs";
        assert_eq!(
            rapidgzip_panic_message(Some(lib)),
            Some("rapidgzip-core 0.3.1 crashed while decoding the reads file. This is a rapidgzip-core problem; rerun with --gzip-decoder legacy.")
        );
        assert_eq!(rapidgzip_panic_message(Some("src/haec_io.rs")), None);
        assert_eq!(rapidgzip_panic_message(None), None);
    }

    // Returns Interrupted `n` times, then one byte per read.
    struct Interrupting(usize);

    impl Read for Interrupting {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.0 > 0 {
                self.0 -= 1;
                return Err(io::ErrorKind::Interrupted.into());
            }
            buf[0] = b'@';
            Ok(1)
        }
    }

    #[test]
    fn interrupted_is_retried_then_final() {
        let error = Mutex::new(None);
        let mut r = TrackedRead {
            inner: Interrupting(3),
            delivered: 0,
            error: &error,
        };
        assert_eq!(r.read(&mut [0; 4]).unwrap(), 1);
        assert!(error.lock().unwrap().is_none());

        r.inner = Interrupting(usize::MAX);
        let e = r.read(&mut [0; 4]).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::Other);
        let (stored, offset) = error.lock().unwrap().take().unwrap();
        assert_eq!(stored.kind(), io::ErrorKind::Interrupted);
        assert_eq!(offset, 1);
    }

    // Known difference: a dynamic block declaring 31 distance codes (HDIST 30), which
    // zlib-rs rejects, is decoded by rapidgzip-core when it takes its parallel (MarkerWindow) path.
    // Test file: one such block (35 MB; 83,890,080 B decoded, CRC32 and ISIZE checked by the
    // library). Ignored: tests the library, depends on its timing-based path choice, needs the file
    // outside the crate.
    // Run: HERRO_HDIST31_GZ=<file> cargo test --lib rapidgzip_accepts_hdist31 -- --ignored
    #[test]
    #[ignore]
    fn rapidgzip_accepts_hdist31() {
        let path =
            std::env::var("HERRO_HDIST31_GZ").expect("set HERRO_HDIST31_GZ to the test file");
        let path = path.as_str();
        assert!(Path::new(path).is_file(), "missing test file {path}");

        let mut legacy = MultiGzDecoder::new(std::fs::File::open(path).unwrap());
        assert!(io::copy(&mut legacy, &mut io::sink()).is_err());

        let decoder = Decoder::builder()
            .decoder_threads(RAPIDGZIP_THREADS)
            .build()
            .unwrap();
        let mut reader = decoder.open(path).unwrap();
        let n = io::copy(&mut reader, &mut io::sink());
        assert_eq!(n.unwrap(), 83_890_080, "path {:?}", reader.handle().stats().path);
    }

    #[test]
    fn encode_sequence1() {
        let sequence = "ACGT";
        let result = encode(sequence.as_bytes());

        assert_eq!(result, Ok((vec![0b11100100], 4)));
    }

    #[test]
    fn encode_invalid_bases() {
        let mut seq = vec![b'A'; 40];
        for (b, pos) in [
            (b'N', 5),
            (b'R', 5),
            (0x80, 5),
            (0xFF, 5),
            (b'N', 0),
            (b'N', 39),
            (b'N', 31),
            (b'N', 32),
        ] {
            seq[pos] = b;
            assert_eq!(encode(&seq), Err(pos));
            seq[pos] = b'A';
        }

        // First invalid base is reported
        let seq = b"ACGTNACGRT";
        assert_eq!(encode(seq), Err(4));
    }

    #[test]
    fn encode_lowercase() {
        let upper: Vec<u8> = (0..100).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect();
        let lower = upper.to_ascii_lowercase();
        assert_eq!(encode(&lower), encode(&upper));
        assert!(encode(&upper).is_ok());
    }

    #[test]
    fn encode_sequence2() {
        let haec_seq = HAECSeq::from("ACGTACG".as_bytes());
        assert_eq!(haec_seq, HAECSeq::new(vec![0b10010011100100], 7));
    }

    #[test]
    fn decode_sequence1() {
        let encoded = vec![0b11100100];
        let length = 4;

        let mut buffer = vec![0; 4];
        decode(&encoded, length, .., false, &mut buffer);
        assert_eq!(&buffer[..4], "ACGT".as_bytes());
    }

    #[test]
    fn decode_sequence2() {
        let haec_seq = HAECSeq::new(vec![0b10010011100100], 7);
        assert_eq!(Vec::from(&haec_seq), "ACGTACG".as_bytes());
    }

    #[test]
    fn test_range1() {
        let haec_seq = HAECSeq::from("ACGTACGTACGT".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_subseq(3..10, &mut buffer);

        assert_eq!(&buffer[..7], "TACGTAC".as_bytes())
    }

    #[test]
    fn test_range2() {
        let haec_seq = HAECSeq::from("ACGTACGTACGT".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_subseq(3.., &mut buffer);

        assert_eq!(&buffer[..9], "TACGTACGT".as_bytes())
    }

    #[test]
    fn test_range3() {
        let haec_seq = HAECSeq::from("ACGTACGTACGT".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_subseq(3..haec_seq.len(), &mut buffer);

        assert_eq!(&buffer[..9], "TACGTACGT".as_bytes())
    }

    #[test]
    fn test_range4() {
        let haec_seq = HAECSeq::from("ACGTACGTACGT".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_subseq(.., &mut buffer);

        assert_eq!(&buffer[..12], "ACGTACGTACGT".as_bytes())
    }

    #[test]
    fn test_range5() {
        let haec_seq = HAECSeq::from("ACGTACGTACGT".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_subseq(8..11, &mut buffer);

        assert_eq!(&buffer[..3], "ACG".as_bytes())
    }

    #[test]
    fn test_rc1() {
        let haec_seq = HAECSeq::from("ATCGATCGATCG".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_rc_subseq(.., &mut buffer);

        assert_eq!(&buffer[..12], "CGATCGATCGAT".as_bytes())
    }

    #[test]
    fn test_rc2() {
        let haec_seq = HAECSeq::from("ATCGATCGATCG".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_rc_subseq(3.., &mut buffer);

        assert_eq!(&buffer[..9], "CGATCGATC".as_bytes())
    }

    #[test]
    fn test_rc3() {
        let haec_seq = HAECSeq::from("ATCGATCGATCG".as_bytes());

        let mut buffer = vec![0; 100];
        haec_seq.get_rc_subseq(..9, &mut buffer);

        assert_eq!(&buffer[..9], "TCGATCGAT".as_bytes())
    }
}
