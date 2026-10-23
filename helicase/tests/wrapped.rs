// Wrapped sequences go through the per-line stride path, unwrapped ones never do:
// both must give the same records, and split chunks must match the non-ACTG runs.
use fastx_gen::{Alphabet, FastxGen};
use helicase::config::advanced::*;
use helicase::config::*;
use helicase::input::*;
use helicase::parser::{FastaParser, FastqParser};
use helicase::*;
use std::io::Cursor;

mod common;
use common::{assert_same, bits};

const STRING: Config = ParserOptions::default().compute_dna_len().config();
const KEEP: Config = ParserOptions::default()
    .and_dna_packed()
    .and_dna_columnar()
    .keep_non_actg()
    .compute_mask_non_actg()
    .compute_dna_len()
    .config()
    | COMPUTE_MASK_N;
const KEEP_PACKED_ONLY: Config = ParserOptions::default()
    .dna_packed()
    .keep_non_actg()
    .config();
const SKIP: Config = ParserOptions::default()
    .and_dna_packed()
    .and_dna_columnar()
    .skip_non_actg()
    .compute_dna_len()
    .config();
const SPLIT: Config = ParserOptions::default()
    .and_dna_packed()
    .and_dna_columnar()
    .config();
const SPLIT_PACKED_ONLY: Config = ParserOptions::default().dna_packed().config();
const STRING_SPLIT: Config = ParserOptions::default().split_non_actg().config();

#[derive(Debug, PartialEq, Eq)]
struct Output {
    header: Option<Vec<u8>>,
    dna: Option<Vec<u8>>,
    packed: Option<String>,
    columnar: Option<String>,
    masks: Option<(Vec<bool>, Vec<bool>)>,
    dna_len: Option<usize>,
}

fn parse_all<const CONFIG: Config, I: InputData<'static>>(
    mut parser: FastaParser<'static, CONFIG, I>,
) -> Vec<Output> {
    let set = |flag| CONFIG & flag != 0;
    let mut out = Vec::new();
    while let Some(event) = parser.next() {
        event.unwrap();
        out.push(Output {
            header: set(COMPUTE_HEADER).then(|| parser.get_header().to_vec()),
            dna: set(COMPUTE_DNA_STRING).then(|| parser.get_dna_string().to_vec()),
            packed: set(COMPUTE_DNA_PACKED).then(|| parser.get_dna_packed().to_string()),
            columnar: set(COMPUTE_DNA_COLUMNAR).then(|| parser.get_dna_columnar().to_string()),
            masks: set(COMPUTE_MASK_NON_ACTG)
                .then(|| (bits(parser.get_mask_non_actg()), bits(parser.get_mask_n()))),
            dna_len: set(COMPUTE_DNA_LEN).then(|| parser.get_dna_len()),
        });
    }
    out
}

fn check_config<const CONFIG: Config>(unwrapped: &'static [u8], wrapped: &'static [u8], id: &str) {
    let expected = parse_all(FastaParser::<CONFIG, _>::from_slice(unwrapped).unwrap());
    let from_slice = parse_all(FastaParser::<CONFIG, _>::from_slice(wrapped).unwrap());
    assert_same(&expected, &from_slice, &format!("slice {id}"));
    let from_reader =
        parse_all(FastaParser::<CONFIG, _>::from_reader(Cursor::new(wrapped)).unwrap());
    assert_same(&expected, &from_reader, &format!("reader {id}"));
}

#[test]
fn test_wrapped_matches_unwrapped() {
    for width in [20, 59, 60, 61, 62, 63, 64, 65, 80, 126, 127, 128, 150] {
        for crlf in [false, true] {
            let generator = FastxGen::fasta(60, 3000)
                .with_std_dev(0.5)
                .with_alphabet(Alphabet::DnaIupac)
                .with_non_actg_prob(0.002)
                .with_crlf(crlf)
                .with_seed(width as u64);
            let unwrapped = generator.generate_bytes().leak();
            let wrapped = generator.with_wrap(width).generate_bytes().leak();
            let id = format!("width {width} crlf {crlf}");
            check_config::<STRING>(unwrapped, wrapped, &format!("string {id}"));
            check_config::<KEEP>(unwrapped, wrapped, &format!("keep {id}"));
            check_config::<KEEP_PACKED_ONLY>(unwrapped, wrapped, &format!("keep packed {id}"));
            check_config::<SKIP>(unwrapped, wrapped, &format!("skip {id}"));
            check_config::<SPLIT>(unwrapped, wrapped, &format!("split {id}"));
            check_config::<SPLIT_PACKED_ONLY>(unwrapped, wrapped, &format!("split packed {id}"));
            check_config::<STRING_SPLIT>(unwrapped, wrapped, &format!("string split {id}"));
        }
    }
}

fn chunks<const CONFIG: Config, P>(mut parser: P) -> Vec<String>
where
    P: HelicaseParser + Iterator<Item = parser::Result<parser::Event>>,
{
    let mut out = Vec::new();
    while let Some(event) = parser.next() {
        event.unwrap();
        out.push(if CONFIG & COMPUTE_DNA_STRING != 0 {
            String::from_utf8_lossy(parser.get_dna_string()).to_string()
        } else {
            parser.get_dna_packed().to_string()
        });
    }
    out
}

fn check_split<const CONFIG: Config>(data: &'static [u8], expected: &[String], id: &str) {
    let got = if data[0] == b'>' {
        chunks::<CONFIG, _>(FastaParser::<CONFIG, _>::from_slice(data).unwrap())
    } else {
        chunks::<CONFIG, _>(FastqParser::<CONFIG, _>::from_slice(data).unwrap())
    };
    assert_eq!(got, expected, "{id}");
}

#[test]
fn test_split_chunks() {
    for (width, crlf) in [(4, false), (60, false), (4, true), (60, true)] {
        let line_end: &[u8] = if crlf { b"\r\n" } else { b"\n" };
        for len in [2 * width, 3 * width + 1, 200, 3000] {
            for first in (0..len).step_by(if len > 500 { 97 } else { 1 }) {
                for second in [first, first + 1, first + width, first + width + 1] {
                    let seq: Vec<u8> = (0..len)
                        .map(|i| match (i == first, i == second) {
                            (true, _) => b'N',
                            // A lone `\r` still splits the sequence in CRLF files.
                            (_, true) if crlf => b'\r',
                            (_, true) => b'N',
                            _ => b"ACGT"[(i * 7 + i / 3) % 4],
                        })
                        .collect();
                    let chunks: Vec<String> = seq
                        .split(|&c| c == b'N' || c == b'\r')
                        .filter(|chunk| !chunk.is_empty())
                        .map(|chunk| String::from_utf8_lossy(chunk).to_string())
                        .collect();
                    let with_next = [chunks.clone(), vec!["ACGT".to_string()]].concat();
                    let lines = |lines: &[&[u8]]| -> Vec<u8> {
                        lines.iter().flat_map(|l| [*l, line_end].concat()).collect()
                    };
                    let unwrapped = lines(&[b">a", &seq]);
                    let mut wrapped = lines(&[b">a"]);
                    for line in seq.chunks(width) {
                        wrapped.extend(lines(&[line]));
                    }
                    let qual = vec![b'I'; len];
                    let fastq = lines(&[b"@a", &seq, b"+", &qual]);
                    let next_fasta = lines(&[b">b", b"ACGT"]);
                    let next_fastq = lines(&[b"@b", b"ACGT", b"+", b"IIII"]);
                    for (name, record, next) in [
                        ("unwrapped", unwrapped, &next_fasta),
                        ("wrapped", wrapped, &next_fasta),
                        ("fastq", fastq, &next_fastq),
                    ] {
                        let no_final_lf = record[..record.len() - line_end.len()].to_vec();
                        for (end, data, expected) in [
                            (
                                "next record",
                                [record.clone(), next.clone()].concat(),
                                &with_next,
                            ),
                            ("last record", record, &chunks),
                            ("no final lf", no_final_lf, &chunks),
                        ] {
                            let data = data.leak();
                            let id = format!(
                                "{name} {end} crlf {crlf} width {width} len {len} N at {first} and {second}"
                            );
                            check_split::<STRING_SPLIT>(data, expected, &id);
                            check_split::<SPLIT_PACKED_ONLY>(data, expected, &id);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn test_wrapped_buffers_fit_records() {
    const LEN: usize = 50_000;
    let data = FastxGen::wrapped_fasta(40, LEN, 60)
        .with_std_dev(0.0)
        .generate_bytes();
    let mut parser = FastaParser::<STRING, _>::from_slice(&data).unwrap();
    parser.next().unwrap().unwrap();
    let dna = parser.get_dna_string_owned();
    assert_eq!(dna.len(), LEN);
    assert!(dna.capacity() < 4 * LEN, "capacity {}", dna.capacity());
    let mut parser = FastaParser::<KEEP_PACKED_ONLY, _>::from_slice(&data).unwrap();
    parser.next().unwrap().unwrap();
    let packed = parser.get_dna_packed_owned();
    assert_eq!(packed.len(), LEN);
    // In 64-bit words of 32 bases.
    assert!(
        32 * packed.capacity() < 4 * LEN,
        "capacity {}",
        packed.capacity()
    );
}

#[test]
fn test_wrapped_keeps_lone_carriage_return() {
    // After a stride over full lines, a `\r` not followed by `\n` is a base, like in the block parser.
    let line = "ACGT".repeat(15);
    let last = format!("\r{}", &line[1..]);
    let wrapped = format!(">r\n{}{last}\n", format!("{line}\n").repeat(8));
    let mut parser = FastxParser::<STRING>::from_slice(wrapped.as_bytes()).unwrap();
    parser.next().unwrap().unwrap();
    assert_eq!(parser.get_dna_string(), (line.repeat(8) + &last).as_bytes());
}
