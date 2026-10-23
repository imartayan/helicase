#![cfg(feature = "std")]

use helicase::config::*;
use helicase::dna_format::{ColumnarDNA, PackedDNA};
use helicase::input::*;
use helicase::parser::{Error, ErrorKind, Event};
use helicase::*;
use std::io::Cursor;

const STRING: Config = ParserOptions::default().reject_non_actg().config();
const PACKED: Config = ParserOptions::default()
    .dna_packed()
    .reject_non_actg()
    .config();
const COLUMNAR: Config = ParserOptions::default()
    .dna_columnar()
    .reject_non_actg()
    .config();
const STRING_PACKED: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .reject_non_actg()
    .config();

const SENTINEL: u8 = 0x01;

#[derive(Clone, Copy, Debug)]
enum Format {
    Fasta,
    Wrapped,
    Fastq,
}

// A record parsed as (header, sequence), or the error as (byte, offset).
type Parsed = Result<(Vec<u8>, Vec<u8>), (u8, usize)>;

fn random_seq(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGTacgt"[(state >> 61) as usize]
        })
        .collect()
}

// Headers and quality lines contain non-ACTG bytes that must be accepted.
fn build(seqs: &[Vec<u8>], format: Format, crlf: bool, trailing_newline: bool) -> Vec<u8> {
    let eol: &[u8] = if crlf { b"\r\n" } else { b"\n" };
    let mut out = Vec::new();
    for (i, seq) in seqs.iter().enumerate() {
        let marker = if matches!(format, Format::Fastq) {
            b'@'
        } else {
            b'>'
        };
        out.push(marker);
        out.extend_from_slice(format!("r{i} N-R.*").as_bytes());
        out.extend_from_slice(eol);
        match format {
            Format::Fasta | Format::Fastq => {
                out.extend_from_slice(seq);
                out.extend_from_slice(eol);
            }
            Format::Wrapped => {
                for line in seq.chunks(60) {
                    out.extend_from_slice(line);
                    out.extend_from_slice(eol);
                }
            }
        }
        if matches!(format, Format::Fastq) {
            out.push(b'+');
            out.extend_from_slice(eol);
            out.extend(b"N#I".iter().cycle().take(seq.len()));
            out.extend_from_slice(eol);
        }
    }
    if !trailing_newline {
        out.truncate(out.len() - eol.len());
    }
    out
}

fn collect<P>(mut parser: P, dna: impl Fn(&P) -> Vec<u8>) -> Vec<Parsed>
where
    P: HelicaseParser + Iterator<Item = Result<Event, Error>>,
{
    let mut out = Vec::new();
    while let Some(event) = parser.next() {
        match event {
            Ok(Event::Record(_)) => out.push(Ok((parser.get_header().to_vec(), dna(&parser)))),
            Ok(Event::DnaChunk(_)) => panic!("unexpected chunk"),
            Err(e) => match e.kind() {
                ErrorKind::InvalidBase { byte, pos } => {
                    out.push(Err((*byte, *pos)));
                    assert!(parser.next().is_none(), "parser must stop after the error");
                    break;
                }
                kind => panic!("unexpected error: {kind}"),
            },
        }
    }
    out
}

fn parse_all(input: &[u8]) -> Vec<(&'static str, Vec<Parsed>)> {
    macro_rules! run {
        ($name:literal, $cfg:ident, $dna:expr) => {
            [
                (
                    concat!($name, " slice"),
                    collect(FastxParser::<$cfg>::from_slice(input).unwrap(), $dna),
                ),
                (
                    concat!($name, " reader"),
                    collect(
                        FastxParser::<$cfg>::from_reader(Cursor::new(input)).unwrap(),
                        $dna,
                    ),
                ),
            ]
        };
    }
    let packed = |dna: &PackedDNA| (0..dna.len()).map(|i| dna.get_char(i) as u8).collect();
    let columnar = |dna: &ColumnarDNA| (0..dna.len()).map(|i| dna.get_char(i) as u8).collect();
    [
        run!("string", STRING, |p| p.get_dna_string().to_vec()),
        run!("packed", PACKED, |p| packed(p.get_dna_packed())),
        run!("columnar", COLUMNAR, |p| columnar(p.get_dna_columnar())),
        run!("string+packed", STRING_PACKED, |p| {
            let s = p.get_dna_string().to_vec();
            assert_eq!(s.to_ascii_uppercase(), packed(p.get_dna_packed()));
            s
        }),
    ]
    .concat()
}

fn expected(seqs: &[Vec<u8>], uppercase: bool) -> Vec<Parsed> {
    seqs.iter()
        .enumerate()
        .map(|(i, seq)| {
            let header = format!("r{i} N-R.*").into_bytes();
            let seq = if uppercase {
                seq.to_ascii_uppercase()
            } else {
                seq.clone()
            };
            Ok((header, seq))
        })
        .collect()
}

fn check(seqs: &[Vec<u8>], input: &[u8], error: Option<(u8, usize)>, valid_records: usize) {
    for (name, got) in parse_all(input) {
        let uppercase = name.starts_with("packed") || name.starts_with("columnar");
        let mut want = expected(&seqs[..valid_records], uppercase);
        if let Some(err) = error {
            want.push(Err(err));
        }
        assert_eq!(got.len(), want.len(), "{name}: record count");
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            match (g, w) {
                (Ok(g), Ok(w)) => assert!(g == w, "{name}: record {i} differs"),
                _ => assert_eq!(g, w, "{name}: record {i}"),
            }
        }
    }
}

#[test]
fn test_reject_non_actg_accepts_valid_input() {
    for format in [Format::Fasta, Format::Wrapped, Format::Fastq] {
        for crlf in [false, true] {
            for trailing_newline in [true, false] {
                let seqs: Vec<Vec<u8>> = [150, 1, 64, 15_000, 300_000, 150]
                    .iter()
                    .enumerate()
                    .map(|(i, &len)| random_seq(len, i as u64))
                    .collect();
                let input = build(&seqs, format, crlf, trailing_newline);
                check(&seqs, &input, None, seqs.len());
            }
        }
    }
}

#[test]
fn test_reject_non_actg_reports_first_invalid_base() {
    for format in [Format::Fasta, Format::Wrapped, Format::Fastq] {
        for crlf in [false, true] {
            for len in [150, 15_000, 300_000] {
                let bytes: &[u8] = if len < 100_000 { b"NnR-.*\r" } else { b"N\r" };
                for pos in [0, 59, 60, 63, 64, len / 2, len - 1] {
                    for &byte in bytes {
                        let mut seqs = vec![random_seq(150, 1), random_seq(len, 2)];
                        seqs.push(random_seq(150, 3));
                        seqs[1][pos] = SENTINEL;
                        let mut input = build(&seqs, format, crlf, true);
                        let offset = input.iter().position(|&b| b == SENTINEL).unwrap();
                        input[offset] = byte;
                        let ctx = format!("{format:?} crlf={crlf} len={len} pos={pos} byte={byte}");
                        let result = std::panic::catch_unwind(|| {
                            check(&seqs, &input, Some((byte, offset)), 1)
                        });
                        assert!(result.is_ok(), "{ctx}");
                    }
                }
            }
        }
    }
}

#[test]
fn test_reject_non_actg_at_end_of_input() {
    for format in [Format::Fasta, Format::Wrapped, Format::Fastq] {
        for crlf in [false, true] {
            for trailing_newline in [true, false] {
                let mut seqs = vec![random_seq(150, 1), random_seq(200, 2)];
                seqs[1][199] = SENTINEL;
                let mut input = build(&seqs, format, crlf, trailing_newline);
                let offset = input.iter().position(|&b| b == SENTINEL).unwrap();
                input[offset] = b'N';
                check(&seqs, &input, Some((b'N', offset)), 1);
            }
        }
    }
}
