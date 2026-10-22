// Slice inputs parse short records from a line feed index, reader inputs never do:
// both must give the same output, including on malformed inputs.
use fastx_gen::{Alphabet, FastxGen};
use helicase::config::*;
use helicase::input::*;
use helicase::parser::{FastaParser, FastqParser};
use helicase::*;
use std::io::Cursor;

mod common;
use common::{assert_same, bits};

const STRING: Config = ParserOptions::default()
    .compute_quality()
    .compute_dna_len()
    .config();
const PACKED: Config = ParserOptions::default()
    .compute_quality()
    .compute_dna_len()
    .and_dna_packed()
    .and_dna_columnar()
    .keep_non_actg()
    .compute_mask_non_actg()
    .config()
    | advanced::COMPUTE_MASK_N;

#[derive(Debug, PartialEq, Eq)]
enum Parsed {
    Record {
        event: parser::Event,
        header: Vec<u8>,
        dna: Vec<u8>,
        qual: Option<Vec<u8>>,
        dna_len: usize,
        packed: Option<(String, String, Vec<bool>, Vec<bool>)>,
    },
    Error(String),
}

fn parse_all<const CONFIG: Config, P>(mut parser: P) -> Vec<Parsed>
where
    P: HelicaseParser + Iterator<Item = parser::Result<parser::Event>>,
{
    let mut out = Vec::new();
    while let Some(event) = parser.next() {
        match event {
            Ok(event) => out.push(Parsed::Record {
                event,
                header: parser.get_header_owned(),
                dna: parser.get_dna_string_owned(),
                qual: parser.get_quality_owned(),
                dna_len: parser.get_dna_len(),
                packed: (CONFIG == PACKED).then(|| {
                    (
                        parser.get_dna_packed().to_string(),
                        parser.get_dna_columnar().to_string(),
                        bits(parser.get_mask_non_actg()),
                        bits(parser.get_mask_n()),
                    )
                }),
            }),
            Err(e) => {
                out.push(Parsed::Error(e.to_string()));
                break;
            }
        }
    }
    out
}

fn check_config<const CONFIG: Config>(data: &[u8], test_id: &str) {
    let to_string = |e: parser::Error| e.to_string();
    let (from_slice, from_reader) = if data.first() == Some(&b'>') {
        (
            FastaParser::<CONFIG, _>::from_slice(data).map(parse_all::<CONFIG, _>),
            FastaParser::<CONFIG, _>::from_reader(Cursor::new(data)).map(parse_all::<CONFIG, _>),
        )
    } else {
        (
            FastqParser::<CONFIG, _>::from_slice(data).map(parse_all::<CONFIG, _>),
            FastqParser::<CONFIG, _>::from_reader(Cursor::new(data)).map(parse_all::<CONFIG, _>),
        )
    };
    match (from_slice, from_reader) {
        (Ok(slice), Ok(reader)) => assert_same(&slice, &reader, test_id),
        (slice, reader) => assert_eq!(slice.map_err(to_string), reader.map_err(to_string)),
    }
}

fn check(data: &[u8], test_id: &str) {
    check_config::<STRING>(data, test_id);
    check_config::<PACKED>(data, test_id);
}

// Runs of short records separated by long or wrapped ones, to switch in and out of the index.
fn generate(fasta: bool, seed: u64, crlf: bool) -> Vec<u8> {
    let mut data = Vec::new();
    for (i, (count, len)) in [(40, 150), (2, 2000), (30, 300), (1, 700), (40, 100)]
        .into_iter()
        .enumerate()
    {
        let generator = if fasta {
            FastxGen::fasta(count, len).with_wrap(if i == 3 { 60 } else { 0 })
        } else {
            FastxGen::fastq(count, len)
        };
        let records = generator
            .with_std_dev(0.3)
            .with_alphabet(Alphabet::DnaIupac)
            .with_non_actg_prob(0.01)
            .with_crlf(crlf)
            .with_seed(seed * 8 + i as u64)
            .generate_bytes();
        data.extend_from_slice(&records);
    }
    data
}

fn line_starts(data: &[u8]) -> Vec<usize> {
    let mut starts = vec![0];
    starts.extend(memchr::memchr_iter(b'\n', data).map(|i| i + 1));
    starts.pop();
    starts
}

#[test]
fn test_indexed_mixed_lengths() {
    for fasta in [false, true] {
        for seed in 0..8 {
            for crlf in [false, true] {
                let data = generate(fasta, seed, crlf);
                let id = format!("fasta {fasta} seed {seed} crlf {crlf}");
                check(&data, &id);
                check(&data[..data.len() - 1], &format!("{id} no final lf"));
            }
        }
    }
}

#[test]
fn test_indexed_truncated() {
    for fasta in [false, true] {
        let data = generate(fasta, 1, false);
        for end in data.len() - 600..data.len() {
            check(&data[..end], &format!("fasta {fasta} truncated at {end}"));
        }
    }
}

#[test]
fn test_indexed_malformed_lines() {
    for fasta in [false, true] {
        for crlf in [false, true] {
            let line_end: &[u8] = if crlf { b"\r\n" } else { b"\n" };
            let data = generate(fasta, 2, crlf);
            let starts = line_starts(&data);
            for (k, &start) in starts.iter().enumerate().step_by(5) {
                let id = format!("fasta {fasta} crlf {crlf} line {k}");
                let mut extra = data.clone();
                extra.insert(start + 1, b'I');
                check(&extra, &format!("{id} extra byte"));

                let mut missing = data.clone();
                missing.remove(start + 1);
                check(&missing, &format!("{id} missing byte"));

                let mut carriage_return = data.clone();
                carriage_return.insert(start, b'\r');
                check(&carriage_return, &format!("{id} carriage return"));

                if let Some(&end) = starts.get(k + 3) {
                    let mut no_lines = data.clone();
                    no_lines.drain(starts[k + 2]..end);
                    no_lines.drain(start..starts[k + 1]);
                    check(&no_lines, &format!("{id} lines removed"));
                }

                // A blank line before a FASTQ header panics in the block parser.
                if fasta || k % 4 != 0 {
                    let mut blank = data.clone();
                    blank.splice(start..start, line_end.iter().copied());
                    check(&blank, &format!("{id} blank line before"));
                }
            }
        }
    }
}

#[test]
fn test_indexed_quality_jump() {
    let data = FastxGen::fastq(60, 300).with_seed(3).generate_bytes();
    let starts = line_starts(&data);
    for record in (8..56).step_by(5) {
        // Shorten a quality line so that its expected end falls on the next header's line feed,
        // in a later block than the actual end.
        let next_header = starts[4 * record + 4];
        let mut short = data.clone();
        short.splice(next_header + 1..next_header + 1, [b'x'; 100]);
        let header_len = starts[4 * record + 5] + 100 - next_header - 1;
        let qual_end = next_header - 1;
        short.drain(qual_end - header_len - 1..qual_end);
        check(&short, &format!("short quality line in record {record}"));
    }
}
