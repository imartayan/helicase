// Random inputs, compared against Needletail when valid, and checked for panics and endless
// iteration when corrupted. `HELICASE_RANDOM_CASES` sets the number of cases (default 1000),
// `HELICASE_RANDOM_SEED` offsets them for long runs, see `NOTES/RANDOM_TESTING.md`.
#![cfg(feature = "std")]

use helicase::config::*;
use helicase::input::*;
use helicase::parser::{Error, ErrorKind, Event};
use helicase::*;
use std::io::Cursor;
use std::panic::{AssertUnwindSafe, catch_unwind};

const KEEP: Config = ParserOptions::default()
    .compute_quality()
    .and_dna_packed()
    .keep_non_actg()
    .config();
const SPLIT: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .config();
const SKIP: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .skip_non_actg()
    .config();
const REJECT: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .reject_non_actg()
    .config();

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn pick(&mut self, bytes: &[u8]) -> u8 {
        bytes[self.below(bytes.len())]
    }
}

type Record = (Vec<u8>, Vec<u8>, Option<Vec<u8>>);

fn cases() -> core::ops::Range<u64> {
    let var = |name, default| std::env::var(name).map_or(default, |n: String| n.parse().unwrap());
    let seed = var("HELICASE_RANDOM_SEED", 0);
    seed..seed + var("HELICASE_RANDOM_CASES", 1000)
}

fn is_actg(b: u8) -> bool {
    matches!(b & 0xDF, b'A' | b'C' | b'G' | b'T')
}

fn random_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
    let noise = [0, 0, 1, 50][rng.below(4)];
    (0..len)
        .map(|_| {
            if rng.below(1000) < noise {
                rng.pick(b"NnRy-.*")
            } else {
                rng.pick(b"ACGTACGTACGTacgt")
            }
        })
        .collect()
}

fn push_line(out: &mut Vec<u8>, line: &[u8], eol: &[u8]) {
    out.extend_from_slice(line);
    out.extend_from_slice(eol);
}

// One line ending per file, no lone `\r`, no blank line between FASTQ records.
fn random_input(rng: &mut Rng) -> Vec<u8> {
    let fastq = rng.below(2) == 0;
    let eol: &[u8] = if rng.below(4) == 0 { b"\r\n" } else { b"\n" };
    let big = rng.below(100) == 0;
    let mut out = Vec::new();
    let records = 1 + rng.below(8);
    for i in 0..records {
        let len = match rng.below(8) {
            // Crosses the 256 KiB reader buffers.
            _ if big && i == 0 => 200_000 + rng.below(200_000),
            0 => rng.below(5000),
            _ => rng.below(300),
        };
        // Needletail rejects an empty last record.
        let len = if i + 1 == records { len.max(1) } else { len };
        let seq = random_seq(rng, len);
        let marker = if fastq { "@" } else { ">" };
        push_line(&mut out, format!("{marker}r{i} x").as_bytes(), eol);
        if fastq {
            push_line(&mut out, &seq, eol);
            push_line(
                &mut out,
                if rng.below(4) == 0 { b"+r x" } else { b"+" },
                eol,
            );
            let qual: Vec<u8> = (0..len).map(|_| b'!' + rng.below(42) as u8).collect();
            push_line(&mut out, &qual, eol);
            continue;
        }
        let width = match rng.below(4) {
            0 => len.max(1),
            _ => 1 + rng.below(150),
        };
        let irregular = rng.below(8) == 0;
        let mut rest = &seq[..];
        while !rest.is_empty() {
            let n = if irregular {
                1 + rng.below(2 * width)
            } else {
                width
            };
            let (line, tail) = rest.split_at(n.min(rest.len()));
            push_line(&mut out, line, eol);
            rest = tail;
            if rng.below(50) == 0 {
                out.extend_from_slice(eol);
            }
        }
    }
    if rng.below(4) == 0 {
        out.truncate(out.len() - eol.len());
    }
    out
}

fn needletail_records(input: &[u8]) -> Vec<Record> {
    let mut reader = needletail::parse_fastx_reader(input).unwrap();
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        records.push((
            record.id().to_vec(),
            record.seq().to_vec(),
            record.qual().map(<[u8]>::to_vec),
        ));
    }
    records
}

fn packed_string<P: HelicaseParser>(parser: &P) -> Vec<u8> {
    let dna = parser.get_dna_packed();
    (0..dna.len()).map(|i| dna.get_char(i) as u8).collect()
}

fn actg_runs(seq: &[u8]) -> impl Iterator<Item = &[u8]> {
    seq.split(|&b| !is_actg(b)).filter(|run| !run.is_empty())
}

fn check_keep<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(
    mut parser: P,
    expected: &[Record],
) {
    let mut n = 0;
    while let Some(event) = parser.next() {
        event.unwrap();
        let (header, seq, qual) = &expected[n];
        assert_eq!(parser.get_header(), header, "record {n} header");
        assert_eq!(parser.get_dna_string(), seq, "record {n} seq");
        assert_eq!(parser.get_quality(), qual.as_deref(), "record {n} quality");
        assert_eq!(
            parser.get_dna_packed().len(),
            seq.len(),
            "record {n} packed"
        );
        n += 1;
    }
    assert_eq!(n, expected.len(), "record count");
}

fn check_split<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(
    mut parser: P,
    expected: &[Record],
) {
    let mut runs = expected.iter().flat_map(|(_, seq, _)| actg_runs(seq));
    while let Some(event) = parser.next() {
        event.unwrap();
        let chunk = parser.get_dna_string();
        if chunk.is_empty() {
            continue;
        }
        let run = runs.next().expect("more chunks than ACTG runs");
        assert_eq!(chunk, run, "chunk");
        assert_eq!(
            packed_string(&parser),
            run.to_ascii_uppercase(),
            "packed chunk"
        );
    }
    assert!(runs.next().is_none(), "fewer chunks than ACTG runs");
}

fn check_skip<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(
    mut parser: P,
    expected: &[Record],
) {
    let mut n = 0;
    while let Some(event) = parser.next() {
        event.unwrap();
        let merged: Vec<u8> = actg_runs(&expected[n].1).flatten().copied().collect();
        assert_eq!(parser.get_dna_string(), merged, "record {n} merged");
        assert_eq!(
            packed_string(&parser),
            merged.to_ascii_uppercase(),
            "record {n} packed"
        );
        n += 1;
    }
    assert_eq!(n, expected.len(), "record count");
}

fn check_reject<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(
    mut parser: P,
    expected: &[Record],
    input: &[u8],
) {
    for (n, (header, seq, _)) in expected.iter().enumerate() {
        let event = parser.next().expect("missing record");
        if let Some(&byte) = seq.iter().find(|&&b| !is_actg(b)) {
            let err = event.expect_err("non-ACTG base accepted");
            let ErrorKind::InvalidBase { byte: b, pos } = *err.kind() else {
                panic!("unexpected error: {err}");
            };
            assert_eq!((b, input[pos]), (byte, byte), "record {n} invalid base");
            assert!(parser.next().is_none(), "parser must stop after the error");
            return;
        }
        event.unwrap();
        assert_eq!(parser.get_header(), header, "record {n} header");
        assert_eq!(parser.get_dna_string(), seq, "record {n} seq");
    }
    assert!(parser.next().is_none(), "extra record");
}

macro_rules! both_inputs {
    ($cfg:ident, $input:expr, |$parser:ident| $body:expr) => {{
        let $parser = FastxParser::<$cfg>::from_slice($input).unwrap();
        $body;
        let $parser = FastxParser::<$cfg>::from_reader(Cursor::new($input)).unwrap();
        $body;
    }};
}

fn run_case(case: u64, input: &[u8], check: impl FnOnce()) {
    if catch_unwind(AssertUnwindSafe(check)).is_err() {
        let path = std::env::temp_dir().join(format!("helicase_random_{case}.fx"));
        std::fs::write(&path, input).unwrap();
        panic!("case {case} failed, input written to {}", path.display());
    }
}

#[test]
fn test_random_inputs_match_needletail() {
    for case in cases() {
        let mut rng = Rng(0x9E3779B97F4A7C15 ^ case);
        let input = random_input(&mut rng);
        run_case(case, &input, || {
            let expected = needletail_records(&input);
            both_inputs!(KEEP, &input, |p| check_keep(p, &expected));
            both_inputs!(SPLIT, &input, |p| check_split(p, &expected));
            both_inputs!(SKIP, &input, |p| check_skip(p, &expected));
            both_inputs!(REJECT, &input, |p| check_reject(p, &expected, &input));
        });
    }
}

// Calls every getter after each event, and bounds the number of events by the input length.
fn drain<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(mut parser: P, len: usize) {
    for _ in 0..len + 2 {
        match parser.next() {
            None => return,
            Some(Err(_)) => {}
            Some(Ok(Event::Record(_))) => {
                parser.get_header();
                parser.get_dna_string();
                parser.get_dna_packed();
            }
            Some(Ok(Event::DnaChunk(_))) => {
                parser.get_dna_string();
                parser.get_dna_packed();
            }
        }
    }
    panic!("more events than input bytes");
}

#[test]
fn test_random_corrupted_inputs_do_not_panic() {
    for case in cases() {
        let mut rng = Rng(0xD1B54A32D192ED03 ^ case);
        let mut input = random_input(&mut rng);
        for _ in 0..1 + rng.below(8) {
            let pos = rng.below(input.len() + 1);
            let byte = rng.pick(b">@+\n\n\r ACGTN");
            match rng.below(3) {
                0 if pos < input.len() => input[pos] = byte,
                1 if pos < input.len() => drop(input.remove(pos)),
                _ => input.insert(pos, byte),
            }
        }
        if !matches!(input.first(), Some(b'>' | b'@')) {
            continue;
        }
        run_case(case, &input, || {
            let len = input.len();
            both_inputs!(KEEP, &input, |p| drain(p, len));
            both_inputs!(SPLIT, &input, |p| drain(p, len));
            both_inputs!(SKIP, &input, |p| drain(p, len));
            both_inputs!(REJECT, &input, |p| drain(p, len));
        });
    }
}
