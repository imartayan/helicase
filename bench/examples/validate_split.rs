use helicase::input::*;
use helicase::*;
use needletail::parse_fastx_file;
use regex::bytes::{Regex, RegexBuilder};

use std::sync::LazyLock;

const CONFIG: Config = ParserOptions::default()
    .compute_quality()
    .split_non_actg()
    .config();

static MATCH_N: LazyLock<Regex> = LazyLock::new(|| {
    RegexBuilder::new(r"[N]+")
        .case_insensitive(true)
        .unicode(false)
        .build()
        .unwrap()
});

fn check_mismatch(left: &[u8], right: &[u8]) -> Option<usize> {
    if left == right {
        return None;
    }
    let len = left.len().min(right.len());
    (0..len).find(|&i| left[i] != right[i]).or(Some(len))
}

fn get_scope(slice: &[u8], pos: usize) -> &str {
    let start = pos.saturating_sub(10);
    let stop = (pos + 5).min(slice.len());
    std::str::from_utf8(&slice[start..stop]).unwrap()
}

fn main() {
    let path = std::env::args().nth(1).expect("No input file given");
    let mut reader = parse_fastx_file(&path).expect("Cannot open file");
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
    while let Some(r) = reader.next() {
        let record = r.expect("Invalid record");
        let _event = parser.next();
        let line = record.start_line_number();

        let (left, right) = (record.id(), parser.get_header());
        if let Some(pos) = check_mismatch(left, right) {
            eprintln!("Header mismatch line {line} pos {pos}");
            eprintln!("Needletail: \t{}", get_scope(left, pos));
            eprintln!("Helicase: \t{}", get_scope(right, pos));
            eprintln!("----------------");
            return;
        }

        let seq = &record.seq();
        for left in MATCH_N.split(seq).filter(|&seq| !seq.is_empty()) {
            let right = parser.get_dna_string();
            if let Some(pos) = check_mismatch(left, right) {
                eprintln!("Seq mismatch line {line} pos {pos}");
                eprintln!("Needletail: \t{}", get_scope(left, pos));
                eprintln!("Helicase: \t{}", get_scope(right, pos));
                eprintln!("----------------");
                return;
            }
            let _event = parser.next();
        }

        let (left, right) = (
            record.qual().unwrap_or(b""),
            parser.get_quality().unwrap_or(b""),
        );
        if let Some(pos) = check_mismatch(left, right) {
            eprintln!("Quality mismatch line {line} pos {pos}");
            eprintln!("Needletail: \t{}", get_scope(left, pos));
            eprintln!("Helicase: \t{}", get_scope(right, pos));
            eprintln!("----------------");
            return;
        }
    }
    eprintln!("Identical results ✓")
}
