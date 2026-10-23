// Regressions found by `tests/random_inputs.rs`.
use helicase::config::*;
use helicase::input::*;
use helicase::parser::{Error, Event};
use helicase::*;
use std::io::Cursor;

const KEEP: Config = ParserOptions::default().compute_quality().config();
const NO_QUALITY: Config = ParserOptions::default().config();
const SPLIT: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .config();
const SPLIT_RECORDS: Config = ParserOptions::default()
    .dna_string()
    .and_dna_packed()
    .return_record(true)
    .return_dna_chunk(false)
    .config();
const REJECT: Config = ParserOptions::default().reject_non_actg().config();

// Sequences of the events, or the error message, for a slice and a reader.
fn parse<const CONFIG: Config>(input: &[u8]) -> Vec<Result<String, String>> {
    fn collect<P: HelicaseParser + Iterator<Item = Result<Event, Error>>>(
        mut parser: P,
    ) -> Vec<Result<String, String>> {
        let mut out = Vec::new();
        while let Some(event) = parser.next() {
            assert!(out.len() < 100, "endless iteration");
            out.push(
                event
                    .map(|_| String::from_utf8(parser.get_dna_string().to_vec()).unwrap())
                    .map_err(|e| e.to_string()),
            );
        }
        out
    }
    let slice = collect(FastxParser::<CONFIG>::from_slice(input).unwrap());
    let reader = collect(FastxParser::<CONFIG>::from_reader(Cursor::new(input)).unwrap());
    assert_eq!(slice, reader, "slice and reader differ");
    slice
}

fn ok(seqs: &[&str]) -> Vec<Result<String, String>> {
    seqs.iter().map(|s| Ok(s.to_string())).collect()
}

#[test]
fn test_fastq_blank_lines_between_records() {
    let input = b"@a\nAC\n+\nII\n\n\n@b\nGT\n+\nII\n\n";
    assert_eq!(parse::<KEEP>(input), ok(&["AC", "GT"]));
    let input = b"@a\r\nAC\r\n+\r\nII\r\n\r\n@b\r\nGT\r\n+\r\nII\r\n\r\n";
    assert_eq!(parse::<KEEP>(input), ok(&["AC", "GT"]));
    let records = parse::<KEEP>(b"@a\nAC\n+\nII\n\nx\n@b\nGT\n+\nII\n");
    assert_eq!(records[0], Ok("AC".into()));
    assert!(records[1].as_ref().unwrap_err().contains("expected '@'"));
    assert_eq!(records.len(), 2);
}

#[test]
fn test_fastq_blank_line_with_header_on_plus_line() {
    // Shifted by the blank line, the plus line has the length of the next sequence line.
    // The line index is active after 16 short records.
    let record = |i: usize| format!("@r{i} x\n{}\n+r{i} x\nIIII\n", "ACGT");
    let mut input: String = (0..50).map(record).collect();
    input.insert(input.find("@r30").unwrap(), '\n');
    assert_eq!(parse::<KEEP>(input.as_bytes()), ok(&["ACGT"; 50]));
}

#[test]
fn test_split_after_non_actg_at_line_end_and_blank_line() {
    let input = b">a\nACGN\n\nCGA\n";
    assert_eq!(parse::<SPLIT>(input), ok(&["ACG", "CGA"]));
}

#[test]
fn test_fastq_split_truncated_record_ends_iteration() {
    let records = parse::<SPLIT>(b"@a\nACGT\n+");
    assert_eq!(records.len(), 2);
    assert!(records[1].is_err());
}

#[test]
fn test_crlf_without_final_line_feed() {
    assert_eq!(parse::<KEEP>(b">a\r\nACGT"), ok(&["ACGT"]));
    assert_eq!(parse::<KEEP>(b">a\r\nACGT\r\nAC"), ok(&["ACGTAC"]));
    assert_eq!(parse::<KEEP>(b"@a\r\nACGT\r\n+\r\nIIII"), ok(&["ACGT"]));
}

#[test]
fn test_fastq_truncated_quality_without_quality() {
    assert!(parse::<NO_QUALITY>(b"@a\nACGT\n+\nII")[0].is_err());
}

#[test]
fn test_fastq_empty_sequence_with_records() {
    let input = b"@a\nA\n+\nI\n@b\n\n+\n\n@c\nAC\n+\nII\n";
    assert_eq!(parse::<SPLIT_RECORDS>(input), ok(&["A", "", "AC"]));
    assert_eq!(parse::<REJECT>(input), ok(&["A", "", "AC"]));
}

#[test]
fn test_line_feed_without_carriage_return_in_crlf_input() {
    // The output is unspecified with mixed line endings, but must not panic.
    parse::<KEEP>(b">a\r\nAC\r\n>\nGT\r\n");
    parse::<KEEP>(b"@a\r\nAC\r\n+\r\nII\r\n@\n\n+\n\n");
}

#[test]
fn test_wrapped_record_ending_after_stride() {
    // Full lines of 63 bases up to the end of input, within the last block or at its end.
    let line = format!("{}AGC", "ACGT".repeat(15));
    for header in [
        "r",
        "r0123456789012345678901234567890123456789012345678901234567890",
    ] {
        let input = format!(">{header}\n{}", format!("{line}\n").repeat(4));
        assert_eq!(parse::<KEEP>(input.as_bytes()), ok(&[&line.repeat(4)]));
        assert_eq!(parse::<REJECT>(input.as_bytes()), ok(&[&line.repeat(4)]));
    }
}
