#![cfg(feature = "std")]

use fastx_gen::{Corrupt, FastxGen};
#[cfg(feature = "gz")]
use flate2::Compression;
#[cfg(feature = "gz")]
use flate2::write::GzEncoder;
use helicase::config::*;
use helicase::input::*;
use helicase::parser::{Error, ErrorKind, Event};
use helicase::*;
#[cfg(feature = "gz")]
use std::io::Write;

const CONFIG: Config = ParserOptions::default().compute_quality().config();

fn count_events<P: Iterator<Item = Result<Event, Error>>>(parser: P) -> usize {
    parser.collect::<Result<Vec<_>, _>>().unwrap().len()
}

#[test]
fn test_empty_input_yields_no_records() {
    const PACKED: Config = ParserOptions::default().dna_packed().config();
    let empty: &[u8] = b"";
    assert_eq!(
        count_events(FastxParser::<CONFIG>::from_slice(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(FastxParser::<PACKED>::from_slice(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(FastxParser::<CONFIG>::from_reader(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(FastxParser::<PACKED>::from_reader(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(parser::FastaParser::<CONFIG, _>::from_slice(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(parser::FastqParser::<CONFIG, _>::from_slice(empty).unwrap()),
        0
    );
    assert_eq!(
        count_events(parser::FastqParser::<PACKED, _>::from_reader(empty).unwrap()),
        0
    );
}

#[test]
fn test_error_invalid_first_byte() {
    let invalid_inputs: [&[u8]; 4] = [
        b"not a fasta or fastq",
        b"\n>valid_later\nACGT\n",
        b" 12345\n",
        b"!record_0\nACGT\n+\nIIII\n",
    ];

    for input in invalid_inputs {
        match FastxParser::<CONFIG>::from_slice(input) {
            Err(e) => assert!(matches!(
                e.kind(),
                ErrorKind::InvalidRecordStart {
                    expected: _,
                    found: _
                }
            )),
            Ok(_) => panic!(
                "Expected error for input starting with byte {:?}, but got Ok",
                input.first()
            ),
        }
    }
}

#[test]
fn test_error_fasta_wrong_start_byte() {
    match parser::FastaParser::<CONFIG, _>::from_slice(b"@fastq_header\nACGT\n") {
        Err(e) => assert!(matches!(
            e.kind(),
            ErrorKind::InvalidRecordStart {
                expected: "'>'",
                found: b'@'
            }
        )),
        Ok(_) => panic!("Expected error for wrong start byte"),
    }
}

#[test]
fn test_error_fastq_wrong_start_byte() {
    match parser::FastqParser::<CONFIG, _>::from_slice(b">fasta_header\nACGT\n") {
        Err(e) => assert!(matches!(
            e.kind(),
            ErrorKind::InvalidRecordStart {
                expected: "'@'",
                found: b'>'
            }
        )),
        Ok(_) => panic!("Expected error for wrong start byte"),
    }
}

#[test]
fn test_error_size_is_pointer_sized() {
    assert_eq!(std::mem::size_of::<Error>(), std::mem::size_of::<usize>());
    assert_eq!(
        std::mem::size_of::<Result<Event, Error>>(),
        std::mem::size_of::<Event>()
    );
}

#[test]
fn test_truncated_record_error() {
    // Truncated in header (no trailing newline)
    let raw = b"@header_no_newline";
    let mut parser = parser::FastqParser::<CONFIG, _>::from_slice(raw).unwrap();
    let res = parser.next().unwrap();
    assert!(res.is_err());
    assert!(matches!(
        res.unwrap_err().kind(),
        ErrorKind::TruncatedRecord
    ));
    assert!(parser.next().is_none());

    // Truncated in sequence (seq ends without newline or plus)
    let raw = b"@h\nACGT";
    let mut parser = parser::FastqParser::<CONFIG, _>::from_slice(raw).unwrap();
    let res = parser.next().unwrap();
    assert!(res.is_err());
    assert!(matches!(
        res.unwrap_err().kind(),
        ErrorKind::TruncatedRecord
    ));
    assert!(parser.next().is_none());

    // Truncated in plus line
    let raw = b"@h\nACGT\n+";
    let mut parser = parser::FastqParser::<CONFIG, _>::from_slice(raw).unwrap();
    let res = parser.next().unwrap();
    assert!(res.is_err());
    assert!(matches!(
        res.unwrap_err().kind(),
        ErrorKind::TruncatedRecord
    ));
    assert!(parser.next().is_none());

    // Truncated in quality (seq has 4 bp, but qual only has 2 bytes before EOF)
    let raw = b"@h\nACGT\n+\nII";
    let mut parser = parser::FastqParser::<CONFIG, _>::from_slice(raw).unwrap();
    let res = parser.next().unwrap();
    assert!(res.is_err());
    assert!(matches!(
        res.unwrap_err().kind(),
        ErrorKind::TruncatedRecord
    ));
    assert!(parser.next().is_none());
}

#[test]
fn test_truncated_mid_record_terminates_gracefully() {
    // Generate valid records then truncate midway through the last record
    let bytes = FastxGen::illumina(5)
        .with_mean_len(100)
        .with_seed(111)
        .generate_bytes();

    // Truncate at 70% of the file length
    let truncated = &bytes[..bytes.len() * 7 / 10];
    let parser_res = FastxParser::<CONFIG>::from_slice(truncated);

    if let Ok(mut parser) = parser_res {
        // Parsing should terminate with None and not panic or loop infinitely
        let mut count = 0;
        for _event in parser.by_ref() {
            count += 1;
            assert!(
                count <= 5,
                "Parser should not yield more records than exist"
            );
        }
    }
}

#[test]
fn test_corrupt_generator_kind_invalid_header_prefix() {
    let bytes = FastxGen::illumina(1)
        .with_mean_len(50)
        .with_corrupt(Corrupt::InvalidHeaderPrefix(1.0))
        .with_seed(222)
        .generate_bytes();
    assert_eq!(bytes[0], b'!');
    assert!(FastxParser::<CONFIG>::from_slice(&bytes).is_err());
}

#[test]
fn test_empty_sequence_fasta_blank_line() {
    let raw = b">empty1\n\n>empty2\n\n>valid\nACGT\n";
    let mut parser = FastxParser::<CONFIG>::from_slice(raw).unwrap();
    let mut headers = Vec::new();
    let mut seqs = Vec::new();
    while let Some(_event) = parser.next() {
        headers.push(parser.get_header_owned());
        seqs.push(parser.get_dna_string_owned());
    }
    assert_eq!(headers.len(), 3);
    assert_eq!(headers[0], b"empty1");
    assert_eq!(seqs[0], b"");
    assert_eq!(headers[1], b"empty2");
    assert_eq!(seqs[1], b"");
    assert_eq!(headers[2], b"valid");
    assert_eq!(seqs[2], b"ACGT");
}

#[test]
fn test_empty_sequence_fasta_trailing_header() {
    let raw = b">valid\nACGT\n>trailing_empty\n";
    let mut parser = FastxParser::<CONFIG>::from_slice(raw).unwrap();
    let mut headers = Vec::new();
    let mut seqs = Vec::new();
    while let Some(_event) = parser.next() {
        headers.push(parser.get_header_owned());
        seqs.push(parser.get_dna_string_owned());
    }
    assert_eq!(headers.len(), 2);
    assert_eq!(headers[0], b"valid");
    assert_eq!(seqs[0], b"ACGT");
    assert_eq!(headers[1], b"trailing_empty");
    assert_eq!(seqs[1], b"");
}

#[test]
fn test_empty_sequence_fastq() {
    let raw = b"@empty\n\n+\n\n@valid\nACGT\n+\nIIII\n";
    let mut parser = FastxParser::<CONFIG>::from_slice(raw).unwrap();
    let mut records = Vec::new();
    while let Some(_event) = parser.next() {
        records.push((
            parser.get_header_owned(),
            parser.get_dna_string_owned(),
            parser.get_quality_owned(),
        ));
    }
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, b"empty");
    assert_eq!(records[0].1, b"");
    assert_eq!(records[0].2, Some(b"".to_vec()));
    assert_eq!(records[1].0, b"valid");
    assert_eq!(records[1].1, b"ACGT");
    assert_eq!(records[1].2, Some(b"IIII".to_vec()));
}

// The error may surface when building the parser, since it reads the first block
#[cfg(feature = "gz")]
fn assert_checksum_error(parser: Result<FastxParser<CONFIG>, Error>) {
    let err = match parser {
        Ok(mut parser) => {
            let err = parser
                .by_ref()
                .find_map(Result::err)
                .expect("missing error");
            assert!(
                parser.next().is_none(),
                "subsequent next() should return None"
            );
            err
        }
        Err(err) => err,
    };
    let msg = err.to_string().to_lowercase();
    assert!(msg.contains("checksum") || msg.contains("corrupt"), "{msg}");
}

#[cfg(feature = "gz")]
#[test]
fn test_corrupted_gzip_trailer_crc() {
    let fastq = FastxGen::illumina(20).with_mean_len(100).with_seed(42);
    for bytes in [
        fastq.generate_bytes(),
        FastxGen::fasta(15, 120).with_seed(123).generate_bytes(),
    ] {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&bytes).unwrap();
        let mut gz_bytes = encoder.finish().unwrap();
        // The gzip trailer is a 4-byte CRC32 then a 4-byte size
        let crc_pos = gz_bytes.len() - 8;
        gz_bytes[crc_pos] ^= 0xFF;

        assert_checksum_error(FastxParser::from_reader(std::io::Cursor::new(
            gz_bytes.clone(),
        )));
        let tmp_file = std::env::temp_dir().join(format!(
            "helicase_corrupted_crc_{}_{}.gz",
            std::process::id(),
            bytes[0]
        ));
        std::fs::write(&tmp_file, &gz_bytes).unwrap();
        assert_checksum_error(FastxParser::from_file(&tmp_file));
        let _ = std::fs::remove_file(&tmp_file);
    }
}

#[test]
fn test_reader_input_fused_on_error() {
    struct FailingReader {
        data: Vec<u8>,
        pos: usize,
    }
    impl std::io::Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.pos < self.data.len() {
                let n = (self.data.len() - self.pos).min(buf.len());
                buf[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
                self.pos += n;
                Ok(n)
            } else {
                Err(std::io::Error::other("simulated I/O failure"))
            }
        }
    }

    let mut input = ReaderInput::new(FailingReader {
        data: b">test\nACGT\n".to_vec(),
        pos: 0,
    })
    .unwrap();
    // First block: Ok
    let first = input.next();
    assert!(first.is_some() && first.unwrap().is_ok());

    // Next block: Err
    let err = input.next();
    assert!(err.is_some() && err.unwrap().is_err());

    // Fused: all subsequent calls must be None, without attempting read again
    assert!(input.next().is_none());
    assert!(input.next().is_none());
}
