use fastx_gen::{Corrupt, FastxGen};
#[cfg(feature = "gz")]
use flate2::Compression;
#[cfg(feature = "gz")]
use flate2::write::GzEncoder;
use helicase::config::*;
use helicase::input::*;
use helicase::parser::Event;
use helicase::*;
#[cfg(feature = "gz")]
use std::io::Write;

const CONFIG: Config = ParserOptions::default().compute_quality().config();

#[test]
fn test_error_empty_slice() {
    // Empty input should fail initialization cleanly
    assert!(FastxParser::<CONFIG>::from_slice(b"").is_err());
    assert!(FastaParser::<CONFIG, _>::from_slice(b"").is_err());
    assert!(FastqParser::<CONFIG, _>::from_slice(b"").is_err());
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
        let res = FastxParser::<CONFIG>::from_slice(input);
        assert!(
            res.is_err(),
            "Expected error for input starting with byte {:?}, but got Ok",
            input.first()
        );
    }
}

#[test]
fn test_error_fasta_wrong_start_byte() {
    let res = FastaParser::<CONFIG, _>::from_slice(b"@fastq_header\nACGT\n");
    assert!(res.is_err());
}

#[test]
fn test_error_fastq_wrong_start_byte() {
    let res = FastqParser::<CONFIG, _>::from_slice(b">fasta_header\nACGT\n");
    assert!(res.is_err());
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

#[cfg(feature = "gz")]
#[test]
fn test_corrupted_gzip_trailer_crc() {
    let fastq_bytes = FastxGen::illumina(20)
        .with_mean_len(100)
        .with_seed(42)
        .generate_bytes();

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(&fastq_bytes).unwrap();
    let mut gz_bytes = encoder.finish().unwrap();

    // Corrupt one byte of the 4-byte CRC32 trailer (gzip trailer is 8 bytes: 4 bytes CRC32, 4 bytes ISIZE)
    let trailer_crc_pos = gz_bytes.len() - 8;
    gz_bytes[trailer_crc_pos] ^= 0xFF;

    // Test with ReaderInput via from_reader
    {
        let cursor = std::io::Cursor::new(gz_bytes.clone());
        let mut parser = FastxParser::<CONFIG>::from_reader(cursor).unwrap();
        let mut records = 0;
        let mut error_seen = false;

        for event in parser.by_ref() {
            match event {
                Event::Record(_) => records += 1,
                Event::IoError(err) => {
                    error_seen = true;
                    assert!(
                        err.to_string().to_lowercase().contains("checksum")
                            || err.to_string().to_lowercase().contains("corrupt"),
                        "Expected CRC checksum error, got: {err}"
                    );
                }
                _ => {}
            }
        }

        assert_eq!(
            records, 20,
            "Should yield all valid records before trailer error"
        );
        assert!(
            error_seen,
            "Should emit Event::IoError for corrupted trailer CRC"
        );
        assert!(
            parser.next().is_none(),
            "subsequent next() should return None"
        );
    }

    // Test with FileInput via from_file
    {
        let tmp_file = std::env::temp_dir().join(format!(
            "helicase_corrupted_crc_{}_{}.fastq.gz",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp_file, &gz_bytes).unwrap();

        let mut parser = FastxParser::<CONFIG>::from_file(&tmp_file).unwrap();
        let mut records = 0;
        let mut error_seen = false;

        for event in parser.by_ref() {
            match event {
                Event::Record(_) => records += 1,
                Event::IoError(err) => {
                    error_seen = true;
                    assert!(
                        err.to_string().to_lowercase().contains("checksum")
                            || err.to_string().to_lowercase().contains("corrupt"),
                        "Expected CRC checksum error, got: {err}"
                    );
                }
                _ => {}
            }
        }

        assert_eq!(
            records, 20,
            "Should yield all valid records before trailer error"
        );
        assert!(
            error_seen,
            "Should emit Event::IoError for corrupted trailer CRC"
        );
        assert!(
            parser.next().is_none(),
            "subsequent next() should return None"
        );

        let _ = std::fs::remove_file(&tmp_file);
    }
}

#[cfg(feature = "gz")]
#[test]
fn test_corrupted_gzip_trailer_crc_fasta() {
    let fasta_bytes = FastxGen::fasta(15, 120).with_seed(123).generate_bytes();

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(&fasta_bytes).unwrap();
    let mut gz_bytes = encoder.finish().unwrap();

    let trailer_crc_pos = gz_bytes.len() - 8;
    gz_bytes[trailer_crc_pos] ^= 0xAA;

    let cursor = std::io::Cursor::new(gz_bytes);
    let mut parser = FastxParser::<CONFIG>::from_reader(cursor).unwrap();
    let mut records = 0;
    let mut error_seen = false;

    for event in parser.by_ref() {
        match event {
            Event::Record(_) => records += 1,
            Event::IoError(_) => {
                error_seen = true;
            }
            _ => {}
        }
    }

    assert_eq!(records, 15);
    assert!(error_seen);
    assert!(parser.next().is_none());
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
