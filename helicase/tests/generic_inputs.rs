#![cfg(feature = "std")]

use fastx_gen::{Alphabet, FastxGen};
use flate2::Compression;
use flate2::write::GzEncoder;
use helicase::config::*;
use helicase::input::*;
use helicase::parser::Event;
use helicase::*;
use std::fs::File;
use std::io::{Cursor, Write};

const CONFIG: Config = ParserOptions::default()
    .and_dna_packed()
    .and_dna_columnar()
    .keep_non_actg()
    .compute_quality()
    .config();

#[derive(Clone, Debug, PartialEq, Eq)]
struct ParsedRecord {
    header: Vec<u8>,
    dna: Vec<u8>,
    packed: Vec<u8>,
    columnar: Vec<u8>,
    qual: Option<Vec<u8>>,
}

fn parse_all<P: HelicaseParser + Iterator<Item = parser::Result<Event>>>(
    mut parser: P,
) -> Vec<ParsedRecord> {
    let mut records = Vec::new();
    while let Some(event) = parser.next() {
        event.unwrap();
        records.push(ParsedRecord {
            header: parser.get_header_owned(),
            dna: parser.get_dna_string_owned(),
            packed: parser.get_dna_packed().to_string().into_bytes(),
            columnar: parser.get_dna_columnar().to_string().into_bytes(),
            qual: parser.get_quality_owned(),
        });
    }
    records
}

fn test_across_all_inputs(generator: FastxGen, test_id: &str) {
    let raw_bytes = generator.generate_bytes();

    // 1. SliceInput
    let from_slice = parse_all(FastxParser::<CONFIG>::from_slice(&raw_bytes).unwrap());

    // 2. ReaderInput (Cursor)
    let from_reader =
        parse_all(FastxParser::<CONFIG>::from_reader(Cursor::new(raw_bytes.clone())).unwrap());
    assert_eq!(
        from_slice, from_reader,
        "ReaderInput mismatch for {test_id}"
    );

    // Setup temporary files on disk
    let tmp_dir =
        std::env::temp_dir().join(format!("helicase_test_{}_{}", test_id, std::process::id()));
    std::fs::create_dir_all(&tmp_dir).unwrap();

    let raw_file_path = tmp_dir.join("input.raw");
    let gz_file_path = tmp_dir.join("input.raw.gz");

    // Write raw file
    {
        let mut f = File::create(&raw_file_path).unwrap();
        f.write_all(&raw_bytes).unwrap();
    }

    // Write gzip-compressed file
    {
        let f = File::create(&gz_file_path).unwrap();
        let mut encoder = GzEncoder::new(f, Compression::default());
        encoder.write_all(&raw_bytes).unwrap();
        encoder.finish().unwrap();
    }

    // 3. MmapInput
    let from_mmap = parse_all(FastxParser::<CONFIG>::from_file_mmap(&raw_file_path).unwrap());
    assert_eq!(from_slice, from_mmap, "MmapInput mismatch for {test_id}");

    // 5. FileInput (uncompressed)
    let from_file = parse_all(FastxParser::<CONFIG>::from_file(&raw_file_path).unwrap());
    assert_eq!(
        from_slice, from_file,
        "FileInput (uncompressed) mismatch for {test_id}"
    );

    // 6. FileInput (gzip compressed)
    if cfg!(feature = "gz") {
        let from_gz_file = parse_all(FastxParser::<CONFIG>::from_file(&gz_file_path).unwrap());
        assert_eq!(
            from_slice, from_gz_file,
            "FileInput (gzip) mismatch for {test_id}"
        );
    }

    // Clean up
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

#[test]
fn test_all_inputs_fastq_short() {
    test_across_all_inputs(FastxGen::illumina(30).with_seed(1234), "fastq_short");
}

#[test]
fn test_all_inputs_fastq_buffer_spill() {
    // Large dataset that forces ReaderInput buffer refill (>512KB buffer)
    let generator = FastxGen::illumina(1500)
        .with_mean_len(300)
        .with_std_dev(0.15)
        .with_seed(5678);
    test_across_all_inputs(generator, "fastq_buffer_spill");
}

#[test]
fn test_all_inputs_fasta_unwrapped() {
    test_across_all_inputs(FastxGen::fasta(50, 200).with_seed(9012), "fasta_unwrapped");
}

#[test]
fn test_all_inputs_fasta_wrapped() {
    test_across_all_inputs(
        FastxGen::wrapped_fasta(10, 2000, 60).with_seed(3456),
        "fasta_wrapped",
    );
}

#[test]
fn test_all_inputs_fasta_buffer_spill() {
    // Multi-line FASTA spanning multiple buffer boundaries (>128KB buffer)
    let generator = FastxGen::wrapped_fasta(500, 1500, 70)
        .with_alphabet(Alphabet::DnaWithN)
        .with_std_dev(0.15)
        .with_seed(7890);
    test_across_all_inputs(generator, "fasta_buffer_spill");
}

#[test]
fn test_all_inputs_mega_reads_fasta() {
    // Single unwrapped and wrapped records spanning multiple 128KB buffer refills
    for (generator, id) in [
        (
            FastxGen::fasta(2, 600_000)
                .with_std_dev(0.1)
                .with_seed(1122),
            "mega_reads_fasta",
        ),
        (
            FastxGen::wrapped_fasta(2, 400_000, 60)
                .with_std_dev(0.1)
                .with_seed(2233),
            "mega_reads_wrapped_fasta",
        ),
    ] {
        let raw_bytes = generator.clone().generate_bytes();
        test_across_all_inputs(generator, id);

        // Also test packed-only config (keep_anchor = false, constant 128KB buffer)
        const PACKED_CFG: Config = ParserOptions::default().dna_packed().config();
        let mut p_slc = FastxParser::<PACKED_CFG>::from_slice(&raw_bytes).unwrap();
        let mut p_rdr = FastxParser::<PACKED_CFG>::from_reader(Cursor::new(&raw_bytes)).unwrap();
        while let (Some(_), Some(_)) = (p_slc.next(), p_rdr.next()) {
            assert_eq!(p_slc.get_header(), p_rdr.get_header());
            assert_eq!(
                p_slc.get_dna_packed().to_string(),
                p_rdr.get_dna_packed().to_string()
            );
        }
        assert!(p_slc.next().is_none() && p_rdr.next().is_none());
    }
}

#[test]
fn test_reader_long_reads_without_headers() {
    // Without headers the reader keeps no anchor, long sequences resume at the start of a buffer.
    const KEEP_CFG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .keep_non_actg()
        .config();
    let raw_bytes = FastxGen::pacbio(100).with_seed(42).generate_bytes();
    let mut p_slc = FastxParser::<KEEP_CFG>::from_slice(&raw_bytes).unwrap();
    let mut p_rdr = FastxParser::<KEEP_CFG>::from_reader(Cursor::new(&raw_bytes)).unwrap();
    while let (Some(e_slc), Some(e_rdr)) = (p_slc.next(), p_rdr.next()) {
        e_slc.unwrap();
        e_rdr.unwrap();
        assert_eq!(
            p_slc.get_dna_packed().to_string(),
            p_rdr.get_dna_packed().to_string()
        );
    }
    assert!(p_slc.next().is_none() && p_rdr.next().is_none());
}

#[test]
fn test_all_inputs_mega_reads_fastq() {
    // Single FASTQ records spanning multiple 128KB buffer refills
    let generator = FastxGen::fastq(2, 600_000)
        .with_std_dev(0.1)
        .with_seed(3344);
    let raw_bytes = generator.generate_bytes();
    test_across_all_inputs(generator, "mega_reads_fastq");

    // Also test packed-only + ignore_quality (keep_anchor = false, O(1) multi-buffer quality skip)
    const PACKED_CFG: Config = ParserOptions::default().dna_packed().config();
    let mut p_slc = FastxParser::<PACKED_CFG>::from_slice(&raw_bytes).unwrap();
    let mut p_rdr = FastxParser::<PACKED_CFG>::from_reader(Cursor::new(&raw_bytes)).unwrap();
    while let (Some(_), Some(_)) = (p_slc.next(), p_rdr.next()) {
        assert_eq!(p_slc.get_header(), p_rdr.get_header());
        assert_eq!(
            p_slc.get_dna_packed().to_string(),
            p_rdr.get_dna_packed().to_string()
        );
    }
    assert!(p_slc.next().is_none() && p_rdr.next().is_none());
}

#[test]
fn test_borrowed_vs_owned_equivalence() {
    let fastq_bytes = FastxGen::illumina(50).with_seed(123).generate_bytes();
    let fasta_bytes = FastxGen::wrapped_fasta(30, 250, 60)
        .with_seed(456)
        .generate_bytes();

    for bytes in [&fastq_bytes, &fasta_bytes] {
        let mut p_borrow = FastxParser::<CONFIG>::from_slice(bytes).unwrap();
        let mut recs_borrow = Vec::new();
        while let Some(_) = p_borrow.next() {
            recs_borrow.push(ParsedRecord {
                header: p_borrow.get_header().to_vec(),
                dna: p_borrow.get_dna_string().to_vec(),
                packed: p_borrow.get_dna_packed().to_string().into_bytes(),
                columnar: p_borrow.get_dna_columnar().to_string().into_bytes(),
                qual: p_borrow.get_quality().map(|q| q.to_vec()),
            });
        }

        let p_owned = FastxParser::<CONFIG>::from_slice(bytes).unwrap();
        let recs_owned = parse_all(p_owned);

        assert_eq!(recs_borrow, recs_owned);
    }
}

#[test]
fn test_fastx_parser_is_send() {
    fn assert_send<T: Send>(_: &T) {}
    let data = b">r\nACGT\n";
    assert_send(&FastxParser::<CONFIG>::from_slice(data).unwrap());
    assert_send(&FastxParser::<CONFIG>::from_reader(Cursor::new(data)).unwrap());
    let parser = FastxParser::<CONFIG>::from_slice(data).unwrap();
    let count = std::thread::spawn(move || parser.count()).join().unwrap();
    assert_eq!(count, 1);
}

// `from_file` reads ahead in a background thread, `from_reader` reads on the parser thread.
#[test]
fn test_read_ahead_matches_direct_reader() {
    let tmp_dir = std::env::temp_dir().join(format!("helicase_read_ahead_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).unwrap();
    let datasets = [
        ("short", FastxGen::illumina(20_000).with_seed(1)),
        (
            "wrapped",
            FastxGen::wrapped_fasta(200, 30_000, 60).with_seed(2),
        ),
        ("long", FastxGen::fasta(40, 200_000).with_seed(3)),
    ];
    for (name, generator) in datasets {
        let raw = generator.generate_bytes();
        let mut gz = GzEncoder::new(Vec::new(), Compression::fast());
        gz.write_all(&raw).unwrap();
        let gz = gz.finish().unwrap();
        let truncated_gz = gz[..gz.len() * 2 / 3].to_vec();
        for (kind, bytes) in [("raw", raw), ("gz", gz), ("truncated gz", truncated_gz)] {
            if kind != "raw" && !cfg!(feature = "gz") {
                continue;
            }
            let path = tmp_dir.join(format!("{name}.{}", kind.replace(' ', "_")));
            std::fs::write(&path, &bytes).unwrap();
            let read_ahead = parse_until_error(FastxParser::<CONFIG>::from_file(&path).unwrap());
            let direct = parse_until_error(
                FastxParser::<CONFIG>::from_reader(File::open(&path).unwrap()).unwrap(),
            );
            assert!(read_ahead == direct, "{name} {kind}");
        }
    }
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

fn parse_until_error<P: HelicaseParser + Iterator<Item = parser::Result<Event>>>(
    mut parser: P,
) -> (Vec<ParsedRecord>, Option<String>) {
    let mut records = Vec::new();
    while let Some(event) = parser.next() {
        if let Err(e) = event {
            return (records, Some(e.to_string()));
        }
        records.push(ParsedRecord {
            header: parser.get_header_owned(),
            dna: parser.get_dna_string_owned(),
            packed: parser.get_dna_packed().to_string().into_bytes(),
            columnar: parser.get_dna_columnar().to_string().into_bytes(),
            qual: parser.get_quality_owned(),
        });
    }
    (records, None)
}
