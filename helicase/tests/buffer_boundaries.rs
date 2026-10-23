#![cfg(feature = "std")]

use fastx_gen::{FastxGen, Format};
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use std::io::Cursor;

const FASTQ_CONFIG: Config = ParserOptions::default()
    .and_dna_columnar()
    .keep_non_actg()
    .compute_quality()
    .config();
const PACKED_KEEP_CONFIG: Config = ParserOptions::default()
    .dna_packed()
    .keep_non_actg()
    .config();

fn verify_record_lengths(format: Format, lengths: &[usize]) {
    for &target_len in lengths {
        let generator = match format {
            Format::Fasta => FastxGen::fasta(3, target_len),
            Format::Fastq => FastxGen::fastq(3, target_len),
        }
        .with_seed((target_len as u64).wrapping_mul(31));

        let bytes = generator.generate_bytes();
        let ground_truth: Vec<_> = generator.records().collect();

        // 1. Test via SliceInput
        let mut slice_parser = FastxParser::<FASTQ_CONFIG>::from_slice(&bytes).unwrap();
        let mut slice_packed = FastxParser::<PACKED_KEEP_CONFIG>::from_slice(&bytes).unwrap();
        let mut slice_records = Vec::new();
        while let Some(event) = slice_parser.next() {
            event.unwrap();
            slice_packed.next().unwrap().unwrap();
            slice_records.push((
                slice_parser.get_header_owned(),
                slice_parser.get_dna_string_owned(),
                slice_packed.get_dna_packed().to_string().into_bytes(),
                slice_parser.get_dna_columnar().to_string().into_bytes(),
                slice_parser.get_quality_owned(),
            ));
        }

        // 2. Test via ReaderInput (Cursor)
        let mut reader_parser =
            FastxParser::<FASTQ_CONFIG>::from_reader(Cursor::new(&bytes)).unwrap();
        let mut reader_packed =
            FastxParser::<PACKED_KEEP_CONFIG>::from_reader(Cursor::new(&bytes)).unwrap();
        let mut reader_records = Vec::new();
        while let Some(event) = reader_parser.next() {
            event.unwrap();
            reader_packed.next().unwrap().unwrap();
            reader_records.push((
                reader_parser.get_header_owned(),
                reader_parser.get_dna_string_owned(),
                reader_packed.get_dna_packed().to_string().into_bytes(),
                reader_parser.get_dna_columnar().to_string().into_bytes(),
                reader_parser.get_quality_owned(),
            ));
        }

        assert_eq!(
            slice_records.len(),
            ground_truth.len(),
            "Slice record count mismatch for length {target_len}"
        );
        assert_eq!(
            reader_records.len(),
            ground_truth.len(),
            "Reader record count mismatch for length {target_len}"
        );

        for (i, expected) in ground_truth.iter().enumerate() {
            assert_eq!(
                slice_records[i].1.len(),
                target_len,
                "Slice sequence length mismatch at record {i} for target {target_len}"
            );
            assert_eq!(
                slice_records[i].1, expected.seq,
                "Slice sequence content mismatch at record {i} for target {target_len}"
            );
            assert_eq!(
                slice_records[i].2, expected.seq,
                "Slice packed DNA content mismatch at record {i} for target {target_len}"
            );
            assert_eq!(
                reader_records[i], slice_records[i],
                "Reader vs Slice record mismatch at record {i} for target {target_len}"
            );
        }
    }
}

#[test]
fn test_boundary_lengths_around_64_bytes() {
    // Exact block boundaries: 64-byte SIMD vector boundary
    let lengths = [1, 2, 62, 63, 64, 65, 66, 126, 127, 128, 129, 130];
    verify_record_lengths(Format::Fasta, &lengths);
    verify_record_lengths(Format::Fastq, &lengths);
}

#[test]
fn test_boundary_lengths_around_1024_bytes() {
    // 1024 bytes is the long-DNA & FASTQ quality fast-path threshold
    let lengths = [1021, 1022, 1023, 1024, 1025, 1026];
    verify_record_lengths(Format::Fasta, &lengths);
    verify_record_lengths(Format::Fastq, &lengths);
}

#[test]
fn test_boundary_lengths_around_buffer_size() {
    for buffer_size in [1 << 17, 1 << 19] {
        let lengths = [
            buffer_size - 65,
            buffer_size - 1,
            buffer_size,
            buffer_size + 1,
            buffer_size + 65,
        ];
        verify_record_lengths(Format::Fasta, &lengths);
        verify_record_lengths(Format::Fastq, &lengths);
    }
}
