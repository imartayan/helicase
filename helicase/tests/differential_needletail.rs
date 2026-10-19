use fastx_gen::{Alphabet, FastxGen};
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use needletail::parse_fastx_reader;

const FASTQ_CONFIG: Config = ParserOptions::default().compute_quality().config();

fn run_differential_test(generator: FastxGen) {
    let raw_bytes = generator.generate_bytes();

    let mut nt_reader =
        parse_fastx_reader(raw_bytes.as_slice()).expect("Needletail failed to parse");
    let mut nt_records = Vec::new();
    while let Some(r) = nt_reader.next() {
        let rec = r.expect("Needletail error on valid record");
        nt_records.push((
            rec.id().to_vec(),
            rec.seq().to_vec(),
            rec.qual().map(|q| q.to_vec()),
        ));
    }

    // 1. Verify borrowed views against Needletail
    let mut parser_borrowed = FastxParser::<FASTQ_CONFIG>::from_slice(&raw_bytes)
        .expect("Helicase failed to initialize parser");
    let mut idx = 0;
    while let Some(_event) = parser_borrowed.next() {
        assert!(
            idx < nt_records.len(),
            "Helicase (borrowed) found more records than Needletail"
        );
        let (nt_hdr, nt_seq, nt_qual) = &nt_records[idx];
        assert_eq!(
            parser_borrowed.get_header(),
            nt_hdr.as_slice(),
            "Record {idx} borrowed header mismatch"
        );
        assert_eq!(
            parser_borrowed.get_dna_string(),
            nt_seq.as_slice(),
            "Record {idx} borrowed seq mismatch"
        );
        assert_eq!(
            parser_borrowed.get_quality(),
            nt_qual.as_deref(),
            "Record {idx} borrowed qual mismatch"
        );
        idx += 1;
    }
    assert_eq!(
        idx,
        nt_records.len(),
        "Helicase (borrowed) record count mismatch: Needletail found {}, Helicase found {idx}",
        nt_records.len()
    );

    // 2. Verify owned variants against Needletail
    let mut parser_owned = FastxParser::<FASTQ_CONFIG>::from_slice(&raw_bytes)
        .expect("Helicase failed to initialize parser");
    let mut idx = 0;
    while let Some(_event) = parser_owned.next() {
        assert!(
            idx < nt_records.len(),
            "Helicase (owned) found more records than Needletail"
        );
        let (nt_hdr, nt_seq, nt_qual) = &nt_records[idx];
        assert_eq!(
            &parser_owned.get_header_owned(),
            nt_hdr,
            "Record {idx} owned header mismatch"
        );
        assert_eq!(
            &parser_owned.get_dna_string_owned(),
            nt_seq,
            "Record {idx} owned seq mismatch"
        );
        assert_eq!(
            &parser_owned.get_quality_owned(),
            nt_qual,
            "Record {idx} owned qual mismatch"
        );
        idx += 1;
    }
    assert_eq!(
        idx,
        nt_records.len(),
        "Helicase (owned) record count mismatch: Needletail found {}, Helicase found {idx}",
        nt_records.len()
    );
}

#[test]
fn test_diff_fasta_short_unwrapped() {
    run_differential_test(FastxGen::fasta(50, 120).with_std_dev(0.15).with_seed(101));
}

#[test]
fn test_diff_fasta_protein() {
    run_differential_test(FastxGen::protein(20, 150).with_seed(999));
}

#[test]
fn test_diff_fasta_wrapped_50() {
    run_differential_test(
        FastxGen::wrapped_fasta(30, 250, 50)
            .with_alphabet(Alphabet::DnaWithN)
            .with_std_dev(0.1)
            .with_seed(101),
    );
}

#[test]
fn test_diff_fasta_wrapped_60() {
    run_differential_test(
        FastxGen::wrapped_fasta(30, 300, 60)
            .with_alphabet(Alphabet::DnaWithN)
            .with_std_dev(0.15)
            .with_seed(202),
    );
}

#[test]
fn test_diff_fasta_wrapped_80() {
    run_differential_test(
        FastxGen::wrapped_fasta(20, 500, 80)
            .with_std_dev(0.2)
            .with_seed(303),
    );
}

#[test]
fn test_diff_fasta_wrapped_120() {
    run_differential_test(
        FastxGen::wrapped_fasta(15, 600, 120)
            .with_std_dev(0.15)
            .with_seed(350),
    );
}

#[test]
fn test_diff_fasta_long_reads() {
    run_differential_test(
        FastxGen::wrapped_fasta(5, 10_000, 70)
            .with_std_dev(0.1)
            .with_seed(404),
    );
}

#[test]
fn test_diff_fastq_illumina_150bp() {
    run_differential_test(FastxGen::illumina(50).with_seed(505));
}

#[test]
fn test_diff_fastq_variable_len() {
    run_differential_test(
        FastxGen::illumina(40)
            .with_mean_len(250)
            .with_alphabet(Alphabet::DnaWithN)
            .with_std_dev(0.25)
            .with_seed(606),
    );
}

#[test]
fn test_diff_fastq_pacbio_long_reads() {
    run_differential_test(FastxGen::pacbio(5).with_seed(707));
}

#[test]
fn test_diff_fastq_at_in_qualities() {
    run_differential_test(
        FastxGen::illumina(100)
            .with_qual_range(b'@', b'@')
            .with_seed(808),
    );
}
