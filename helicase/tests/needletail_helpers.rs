#![cfg(feature = "needletail")]

use helicase::parser::Format;
use helicase::*;

#[test]
fn test_needletail_record_fastq() {
    let data = b"@read1 some description\nACGTACGT\n+\nIIIIIIII\n";
    let mut reader = parse_fastx_slice(data).expect("failed to parse");

    let rec = reader.next().unwrap().unwrap();
    assert_eq!(rec.id(), b"read1 some description");
    assert_eq!(&*rec.seq(), b"ACGTACGT");
    assert_eq!(rec.raw_seq(), b"ACGTACGT");
    assert_eq!(rec.qual(), Some(&b"IIIIIIII"[..]));
    assert_eq!(rec.num_bases(), 8);
    assert_eq!(rec.format(), Format::Fastq);

    assert!(reader.next().is_none());
}

#[test]
fn test_needletail_record_fasta() {
    let data = b">contig_1\nACGT\nACGT\n>contig_2\nTT\n";
    let mut reader = parse_fastx_slice(data).expect("failed to parse");

    let rec = reader.next().unwrap().unwrap();
    assert_eq!(rec.id(), b"contig_1");
    assert_eq!(&*rec.seq(), b"ACGTACGT");
    assert_eq!(rec.raw_seq(), b"ACGTACGT");
    assert_eq!(rec.qual(), None);
    assert_eq!(rec.num_bases(), 8);
    assert_eq!(rec.format(), Format::Fasta);

    let rec = reader.next().unwrap().unwrap();
    assert_eq!(rec.id(), b"contig_2");
    assert_eq!(&*rec.seq(), b"TT");

    assert!(reader.next().is_none());
}

#[test]
fn test_needletail_parse_reader() {
    let data = b"@read1\nACGT\n+\nIIII\n";
    let cursor = std::io::Cursor::new(data);
    let mut reader = parse_fastx_reader(cursor).expect("failed to parse reader");

    let rec = reader.next().unwrap().unwrap();
    assert_eq!(rec.id(), b"read1");
    assert_eq!(&*rec.seq(), b"ACGT");
    assert_eq!(rec.num_bases(), 4);
    assert!(reader.next().is_none());
}

#[test]
fn test_needletail_parse_file() {
    use std::io::Write;

    let path = std::env::temp_dir().join("test_needletail_helpers.fq");
    let mut f = std::fs::File::create(&path).expect("failed to create temp file");
    f.write_all(b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nAC\n+\nII\n")
        .expect("failed to write temp file");
    drop(f);

    let mut reader = parse_fastx_file(&path).expect("failed to parse file");
    let mut ids = Vec::new();
    let mut total_bases = 0;
    while let Some(record) = reader.next() {
        let rec = record.expect("invalid record");
        ids.push(rec.id().to_vec());
        total_bases += rec.num_bases();
    }
    assert_eq!(ids, [b"read1".to_vec(), b"read2".to_vec()]);
    assert_eq!(total_bases, 10);

    let _ = std::fs::remove_file(&path);
}

#[test]
fn test_needletail_truncated_fastq_errors() {
    let mut reader = parse_fastx_slice(b"@read1\nACGT\n+\n").expect("failed to parse");
    assert!(reader.next().unwrap().is_err());
}
