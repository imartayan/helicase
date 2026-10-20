#![cfg(feature = "needletail")]

use helicase::*;

#[test]
fn test_needletail_method_aliases_fastq() {
    let data = b"@read1 some description\nACGTACGT\n+\nIIIIIIII\n";
    let mut parser = parse_fastx_slice(data).expect("failed to parse");

    assert!(parser.next().is_some());
    // Test Needletail aliases
    assert_eq!(parser.id(), b"read1 some description");
    assert_eq!(&*parser.seq(), b"ACGTACGT");
    assert_eq!(parser.raw_seq(), b"ACGTACGT");
    assert_eq!(parser.qual(), Some(&b"IIIIIIII"[..]));
    assert_eq!(parser.num_bases(), 8);

    // Verify aliases match canonical methods
    assert_eq!(parser.id(), parser.get_header());
    assert_eq!(parser.raw_seq(), parser.get_dna_string());
    assert_eq!(parser.qual(), parser.get_quality());
    assert_eq!(parser.num_bases(), parser.get_dna_len());

    assert!(parser.next().is_none());
}

#[test]
fn test_needletail_method_aliases_fasta() {
    let data = b">contig_1\nACGT\nACGT\n";
    let mut parser = parse_fastx_slice(data).expect("failed to parse");

    assert!(parser.next().is_some());
    assert_eq!(parser.id(), b"contig_1");
    assert_eq!(&*parser.seq(), b"ACGTACGT");
    assert_eq!(parser.raw_seq(), b"ACGTACGT");
    assert_eq!(parser.qual(), None);
    assert_eq!(parser.num_bases(), 8);

    assert!(parser.next().is_none());
}

#[test]
fn test_needletail_parse_reader() {
    let data = b"@read1\nACGT\n+\nIIII\n";
    let cursor = std::io::Cursor::new(data);
    let mut parser = parse_fastx_reader(cursor).expect("failed to parse reader");

    assert!(parser.next().is_some());
    assert_eq!(parser.id(), b"read1");
    assert_eq!(&*parser.seq(), b"ACGT");
    assert_eq!(parser.num_bases(), 4);
    assert!(parser.next().is_none());
}

#[test]
fn test_needletail_parse_file() {
    use std::io::Write;

    let path = std::env::temp_dir().join("test_needletail_helpers.fq");
    let mut f = std::fs::File::create(&path).expect("failed to create temp file");
    f.write_all(b"@read1\nACGTACGT\n+\nIIIIIIII\n")
        .expect("failed to write temp file");
    drop(f);

    let mut parser = parse_fastx_file(&path).expect("failed to parse file");
    assert!(parser.next().is_some());
    assert_eq!(parser.id(), b"read1");
    assert_eq!(&*parser.seq(), b"ACGTACGT");
    assert_eq!(parser.qual(), Some(&b"IIIIIIII"[..]));
    assert_eq!(parser.num_bases(), 8);
    assert!(parser.next().is_none());

    let _ = std::fs::remove_file(&path);
}
