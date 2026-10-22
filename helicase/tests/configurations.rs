use fastx_gen::{Alphabet, FastxGen};
use helicase::config::*;
use helicase::input::*;
use helicase::parser::Event;
use helicase::*;

#[test]
fn test_config_packed_and_columnar_against_string() {
    // Generate valid DNA records with varied lengths
    let bytes = FastxGen::illumina(25)
        .with_mean_len(200)
        .with_std_dev(0.2)
        .with_seed(4242)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .and_dna_packed()
        .and_dna_columnar()
        .keep_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string();
        let p = parser.get_dna_packed();
        let c = parser.get_dna_columnar();

        assert_eq!(p.len(), s.len(), "Packed len mismatch");
        assert_eq!(c.len(), s.len(), "Columnar len mismatch");

        for (i, &b) in s.iter().enumerate() {
            let expected = b.to_ascii_uppercase() as char;
            assert_eq!(p.get_char(i), expected, "Packed base mismatch at index {i}");
            assert_eq!(
                c.get_char(i),
                expected,
                "Columnar base mismatch at index {i}"
            );
        }
    }
}

#[test]
fn test_config_packed_and_columnar_wrapped_fasta() {
    for &width in &[50, 60, 80, 150] {
        let bytes = FastxGen::wrapped_fasta(10, 350, width)
            .with_std_dev(0.15)
            .with_seed(1234 + width as u64)
            .generate_bytes();

        const CFG: Config = ParserOptions::default()
            .dna_string()
            .and_dna_packed()
            .and_dna_columnar()
            .keep_non_actg()
            .config();

        let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
        while let Some(_event) = parser.next() {
            let s = parser.get_dna_string();
            let p = parser.get_dna_packed();
            let c = parser.get_dna_columnar();

            assert_eq!(
                p.len(),
                s.len(),
                "Packed len mismatch in wrapped FASTA (width {width})"
            );
            assert_eq!(
                c.len(),
                s.len(),
                "Columnar len mismatch in wrapped FASTA (width {width})"
            );

            for (i, &b) in s.iter().enumerate() {
                let expected = b.to_ascii_uppercase() as char;
                assert_eq!(
                    p.get_char(i),
                    expected,
                    "Packed base mismatch at index {i} (width {width})"
                );
                assert_eq!(
                    c.get_char(i),
                    expected,
                    "Columnar base mismatch at index {i} (width {width})"
                );
            }
        }
    }
}

#[test]
fn test_config_packed_and_columnar_owned_against_string() {
    let bytes = FastxGen::illumina(25)
        .with_mean_len(200)
        .with_std_dev(0.2)
        .with_seed(4242)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .and_dna_packed()
        .and_dna_columnar()
        .keep_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string_owned();
        let p = parser.get_dna_packed_owned();
        let c = parser.get_dna_columnar_owned();

        assert_eq!(p.len(), s.len(), "Packed len mismatch");
        assert_eq!(c.len(), s.len(), "Columnar len mismatch");

        for (i, &b) in s.iter().enumerate() {
            let expected = b.to_ascii_uppercase() as char;
            assert_eq!(p.get_char(i), expected, "Packed base mismatch at index {i}");
            assert_eq!(
                c.get_char(i),
                expected,
                "Columnar base mismatch at index {i}"
            );
        }
    }
}

#[test]
fn test_config_split_non_actg_yields_pure_actg_chunks() {
    // Generate records with IUPAC bases (containing non-ACTG characters)
    let bytes = FastxGen::fasta(15, 150)
        .with_alphabet(Alphabet::DnaIupac)
        .with_seed(1337)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .split_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    let mut chunk_count = 0;
    while let Some(Ok(event)) = parser.next() {
        if let Event::DnaChunk(_) = event {
            chunk_count += 1;
            let chunk = parser.get_dna_string();
            assert!(!chunk.is_empty(), "DnaChunk should not be empty");
            for &b in chunk {
                assert!(
                    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'),
                    "DnaChunk contained non-ACTG base: '{}'",
                    b as char
                );
            }
        }
    }
    assert!(chunk_count > 0, "Should have yielded DnaChunk events");
}

#[test]
fn test_config_split_non_actg_long_runs() {
    let mut state = 42u64;
    let mut seq: Vec<u8> = (0..600_000)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGTacgt"[(state >> 61) as usize]
        })
        .collect();
    // Long runs cross several line feed search windows and reader buffers.
    for pos in [
        1100, 1279, 1280, 3000, 3001, 3002, 5183, 9999, 19_990, 23_000, 262_100, 262_200, 524_300,
        599_990,
    ] {
        seq[pos] = b'N';
    }
    let expected: Vec<&[u8]> = seq
        .split(|b| !b"ACGTacgt".contains(b))
        .filter(|run| !run.is_empty())
        .collect();
    let fasta = [b">r\n".as_slice(), &seq, b"\n"].concat();
    let wrapped: Vec<u8> = std::iter::once(b">r".as_slice())
        .chain(seq.chunks(60))
        .flat_map(|line| [line, b"\n"].concat())
        .collect();
    let qual = vec![b'I'; seq.len()];
    let fastq = [b"@r\n".as_slice(), &seq, b"\n+\n", &qual, b"\n"].concat();

    const CFG: Config = ParserOptions::default()
        .dna_packed()
        .and_dna_columnar()
        .config();
    fn check<P: HelicaseParser + Iterator<Item = Result<Event, parser::Error>>>(
        mut parser: P,
        expected: &[&[u8]],
    ) {
        let mut chunks = 0;
        while let Some(event) = parser.next() {
            if let Event::DnaChunk(_) = event.unwrap() {
                let run = expected[chunks];
                let packed = parser.get_dna_packed();
                let columnar = parser.get_dna_columnar();
                assert_eq!(packed.len(), run.len(), "chunk {chunks}");
                assert_eq!(columnar.len(), run.len(), "chunk {chunks}");
                for (i, &b) in run.iter().enumerate() {
                    assert_eq!(packed.get_char(i), b.to_ascii_uppercase() as char);
                    assert_eq!(columnar.get_char(i), b.to_ascii_uppercase() as char);
                }
                chunks += 1;
            }
        }
        assert_eq!(chunks, expected.len());
    }

    const SKIP: Config = ParserOptions::default()
        .dna_packed()
        .and_dna_columnar()
        .skip_non_actg()
        .config();
    fn check_skip<P: HelicaseParser + Iterator<Item = Result<Event, parser::Error>>>(
        mut parser: P,
        merged: &[u8],
    ) {
        assert!(matches!(parser.next(), Some(Ok(Event::Record(_)))));
        let packed = parser.get_dna_packed();
        let columnar = parser.get_dna_columnar();
        assert_eq!(packed.len(), merged.len());
        assert_eq!(columnar.len(), merged.len());
        for (i, &b) in merged.iter().enumerate() {
            assert_eq!(
                packed.get_char(i),
                b.to_ascii_uppercase() as char,
                "base {i}"
            );
            assert_eq!(
                columnar.get_char(i),
                b.to_ascii_uppercase() as char,
                "base {i}"
            );
        }
        assert!(parser.next().is_none());
    }

    let merged = expected.concat();
    for data in [&fasta, &wrapped, &fastq] {
        check(FastxParser::<CFG>::from_slice(data).unwrap(), &expected);
        check(
            FastxParser::<CFG>::from_reader(data.as_slice()).unwrap(),
            &expected,
        );
        check_skip(FastxParser::<SKIP>::from_slice(data).unwrap(), &merged);
        check_skip(
            FastxParser::<SKIP>::from_reader(data.as_slice()).unwrap(),
            &merged,
        );
    }
}

#[test]
fn test_config_keep_non_actg_with_masks() {
    // Generate records containing Ns and IUPAC bases
    let bytes = FastxGen::fasta(10, 120)
        .with_alphabet(Alphabet::DnaIupac)
        .with_seed(8888)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .keep_non_actg()
        .compute_mask_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string();
        let mask_non_actg = parser.get_mask_non_actg();

        assert_eq!(
            mask_non_actg.len(),
            s.len(),
            "mask_non_actg length mismatch"
        );

        for (i, &b) in s.iter().enumerate() {
            let is_actg = matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');

            assert_eq!(
                mask_non_actg.get(i),
                !is_actg,
                "mask_non_actg mismatch at index {i} for base '{}'",
                b as char
            );
        }
    }
}

#[test]
fn test_config_keep_non_actg_with_masks_owned() {
    let bytes = FastxGen::fasta(10, 120)
        .with_alphabet(Alphabet::DnaIupac)
        .with_seed(8888)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .keep_non_actg()
        .compute_mask_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string_owned();
        let mask_non_actg = parser.get_mask_non_actg_owned();

        assert_eq!(
            mask_non_actg.len(),
            s.len(),
            "mask_non_actg length mismatch"
        );

        for (i, &b) in s.iter().enumerate() {
            let is_actg = matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');

            assert_eq!(
                mask_non_actg.get(i),
                !is_actg,
                "mask_non_actg mismatch at index {i} for base '{}'",
                b as char
            );
        }
    }
}

#[test]
fn test_config_keep_non_actg_with_mask_n() {
    let bytes = FastxGen::fasta(10, 120)
        .with_alphabet(Alphabet::DnaIupac)
        .with_seed(8888)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .keep_non_actg()
        .compute_mask_non_actg()
        .config()
        | advanced::COMPUTE_MASK_N;

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string();
        let mask_non_actg = parser.get_mask_non_actg();
        let mask_n = parser.get_mask_n();

        assert_eq!(
            mask_non_actg.len(),
            s.len(),
            "mask_non_actg length mismatch"
        );
        assert_eq!(mask_n.len(), s.len(), "mask_n length mismatch");

        for (i, &b) in s.iter().enumerate() {
            let is_actg = matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
            let is_n = matches!(b, b'N' | b'n');

            assert_eq!(
                mask_non_actg.get(i),
                !is_actg,
                "mask_non_actg mismatch at index {i} for base '{}'",
                b as char
            );
            assert_eq!(
                mask_n.get(i),
                is_n,
                "mask_n mismatch at index {i} for base '{}'",
                b as char
            );
        }
    }
}

#[test]
fn test_config_keep_non_actg_with_mask_n_owned() {
    let bytes = FastxGen::fasta(10, 120)
        .with_alphabet(Alphabet::DnaIupac)
        .with_seed(8888)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .dna_string()
        .keep_non_actg()
        .compute_mask_non_actg()
        .config()
        | advanced::COMPUTE_MASK_N;

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let s = parser.get_dna_string_owned();
        let mask_non_actg = parser.get_mask_non_actg_owned();
        let mask_n = parser.get_mask_n_owned();

        assert_eq!(
            mask_non_actg.len(),
            s.len(),
            "mask_non_actg length mismatch"
        );
        assert_eq!(mask_n.len(), s.len(), "mask_n length mismatch");

        for (i, &b) in s.iter().enumerate() {
            let is_actg = matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
            let is_n = matches!(b, b'N' | b'n');

            assert_eq!(
                mask_non_actg.get(i),
                !is_actg,
                "mask_non_actg mismatch at index {i} for base '{}'",
                b as char
            );
            assert_eq!(
                mask_n.get(i),
                is_n,
                "mask_n mismatch at index {i} for base '{}'",
                b as char
            );
        }
    }
}

#[test]
fn test_config_skip_non_actg_merges_chunks() {
    let raw = b">test\nACGTNNNACGT\n";
    const CFG: Config = ParserOptions::default()
        .dna_string()
        .skip_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(raw).unwrap();
    let mut count = 0;
    while let Some(Ok(event)) = parser.next() {
        if let Event::Record(_) = event {
            count += 1;
            assert_eq!(parser.get_dna_string(), b"ACGTACGT");
        }
    }
    assert_eq!(count, 1);
}

#[test]
fn test_config_dna_len() {
    let bytes = FastxGen::illumina(20)
        .with_mean_len(180)
        .with_seed(5555)
        .generate_bytes();

    const CFG: Config = ParserOptions::default()
        .compute_dna_len()
        .dna_string()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        assert_eq!(parser.get_dna_len(), parser.get_dna_string().len());
    }
}

#[test]
fn test_config_dna_len_merge_records() {
    use helicase::config::advanced::*;

    const CFG: Config = COMPUTE_DNA_LEN | SPLIT_NON_ACTG | MERGE_DNA_CHUNKS | MERGE_RECORDS;

    // FASTA: 3 records of 8 bp = 24 bp
    let fasta_bytes = b">r1\nACGTACGT\n>r2\nTGCAACGT\n>r3\nACGTACGT\n";
    let mut parser = FastxParser::<CFG>::from_slice(fasta_bytes).unwrap();
    parser.next();
    assert_eq!(parser.get_dna_len(), 24);

    // FASTQ: 3 records of 8 bp = 24 bp
    let fastq_bytes =
        b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTGCAACGT\n+\nIIIIIIII\n@r3\nACGTACGT\n+\nIIIIIIII\n";
    let mut parser = FastxParser::<CFG>::from_slice(fastq_bytes).unwrap();
    parser.next();
    assert_eq!(parser.get_dna_len(), 24);
}

#[test]
fn test_config_ignore_headers() {
    let bytes = FastxGen::illumina(20).with_seed(101).generate_bytes();
    const CFG: Config = ParserOptions::default()
        .compute_quality()
        .ignore_headers()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    let mut count = 0;
    while let Some(_event) = parser.next() {
        count += 1;
        assert_eq!(parser.get_dna_string().len(), 150);
    }
    assert_eq!(count, 20);
}

#[test]
fn test_config_ignore_quality_fastq() {
    let bytes = FastxGen::illumina(20).with_seed(202).generate_bytes();
    // Default config ignores quality
    const CFG: Config = ParserOptions::default().config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    let mut count = 0;
    while let Some(_event) = parser.next() {
        let expected_header = format!("record_{count}");
        assert_eq!(parser.get_header(), expected_header.as_bytes());
        count += 1;
    }
    assert_eq!(count, 20);
}

#[test]
fn test_config_ignore_dna() {
    let bytes = FastxGen::illumina(25).with_seed(303).generate_bytes();
    const CFG: Config = ParserOptions::default()
        .compute_quality()
        .ignore_dna()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    let mut count = 0;
    while let Some(Ok(event)) = parser.next() {
        if let Event::Record(_) = event {
            let expected_header = format!("record_{count}");
            assert_eq!(parser.get_header(), expected_header.as_bytes());
            let qual = parser.get_quality().unwrap();
            assert_eq!(qual.len(), 150);
            count += 1;
        }
    }
    assert_eq!(count, 25);
}

#[test]
fn test_config_split_non_actg_return_both_chunk_and_record() {
    let raw = b">read1\nACGTNNNACGT\n";
    const CFG: Config = ParserOptions::default()
        .dna_string()
        .split_non_actg()
        .return_record(true)
        .return_dna_chunk(true)
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(raw).unwrap();
    let mut chunk_count = 0;
    let mut record_count = 0;
    while let Some(Ok(event)) = parser.next() {
        match event {
            Event::DnaChunk(_) => {
                chunk_count += 1;
                assert_eq!(parser.get_dna_string(), b"ACGT");
            }
            Event::Record(_) => {
                record_count += 1;
                assert_eq!(parser.get_header(), b"read1");
            }
        }
    }
    assert_eq!(chunk_count, 2);
    assert_eq!(record_count, 1);
}

#[cfg(feature = "packed-seq")]
#[test]
fn test_config_packed_seq_trait_integration() {
    use packed_seq::{Seq, SeqVec};

    let bytes = FastxGen::illumina(5).with_seed(404).generate_bytes();
    const CFG: Config = ParserOptions::default()
        .dna_packed()
        .keep_non_actg()
        .config();

    let mut parser = FastxParser::<CFG>::from_slice(&bytes).unwrap();
    while let Some(_event) = parser.next() {
        let seq_ref = parser.get_packed_seq();
        assert_eq!(seq_ref.len(), 150);

        let seq_vec = parser.get_packed_seq_vec();
        assert_eq!(seq_vec.len(), 150);
    }
}
