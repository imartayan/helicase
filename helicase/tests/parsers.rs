use helicase::config::*;
use helicase::input::*;
use helicase::*;

mod fasta {
    use super::*;

    const CONFIG_HEADER: Config = ParserOptions::default().ignore_dna().config();
    const CONFIG_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .config();
    const CONFIG_STRING_ACTG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .config();
    const CONFIG_STRING_ACTG_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .skip_non_actg()
        .config();
    const CONFIG_COLUMNAR: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .config();
    const CONFIG_COLUMNAR_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .config();
    const CONFIG_PACKED_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED_KEEP: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .and_dna_packed()
        .keep_non_actg()
        .config();

    static FASTA: &[u8] =
        b">head\nTTTCTtaAAAA\nAGAAAA\nACAAN\n\n>h>h\nCTCTTANNAAA\nCAAAnAGCTTT\n>A B C \nCCAC"
            .as_slice();

    #[test]
    fn test_header() {
        let mut f = FastaParser::<CONFIG_HEADER, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
        }
        assert_eq!(res, vec!["head", "h>h", "A B C ",]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = FastaParser::<CONFIG_STRING, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC",]
        );

        let mut f = FastaParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec![
                "TTTCTtaAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_dna_columnar() {
        let mut f = FastaParser::<CONFIG_COLUMNAR, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_dna_packed() {
        let mut f = FastaParser::<CONFIG_PACKED, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_packed_matches_string_long_sequence() {
        let seq = b"ACGT".repeat(40);
        let fasta = [b">seq\n".as_ref(), seq.as_ref()].concat();
        let mut f = FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
        while let Some(_) = f.next() {
            let s = f.get_dna_string();
            let p = f.get_dna_packed();
            assert_eq!(p.len(), s.len(), "length mismatch");
            for (i, base) in s.iter().enumerate() {
                assert_eq!(
                    p.get_char(i) as u8,
                    base.to_ascii_uppercase(),
                    "mismatch at base {i}"
                );
            }
        }
    }

    #[test]
    fn test_packed_length_keep_non_actg_last_chunk_all_dna() {
        let seq = b"A".repeat(194);
        let fasta = [b">s\n".as_ref(), seq.as_ref()].concat();
        assert_eq!(fasta.len() % 64, 197 % 64);
        let mut f = FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
        while let Some(_) = f.next() {
            let string_len = f.get_dna_string().len();
            let packed_len = f.get_dna_packed().len();
            assert_eq!(packed_len, string_len, "packed length mismatch");
        }
    }
}

mod fastq {
    use super::*;

    const CONFIG_HEADER: Config = ParserOptions::default().ignore_dna().config();
    const CONFIG_QUALITY: Config = ParserOptions::default()
        .ignore_headers()
        .ignore_dna()
        .compute_quality()
        .config();
    const CONFIG_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .config();
    const CONFIG_STRING_ACTG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .config();
    const CONFIG_STRING_ACTG_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .skip_non_actg()
        .config();
    const CONFIG_COLUMNAR: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .config();
    const CONFIG_COLUMNAR_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .config();
    const CONFIG_PACKED_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED_KEEP: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .and_dna_packed()
        .keep_non_actg()
        .config();

    static FASTQ: &[u8] =
        b"@head\nTTTCTtaAAAAAGAAAAACAAN\n+\n123\n@hhh\nCTCTTANNAAACAAAnAGCTTT\n+\nQQ@@++AA\n@A B C \nCCAC\n+\nQUAL"
            .as_slice();

    #[test]
    fn test_header() {
        let mut f = FastqParser::<CONFIG_HEADER, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        let mut c = 0;
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
            c += 1;
            if c > 3 {
                break;
            }
        }
        assert_eq!(res, vec!["head", "hhh", "A B C "]);
    }

    #[test]
    fn test_quality() {
        let mut f = FastqParser::<CONFIG_QUALITY, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        let mut c = 0;
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_quality_owned().unwrap()).unwrap());
            c += 1;
            if c > 3 {
                break;
            }
        }
        assert_eq!(res, vec!["123", "QQ@@++AA", "QUAL"]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = FastqParser::<CONFIG_STRING, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC"]
        );

        let mut f = FastqParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec![
                "TTTCTtaAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC"]
        );
    }

    #[test]
    fn test_dna_columnar() {
        let mut f = FastqParser::<CONFIG_COLUMNAR, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC"]
        );
    }

    #[test]
    fn test_dna_packed() {
        let mut f = FastqParser::<CONFIG_PACKED, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC"]
        );
    }

    #[test]
    fn test_packed_matches_string_long_sequence() {
        let seq = b"ACGT".repeat(40);
        let qual = b"I".repeat(160);
        let fastq = [b"@seq\n".as_ref(), &seq, b"\n+\n".as_ref(), &qual].concat();
        let mut f = FastqParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fastq).unwrap();
        while let Some(_) = f.next() {
            let s = f.get_dna_string();
            let p = f.get_dna_packed();
            assert_eq!(p.len(), s.len(), "length mismatch");
            for (i, base) in s.iter().enumerate() {
                assert_eq!(
                    p.get_char(i) as u8,
                    base.to_ascii_uppercase(),
                    "mismatch at base {i}"
                );
            }
        }
    }

    #[test]
    fn test_packed_length_keep_non_actg_last_chunk_all_dna() {
        let seq = b"A".repeat(194);
        let qual = b"I".repeat(194);
        let fastq = [b"@s\n".as_ref(), &seq, b"\n+\n".as_ref(), &qual].concat();
        let mut f = FastqParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fastq).unwrap();
        while let Some(_) = f.next() {
            let string_len = f.get_dna_string().len();
            let packed_len = f.get_dna_packed().len();
            assert_eq!(packed_len, string_len, "packed length mismatch");
        }
    }
}

#[cfg(feature = "paraseq")]
mod paraseq_tests {
    use super::*;
    use paraseq::Record;

    const CONFIG: Config = ParserOptions::default().dna_string().config();

    #[test]
    fn test_fasta_record_index() {
        let fasta_data = b">r0\nACGT\n>r1\nTGCA\nTGCA\n>r2\nNNNN\n";
        let mut parser = FastaParser::<CONFIG, _>::from_slice(fasta_data).unwrap();
        let mut indices = Vec::new();
        let mut trait_indices = Vec::new();
        while let Some(_) = parser.next() {
            indices.push(parser.record_index());
            trait_indices.push(Record::index(&parser));
        }
        assert_eq!(indices, vec![0, 1, 2]);
        assert_eq!(trait_indices, vec![0, 1, 2]);
    }

    #[test]
    fn test_fastq_record_index() {
        let fastq_data = b"@r0\nACGT\n+\nIIII\n@r1\nTGCA\n+\nIIII\n@r2\nGGCC\n+\nIIII\n";
        let mut parser = FastqParser::<CONFIG, _>::from_slice(fastq_data).unwrap();
        let mut indices = Vec::new();
        let mut trait_indices = Vec::new();
        while let Some(_) = parser.next() {
            indices.push(parser.record_index());
            trait_indices.push(Record::index(&parser));
        }
        assert_eq!(indices, vec![0, 1, 2]);
        assert_eq!(trait_indices, vec![0, 1, 2]);
    }

    #[test]
    fn test_fastx_record_index() {
        let fasta_data = b">f0\nACGT\n>f1\nTGCA\n";
        let mut parser = FastxParser::<CONFIG>::from_slice(fasta_data).unwrap();
        let mut indices = Vec::new();
        let mut trait_indices = Vec::new();
        while let Some(_) = parser.next() {
            indices.push(parser.record_index());
            trait_indices.push(Record::index(&parser));
        }
        assert_eq!(indices, vec![0, 1]);
        assert_eq!(trait_indices, vec![0, 1]);

        let fastq_data = b"@q0\nACGT\n+\nIIII\n@q1\nTGCA\n+\nIIII\n";
        let mut parser = FastxParser::<CONFIG>::from_slice(fastq_data).unwrap();
        let mut indices = Vec::new();
        let mut trait_indices = Vec::new();
        while let Some(_) = parser.next() {
            indices.push(parser.record_index());
            trait_indices.push(Record::index(&parser));
        }
        assert_eq!(indices, vec![0, 1]);
        assert_eq!(trait_indices, vec![0, 1]);
    }
}
