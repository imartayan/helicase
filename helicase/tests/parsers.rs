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
        let mut f = parser::FastaParser::<CONFIG_HEADER, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
        }
        assert_eq!(res, vec!["head", "h>h", "A B C ",]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = parser::FastaParser::<CONFIG_STRING, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC",]
        );

        let mut f = parser::FastaParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTA).unwrap();
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

        let mut f = parser::FastaParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTA).unwrap();
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
        let mut f = parser::FastaParser::<CONFIG_COLUMNAR, _>::from_slice(FASTA).unwrap();
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

        let mut f = parser::FastaParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTA).unwrap();
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
        let mut f = parser::FastaParser::<CONFIG_PACKED, _>::from_slice(FASTA).unwrap();
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

        let mut f = parser::FastaParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTA).unwrap();
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
        let mut f = parser::FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
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
        let mut f = parser::FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
        while let Some(_) = f.next() {
            let string_len = f.get_dna_string().len();
            let packed_len = f.get_dna_packed().len();
            assert_eq!(packed_len, string_len, "packed length mismatch");
        }
    }

    #[test]
    fn test_crlf_fasta() {
        let crlf_data = b">seq1\r\nACGTACGT\r\n>seq2\r\nAAAA\r\nCCCC\r\nGGGG\r\n>seq3\r\nTTTT\r\n";
        const CONFIG_ALL: Config = ParserOptions::default()
            .compute_headers()
            .dna_string()
            .and_dna_packed()
            .and_dna_columnar()
            .keep_non_actg()
            .config();

        let mut parser = parser::FastaParser::<CONFIG_ALL, _>::from_slice(crlf_data).unwrap();
        let mut headers = Vec::new();
        let mut seqs = Vec::new();
        while let Some(_event) = parser.next() {
            let hdr = String::from_utf8(parser.get_header().to_vec()).unwrap();
            let dna = String::from_utf8(parser.get_dna_string().to_vec()).unwrap();
            headers.push(hdr);
            assert_eq!(parser.get_dna_packed().to_string(), dna);
            assert_eq!(parser.get_dna_columnar().to_string(), dna);
            seqs.push(dna);
        }
        assert_eq!(headers, vec!["seq1", "seq2", "seq3"]);
        assert_eq!(seqs, vec!["ACGTACGT", "AAAACCCCGGGG", "TTTT"]);
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
        let mut f = parser::FastqParser::<CONFIG_HEADER, _>::from_slice(FASTQ).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_QUALITY, _>::from_slice(FASTQ).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_STRING, _>::from_slice(FASTQ).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC"]
        );

        let mut f = parser::FastqParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTQ).unwrap();
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

        let mut f = parser::FastqParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTQ).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_COLUMNAR, _>::from_slice(FASTQ).unwrap();
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

        let mut f = parser::FastqParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTQ).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_PACKED, _>::from_slice(FASTQ).unwrap();
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

        let mut f = parser::FastqParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTQ).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fastq).unwrap();
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
        let mut f = parser::FastqParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fastq).unwrap();
        while let Some(_) = f.next() {
            let string_len = f.get_dna_string().len();
            let packed_len = f.get_dna_packed().len();
            assert_eq!(packed_len, string_len, "packed length mismatch");
        }
    }

    #[test]
    fn test_crlf_fastq() {
        let crlf_data =
            b"@read1\r\nACGTACGT\r\n+\r\nIIIIIIII\r\n@read2\r\nTTTTCCCC\r\n+\r\nJJJJJJJJ\r\n";
        const CONFIG_ALL: Config = ParserOptions::default()
            .compute_headers()
            .dna_string()
            .and_dna_packed()
            .and_dna_columnar()
            .keep_non_actg()
            .compute_quality()
            .config();

        let mut parser = parser::FastqParser::<CONFIG_ALL, _>::from_slice(crlf_data).unwrap();
        let mut headers = Vec::new();
        let mut seqs = Vec::new();
        let mut quals = Vec::new();
        while let Some(_event) = parser.next() {
            headers.push(String::from_utf8(parser.get_header().to_vec()).unwrap());
            let dna = String::from_utf8(parser.get_dna_string().to_vec()).unwrap();
            quals.push(String::from_utf8(parser.get_quality().unwrap().to_vec()).unwrap());
            assert_eq!(parser.get_dna_packed().to_string(), dna);
            assert_eq!(parser.get_dna_columnar().to_string(), dna);
            seqs.push(dna);
        }
        assert_eq!(headers, vec!["read1", "read2"]);
        assert_eq!(seqs, vec!["ACGTACGT", "TTTTCCCC"]);
        assert_eq!(quals, vec!["IIIIIIII", "JJJJJJJJ"]);
    }

    #[test]
    fn test_crlf_fastq_split_non_actg() {
        let crlf_data = b"@read1\r\nACGTNNACGT\r\n+\r\nIIIIIIIIII\r\n";
        const CONFIG_SPLIT: Config = ParserOptions::default()
            .compute_headers()
            .dna_string()
            .split_non_actg()
            .config();

        let mut parser = parser::FastqParser::<CONFIG_SPLIT, _>::from_slice(crlf_data).unwrap();
        let mut chunks = Vec::new();
        while let Some(_event) = parser.next() {
            chunks.push(String::from_utf8(parser.get_dna_string().to_vec()).unwrap());
        }
        assert_eq!(chunks, vec!["ACGT", "ACGT"]);
    }

    #[test]
    fn test_crlf_block_boundary() {
        // Construct header such that sequence ends with \r at byte 63 and \n at byte 64
        let header_prefix = "A".repeat(57);
        let fastq = format!("@{header_prefix}\r\nACG\r\n+\r\nIII\r\n");
        const CONFIG_ALL: Config = ParserOptions::default()
            .compute_headers()
            .dna_string()
            .and_dna_packed()
            .and_dna_columnar()
            .keep_non_actg()
            .compute_quality()
            .config();

        let mut parser =
            parser::FastqParser::<CONFIG_ALL, _>::from_slice(fastq.as_bytes()).unwrap();
        let mut count = 0;
        while let Some(_event) = parser.next() {
            count += 1;
            assert_eq!(parser.get_header(), header_prefix.as_bytes());
            assert_eq!(parser.get_dna_string(), b"ACG");
            assert_eq!(parser.get_quality(), Some(b"III".as_slice()));
            assert_eq!(parser.get_dna_packed().to_string(), "ACG");
            assert_eq!(parser.get_dna_columnar().to_string(), "ACG");
        }
        assert_eq!(count, 1);
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
        let mut parser = parser::FastaParser::<CONFIG, _>::from_slice(fasta_data).unwrap();
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
        let mut parser = parser::FastqParser::<CONFIG, _>::from_slice(fastq_data).unwrap();
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

#[test]
fn test_find_record_start() {
    use helicase::parser::Format;

    let mut fastq = Vec::new();
    let mut fasta = Vec::new();
    for eol in ["\n", "\r\n"] {
        for seq in ["", "A", "AN"] {
            for qual_start in ['@', '+', 'I'] {
                for plus in ["+", "+@"] {
                    let qual: String = seq
                        .chars()
                        .enumerate()
                        .map(|(i, _)| if i == 0 { qual_start } else { '@' })
                        .collect();
                    fastq.push(format!("@r{eol}{seq}{eol}{plus}{eol}{qual}{eol}"));
                }
            }
        }
        for lines in ["", "A", "A\nN"] {
            fasta.push(format!(
                ">r{eol}{}",
                lines
                    .lines()
                    .map(|l| format!("{l}{eol}"))
                    .collect::<String>()
            ));
        }
    }
    fastq.sort();
    fastq.dedup();

    // Expected result on the window `data[w0..w1]` at local position `pos`.
    let expected = |format, data: &[u8], starts: &[usize], w0: usize, w1: usize, pos: usize| {
        let start = *starts.iter().find(|&&s| s > w0 && s >= w0 + pos)?;
        let fits = match format {
            Format::Fasta => start < w1,
            Format::Fastq => start + data[start..].iter().position(|&b| b == b'\n')? + 1 < w1,
        };
        fits.then_some(start - w0)
    };
    for (records, format) in [(&fastq, Format::Fastq), (&fasta, Format::Fasta)] {
        for num_records in [2, 3] {
            let mut combination = vec![0; num_records];
            loop {
                let mut data = Vec::new();
                let mut starts = Vec::new();
                for &r in &combination {
                    starts.push(data.len());
                    data.extend_from_slice(records[r].as_bytes());
                }
                // All windows on 2 records, whole input on 3 records.
                let windows: Vec<(usize, usize)> = if num_records == 2 {
                    (0..=data.len())
                        .flat_map(|w0| (w0..=data.len()).map(move |w1| (w0, w1)))
                        .collect()
                } else {
                    vec![(0, data.len())]
                };
                for (w0, w1) in windows {
                    for pos in 0..=w1 - w0 + 1 {
                        assert_eq!(
                            format.find_record_start(&data[w0..w1], pos),
                            expected(format, &data, &starts, w0, w1, pos),
                            "{format:?} {:?} window {w0}..{w1} at {pos}",
                            String::from_utf8_lossy(&data)
                        );
                    }
                }
                let Some(i) = combination.iter().rposition(|&r| r + 1 < records.len()) else {
                    break;
                };
                combination[i] += 1;
                combination[i + 1..].fill(0);
            }
        }
    }
}
