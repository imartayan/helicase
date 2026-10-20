#[cfg(feature = "zstd")]
mod tests {
    use helicase::config::*;
    use helicase::input::*;
    use helicase::*;

    const CONFIG: Config = ParserOptions::default().config();

    fn generate_test_fastq(n: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for i in 0..n {
            let rec = format!(
                "@read_{i}\nACGTACGTACGTNNACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
            );
            data.extend_from_slice(rec.as_bytes());
        }
        data
    }

    fn generate_test_fasta(n: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for i in 0..n {
            let rec = format!(">seq_{i}\nACGTACGTACGTNNACGTACGTACGTACGT\n");
            data.extend_from_slice(rec.as_bytes());
        }
        data
    }

    #[test]
    fn test_parallel_fastx_parser_zstd_multiframe_fastq() {
        let chunk1 = generate_test_fastq(1250);
        let chunk2 = generate_test_fastq(1250);
        let chunk3 = generate_test_fastq(1250);
        let chunk4 = generate_test_fastq(1250);

        let mut zst_data = Vec::new();
        for chunk in [&chunk1, &chunk2, &chunk3, &chunk4] {
            let compressed = zstd::bulk::compress(chunk, 1).unwrap();
            zst_data.extend_from_slice(&compressed);
        }

        let temp_file = "/tmp/test_parallel_zstd_multiframe.fastq.zst";
        std::fs::write(temp_file, &zst_data).unwrap();

        for threads in [1, 2, 4, 8] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut count = 0;
                        let mut len = 0;
                        while parser.next().is_some() {
                            count += 1;
                            len += parser.get_dna_len();
                        }
                        (count, len)
                    },
                    |a, b| (a.0 + b.0, a.1 + b.1),
                )
                .unwrap();

            assert_eq!(total_records, 5000, "Failed at threads={threads}");
            assert_eq!(total_len, 5000 * 30, "Failed at threads={threads}");
        }

        // Test single-threaded from_file and from_reader
        {
            let mut parser = FastxParser::<CONFIG>::from_file(temp_file).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(count, 5000, "Single-threaded from_file mismatch");
            assert_eq!(len, 5000 * 30, "Single-threaded from_file len mismatch");
        }
        {
            let cursor = std::io::Cursor::new(zst_data.clone());
            let mut parser = FastxParser::<CONFIG>::from_reader(cursor).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(count, 5000, "Single-threaded from_reader mismatch");
            assert_eq!(len, 5000 * 30, "Single-threaded from_reader len mismatch");
        }

        let _ = std::fs::remove_file(temp_file);
    }

    #[test]
    fn test_parallel_fastx_parser_zstd_multiframe_fasta() {
        let chunk1 = generate_test_fasta(1000);
        let chunk2 = generate_test_fasta(1000);
        let chunk3 = generate_test_fasta(1000);

        let mut zst_data = Vec::new();
        for chunk in [&chunk1, &chunk2, &chunk3] {
            let compressed = zstd::bulk::compress(chunk, 1).unwrap();
            zst_data.extend_from_slice(&compressed);
        }

        let temp_file = "/tmp/test_parallel_zstd_multiframe.fasta.zst";
        std::fs::write(temp_file, &zst_data).unwrap();

        for threads in [1, 2, 4] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut count = 0;
                        let mut len = 0;
                        while parser.next().is_some() {
                            count += 1;
                            len += parser.get_dna_len();
                        }
                        (count, len)
                    },
                    |a, b| (a.0 + b.0, a.1 + b.1),
                )
                .unwrap();

            assert_eq!(total_records, 3000, "Failed at threads={threads}");
            assert_eq!(total_len, 3000 * 30, "Failed at threads={threads}");
        }

        let _ = std::fs::remove_file(temp_file);
    }

    #[test]
    fn test_parallel_fastx_parser_zstd_singleframe() {
        let data = generate_test_fastq(1000);
        let compressed = zstd::bulk::compress(&data, 1).unwrap();

        let temp_file = "/tmp/test_parallel_zstd_singleframe.fastq.zst";
        std::fs::write(temp_file, &compressed).unwrap();

        for threads in [1, 2, 4] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut count = 0;
                        let mut len = 0;
                        while parser.next().is_some() {
                            count += 1;
                            len += parser.get_dna_len();
                        }
                        (count, len)
                    },
                    |a, b| (a.0 + b.0, a.1 + b.1),
                )
                .unwrap();

            assert_eq!(total_records, 1000, "Failed at threads={threads}");
            assert_eq!(total_len, 1000 * 30, "Failed at threads={threads}");
        }

        let _ = std::fs::remove_file(temp_file);
    }

    #[test]
    fn test_parallel_fastx_parser_zstd_records_crossing_frame_boundary() {
        let raw_fastq = generate_test_fastq(3000);

        // Intentionally split raw bytes at arbitrary offsets (e.g. 500 bytes)
        // so that records are cut across frame boundaries.
        let mut zst_data = Vec::new();
        let frame_size = 500;
        let mut offset = 0;
        while offset < raw_fastq.len() {
            let end = (offset + frame_size).min(raw_fastq.len());
            let compressed = zstd::bulk::compress(&raw_fastq[offset..end], 1).unwrap();
            zst_data.extend_from_slice(&compressed);
            offset = end;
        }

        let temp_file = "/tmp/test_parallel_zstd_crossing_fastq.zst";
        std::fs::write(temp_file, &zst_data).unwrap();

        for threads in [1, 2, 4, 8] {
            let (total_records, total_len, records) =
                ParallelFastxParser::<CONFIG>::from_file(temp_file)
                    .unwrap()
                    .threads(threads)
                    .reduce(
                        |mut parser| {
                            let mut count = 0;
                            let mut len = 0;
                            let mut recs = Vec::new();
                            while parser.next().is_some() {
                                count += 1;
                                len += parser.get_dna_len();
                                recs.push((
                                    parser.get_header().to_vec(),
                                    parser.get_dna_string().to_vec(),
                                ));
                            }
                            (count, len, recs)
                        },
                        |mut a, b| {
                            a.0 += b.0;
                            a.1 += b.1;
                            a.2.extend(b.2);
                            a
                        },
                    )
                    .unwrap();

            assert_eq!(total_records, 3000, "Failed at threads={threads}");
            assert_eq!(total_len, 3000 * 30, "Failed at threads={threads}");

            // Verify exact sequential match with original records
            for (i, (hdr, seq)) in records.iter().enumerate() {
                let expected_hdr = format!("read_{i}");
                let expected_seq = b"ACGTACGTACGTNNACGTACGTACGTACGT";
                assert_eq!(
                    hdr,
                    expected_hdr.as_bytes(),
                    "Header mismatch at record {i}"
                );
                assert_eq!(seq, expected_seq, "Sequence mismatch at record {i}");
            }
        }

        let _ = std::fs::remove_file(temp_file);

        // Same test for FASTA
        let raw_fasta = generate_test_fasta(3000);
        let mut zst_fasta_data = Vec::new();
        let mut offset = 0;
        while offset < raw_fasta.len() {
            let end = (offset + frame_size).min(raw_fasta.len());
            let compressed = zstd::bulk::compress(&raw_fasta[offset..end], 1).unwrap();
            zst_fasta_data.extend_from_slice(&compressed);
            offset = end;
        }

        let temp_fasta_file = "/tmp/test_parallel_zstd_crossing_fasta.zst";
        std::fs::write(temp_fasta_file, &zst_fasta_data).unwrap();

        for threads in [1, 2, 4, 8] {
            let (total_records, total_len) =
                ParallelFastxParser::<CONFIG>::from_file(temp_fasta_file)
                    .unwrap()
                    .threads(threads)
                    .reduce(
                        |mut parser| {
                            let mut count = 0;
                            let mut len = 0;
                            while parser.next().is_some() {
                                count += 1;
                                len += parser.get_dna_len();
                            }
                            (count, len)
                        },
                        |a, b| (a.0 + b.0, a.1 + b.1),
                    )
                    .unwrap();

            assert_eq!(total_records, 3000, "FASTA failed at threads={threads}");
            assert_eq!(total_len, 3000 * 30, "FASTA failed at threads={threads}");
        }

        let _ = std::fs::remove_file(temp_fasta_file);
    }

    #[test]
    fn test_parallel_fastx_parser_zstd_skippable_frames() {
        let chunk1 = generate_test_fastq(1000);
        let chunk2 = generate_test_fastq(1000);

        let c1 = zstd::bulk::compress(&chunk1, 1).unwrap();
        let c2 = zstd::bulk::compress(&chunk2, 1).unwrap();

        // Skippable frame 1: magic 0x184D2A50, len 4, data [1, 2, 3, 4]
        let skip1: [u8; 12] = [
            0x50, 0x2a, 0x4d, 0x18, // magic
            0x04, 0x00, 0x00, 0x00, // payload len
            0x01, 0x02, 0x03, 0x04, // payload
        ];

        // Skippable frame 2: magic 0x184D2A55, len 8, data
        let skip2: [u8; 16] = [
            0x55, 0x2a, 0x4d, 0x18, // magic
            0x08, 0x00, 0x00, 0x00, // payload len
            0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x11, 0x22,
        ];

        // Interleave skippable frames: skip1 + c1 + skip2 + c2
        let mut zst_data = Vec::new();
        zst_data.extend_from_slice(&skip1);
        zst_data.extend_from_slice(&c1);
        zst_data.extend_from_slice(&skip2);
        zst_data.extend_from_slice(&c2);

        let temp_file = "/tmp/test_parallel_zstd_skippable.fastq.zst";
        std::fs::write(temp_file, &zst_data).unwrap();

        for threads in [1, 2, 4] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut count = 0;
                        let mut len = 0;
                        while parser.next().is_some() {
                            count += 1;
                            len += parser.get_dna_len();
                        }
                        (count, len)
                    },
                    |a, b| (a.0 + b.0, a.1 + b.1),
                )
                .unwrap();

            assert_eq!(total_records, 2000, "Failed at threads={threads}");
            assert_eq!(total_len, 2000 * 30, "Failed at threads={threads}");
        }

        // Also test single-threaded from_file
        {
            let mut parser = FastxParser::<CONFIG>::from_file(temp_file).unwrap();
            let mut count = 0;
            while parser.next().is_some() {
                count += 1;
            }
            assert_eq!(count, 2000, "Single-threaded with skippable frames failed");
        }

        let _ = std::fs::remove_file(temp_file);
    }
}
