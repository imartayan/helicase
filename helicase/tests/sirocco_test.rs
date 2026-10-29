#[cfg(feature = "gz")]
mod tests {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use helicase::config::*;
    use helicase::input::*;
    use helicase::*;
    use std::io::Write;

    const CONFIG: Config = ParserOptions::default().config();

    // Encodes each chunk as a gzip block with the BGZF extra field.
    fn encode_bgzf(chunks: &[&[u8]]) -> Vec<u8> {
        let mut bgzf_data = Vec::new();
        for chunk in chunks {
            let mut block_data = Vec::new();
            let mut encoder = flate2::GzBuilder::new()
                .extra([b'B', b'C', 2, 0, 0, 0])
                .write(&mut block_data, flate2::Compression::fast());
            encoder.write_all(chunk).unwrap();
            encoder.finish().unwrap();
            // Update BSIZE in header
            let bsize = (block_data.len() - 1) as u16;
            block_data[16..18].copy_from_slice(&bsize.to_le_bytes());
            bgzf_data.extend_from_slice(&block_data);
        }
        bgzf_data
    }

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

    #[test]
    fn test_parallel_fastx_parser_sirocco_standard_gzip() {
        let original = generate_test_fastq(5000);
        let mut gz_data = Vec::new();
        {
            let mut encoder = GzEncoder::new(&mut gz_data, Compression::default());
            encoder.write_all(&original).unwrap();
            encoder.finish().unwrap();
        }

        let temp_file = "/tmp/test_parallel_sirocco_standard.fastq.gz";
        std::fs::write(temp_file, &gz_data).unwrap();

        for threads in [1, 2, 4] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .for_each_with(
                    || (0, 0),
                    |(count, len), mut parser| {
                        while let Some(event) = parser.next() {
                            event?;
                            *count += 1;
                            *len += parser.get_dna_len();
                        }
                        Ok(())
                    },
                )
                .unwrap()
                .into_iter()
                .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));

            assert_eq!(total_records, 5000, "Failed at threads={threads}");
            assert_eq!(total_len, 5000 * 30, "Failed at threads={threads}");
        }

        // Test single-threaded ReaderInput and FileInput with standard gzip
        {
            let mut parser = FastxParser::<CONFIG>::from_file(temp_file).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(
                count, 5000,
                "Single-threaded from_file standard gzip mismatch"
            );
            assert_eq!(
                len,
                5000 * 30,
                "Single-threaded from_file standard gzip len mismatch"
            );
        }
        {
            let cursor = std::io::Cursor::new(gz_data.clone());
            let mut parser = FastxParser::<CONFIG>::from_reader(cursor).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(
                count, 5000,
                "Single-threaded from_reader standard gzip mismatch"
            );
            assert_eq!(
                len,
                5000 * 30,
                "Single-threaded from_reader standard gzip len mismatch"
            );
        }
    }

    #[test]
    fn test_parallel_fastx_parser_sirocco_bgzf() {
        // Construct a synthetic multi-block BGZF file
        let chunk1 = generate_test_fastq(1000);
        let chunk2 = generate_test_fastq(1000);

        let bgzf_data = encode_bgzf(&[&chunk1, &chunk2]);

        let temp_file = "/tmp/test_parallel_sirocco_bgzf.fastq.gz";
        std::fs::write(temp_file, &bgzf_data).unwrap();

        for threads in [1, 2, 4] {
            let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file(temp_file)
                .unwrap()
                .threads(threads)
                .for_each_with(
                    || (0, 0),
                    |(count, len), mut parser| {
                        while let Some(event) = parser.next() {
                            event?;
                            *count += 1;
                            *len += parser.get_dna_len();
                        }
                        Ok(())
                    },
                )
                .unwrap()
                .into_iter()
                .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));

            assert_eq!(total_records, 2000, "Failed at threads={threads}");
            assert_eq!(total_len, 2000 * 30, "Failed at threads={threads}");
        }

        // Test single-threaded ReaderInput and FileInput with BGZF multi-block data
        {
            let mut parser = FastxParser::<CONFIG>::from_file(temp_file).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(count, 2000, "Single-threaded from_file BGZF mismatch");
            assert_eq!(
                len,
                2000 * 30,
                "Single-threaded from_file BGZF len mismatch"
            );
        }
        {
            let cursor = std::io::Cursor::new(bgzf_data.clone());
            let mut parser = FastxParser::<CONFIG>::from_reader(cursor).unwrap();
            let mut count = 0;
            let mut len = 0;
            while parser.next().is_some() {
                count += 1;
                len += parser.get_dna_len();
            }
            assert_eq!(count, 2000, "Single-threaded from_reader BGZF mismatch");
            assert_eq!(
                len,
                2000 * 30,
                "Single-threaded from_reader BGZF len mismatch"
            );
        }
    }

    #[test]
    fn test_parallel_bgzf_custom_chunks() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let block = generate_test_fastq(1000);
        let bgzf_data = encode_bgzf(&[&block[..]; 8]);
        let path = std::env::temp_dir().join("test_helicase_bgzf_chunks.fastq.gz");
        std::fs::write(&path, &bgzf_data).unwrap();
        for (chunks, expected) in [(None, 2), (Some(3), 3), (Some(100), 8)] {
            let parser = ParallelFastxParser::<CONFIG>::from_file(&path)
                .unwrap()
                .threads(2);
            let parser = match chunks {
                Some(n) => parser.chunks(n),
                None => parser,
            };
            let calls = AtomicUsize::new(0);
            parser
                .for_each(|_| {
                    calls.fetch_add(1, Ordering::Relaxed);
                })
                .unwrap();
            assert_eq!(calls.into_inner(), expected, "chunks({chunks:?})");
        }
        let _ = std::fs::remove_file(path);
    }
}
