use fastx_gen::FastxGen;
use helicase::config::*;
use helicase::input::*;
use helicase::*;

const CONFIG_STR: Config = ParserOptions::default().dna_string().config();
const CONFIG_PACKED: Config = ParserOptions::default().dna_packed().config();
const CONFIG_COL: Config = ParserOptions::default().dna_columnar().config();

#[test]
fn test_parallel_fasta_unwrapped() {
    let bytes = FastxGen::fasta(100, 500).with_seed(1234).generate_bytes();

    let mut single_len = 0;
    let mut single = FastxParser::<CONFIG_STR>::from_slice(&bytes).unwrap();
    while single.next().is_some() {
        single_len += single.get_dna_string().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_len = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_string().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(single_len, par_len, "Length mismatch on {threads} threads");
    }
}

#[test]
fn test_parallel_fasta_wrapped() {
    let bytes = FastxGen::wrapped_fasta(50, 1000, 60)
        .with_seed(5678)
        .generate_bytes();

    let mut single_len = 0;
    let mut single = FastxParser::<CONFIG_STR>::from_slice(&bytes).unwrap();
    while single.next().is_some() {
        single_len += single.get_dna_string().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_len = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_string().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(single_len, par_len, "Length mismatch on {threads} threads");
    }
}

#[test]
fn test_parallel_fastq_with_at_qualities() {
    // FASTQ where quality strings contain `@` and `+` chars
    let raw = b"@read1\nACGTACGT\n+\n@@@@++++\n@read2\nTTTTAAAA\n+\n+@+@+@+@\n@read3\nCCCCCCCC\n+\n@@@@@@@@\n";

    let mut single_len = 0;
    let mut single_records = 0;
    let mut single = FastxParser::<CONFIG_STR>::from_slice(raw).unwrap();
    while single.next().is_some() {
        single_records += 1;
        single_len += single.get_dna_string().len();
    }
    assert_eq!(single_records, 3);
    assert_eq!(single_len, 24);

    for threads in [1, 2, 4] {
        let (par_records, par_len) = ParallelFastxParser::<CONFIG_STR>::from_slice(raw)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut recs = 0;
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        recs += 1;
                        len += parser.get_dna_string().len();
                    }
                    (recs, len)
                },
                |a, b| (a.0 + b.0, a.1 + b.1),
            )
            .unwrap();

        assert_eq!(single_records, par_records);
        assert_eq!(single_len, par_len);
    }
}

#[test]
fn test_parallel_packed_and_columnar() {
    let bytes = FastxGen::illumina(200).with_seed(9999).generate_bytes();

    let mut single_packed_len = 0;
    let mut single_col_len = 0;

    let mut p_pack = FastxParser::<CONFIG_PACKED>::from_slice(&bytes).unwrap();
    while p_pack.next().is_some() {
        single_packed_len += p_pack.get_dna_packed().len();
    }

    let mut p_col = FastxParser::<CONFIG_COL>::from_slice(&bytes).unwrap();
    while p_col.next().is_some() {
        single_col_len += p_col.get_dna_columnar().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_packed = ParallelFastxParser::<CONFIG_PACKED>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_packed().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        let par_col = ParallelFastxParser::<CONFIG_COL>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_columnar().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(single_packed_len, par_packed);
        assert_eq!(single_col_len, par_col);
    }
}

#[test]
fn test_parallel_giant_chromosome_spanning_chunks() {
    // 1 giant chromosome (500 KB) followed by 10 short contigs
    let generator = FastxGen::fasta(1, 500_000).with_seed(42);
    let bytes = generator.generate_bytes();

    for threads in [2, 4, 8] {
        let count = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut c = 0;
                    for _ in parser.by_ref() {
                        c += 1;
                    }
                    c
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(
            count, 1,
            "Expected exactly 1 giant record across {threads} threads"
        );
    }
}

#[test]
fn test_parallel_file_parsing() {
    let bytes = FastxGen::illumina(500).with_seed(777).generate_bytes();
    let temp_dir = std::env::temp_dir();
    let path = temp_dir.join("test_helicase_parallel_file.fq");
    std::fs::write(&path, &bytes).unwrap();

    let mut single_len = 0;
    let mut single = FastxParser::<CONFIG_STR>::from_file(&path).unwrap();
    while single.next().is_some() {
        single_len += single.get_dna_string().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_len = ParallelFastxParser::<CONFIG_STR>::from_file_mmap(&path)
            .unwrap()
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_string().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(
            single_len, par_len,
            "Length mismatch on {threads} threads for file"
        );
    }

    let _ = std::fs::remove_file(path);
}

#[test]
fn test_parallel_multi_reader_file_parsing() {
    // 1. FASTQ test
    let fq_bytes = FastxGen::illumina(500).with_seed(888).generate_bytes();
    let temp_dir = std::env::temp_dir();
    let fq_path = temp_dir.join("test_helicase_parallel_reader.fq");
    std::fs::write(&fq_path, &fq_bytes).unwrap();

    let mut single_fq_len = 0;
    let mut single_fq = FastxParser::<CONFIG_STR>::from_file(&fq_path).unwrap();
    while single_fq.next().is_some() {
        single_fq_len += single_fq.get_dna_string().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_len = ParallelFastxParser::<CONFIG_STR>::from_file(&fq_path)
            .unwrap()
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_string().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(
            single_fq_len, par_len,
            "FASTQ length mismatch on {threads} threads for multi-reader"
        );
    }
    let _ = std::fs::remove_file(fq_path);

    // 2. FASTA wrapped test
    let fa_bytes = FastxGen::wrapped_fasta(50, 1000, 60)
        .with_seed(999)
        .generate_bytes();
    let fa_path = temp_dir.join("test_helicase_parallel_reader.fa");
    std::fs::write(&fa_path, &fa_bytes).unwrap();

    let mut single_fa_len = 0;
    let mut single_fa = FastxParser::<CONFIG_STR>::from_file(&fa_path).unwrap();
    while single_fa.next().is_some() {
        single_fa_len += single_fa.get_dna_string().len();
    }

    for threads in [1, 2, 4, 8] {
        let par_len = ParallelFastxParser::<CONFIG_STR>::from_file(&fa_path)
            .unwrap()
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(_) = parser.next() {
                        len += parser.get_dna_string().len();
                    }
                    len
                },
                |a, b| a + b,
            )
            .unwrap();

        assert_eq!(
            single_fa_len, par_len,
            "FASTA length mismatch on {threads} threads for multi-reader"
        );
    }
    let _ = std::fs::remove_file(fa_path);
}

#[test]
fn test_parallel_fastx_parser_builder() {
    let fq_bytes = FastxGen::illumina(200).with_seed(101).generate_bytes();
    let temp_dir = std::env::temp_dir();
    let path = temp_dir.join("test_builder.fq");
    std::fs::write(&path, &fq_bytes).unwrap();

    // 1. from_slice: reduce with manual loop
    let slice_bp = ParallelFastxParser::<CONFIG_STR>::from_slice(&fq_bytes)
        .threads(4)
        .reduce(
            |mut p| {
                let mut len = 0;
                while let Some(_) = p.next() {
                    len += p.get_dna_string().len();
                }
                len
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(slice_bp, 200 * 150);

    // 2. from_file (default multi-reader): reduce
    let file_bp = ParallelFastxParser::<CONFIG_STR>::from_file(&path)
        .unwrap()
        .threads(4)
        .reduce(
            |mut parser| {
                let mut len = 0;
                while let Some(_) = parser.next() {
                    len += parser.get_dna_string().len();
                }
                len
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(file_bp, slice_bp);

    let file_records = ParallelFastxParser::<CONFIG_STR>::from_file(&path)
        .unwrap()
        .threads(4)
        .reduce(
            |mut parser| {
                let mut c = 0;
                for _ in parser.by_ref() {
                    c += 1;
                }
                c
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(file_records, 200);

    // 3. from_file_mmap (zero-copy memory mapping)
    let mmap_bp = ParallelFastxParser::<CONFIG_STR>::from_file_mmap(&path)
        .unwrap()
        .threads(4)
        .reduce(
            |mut parser| {
                let mut len = 0;
                while let Some(_) = parser.next() {
                    len += parser.get_dna_string().len();
                }
                len
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(mmap_bp, slice_bp);

    // 4. from_reader (fallback streaming reader)
    let file = std::fs::File::open(&path).unwrap();
    let stream_bp = ParallelFastxParser::<CONFIG_STR>::from_reader(file)
        .threads(4)
        .reduce(
            |mut parser| {
                let mut len = 0;
                while let Some(_) = parser.next() {
                    len += parser.get_dna_string().len();
                }
                len
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(stream_bp, slice_bp);

    // 5. FastxParser::parallel_file convenience
    let conv_bp = FastxParser::<CONFIG_STR>::parallel_file(&path)
        .unwrap()
        .reduce(
            |mut parser| {
                let mut len = 0;
                while let Some(_) = parser.next() {
                    len += parser.get_dna_string().len();
                }
                len
            },
            |a, b| a + b,
        )
        .unwrap();
    assert_eq!(conv_bp, slice_bp);

    // 6. for_each helper without reduction (reduce on unit)
    use std::sync::atomic::{AtomicUsize, Ordering};
    let shared_counter = AtomicUsize::new(0);
    ParallelFastxParser::<CONFIG_STR>::from_file(&path)
        .unwrap()
        .threads(4)
        .for_each(|mut parser| {
            let mut local = 0;
            for _ in parser.by_ref() {
                local += 1;
            }
            shared_counter.fetch_add(local, Ordering::Relaxed);
        })
        .unwrap();
    assert_eq!(shared_counter.load(Ordering::Relaxed), 200);

    let _ = std::fs::remove_file(path);
}

#[test]
fn test_parallel_records_fewer_than_threads() {
    let temp_dir = std::env::temp_dir();

    // Test cases: (record_count, record_len, is_fastq)
    let cases = [
        (0, 0, false),      // 0 records FASTA
        (0, 0, true),       // 0 records FASTQ
        (1, 100, false),    // 1 short FASTA record (< 32 KiB)
        (1, 100, true),     // 1 short FASTQ record (< 32 KiB)
        (1, 50_000, false), // 1 large FASTA record (> 32 KiB)
        (1, 50_000, true),  // 1 large FASTQ record (> 32 KiB)
        (2, 40_000, false), // 2 large FASTA records (> 32 KiB)
        (2, 40_000, true),  // 2 large FASTQ records (> 32 KiB)
        (3, 150, false),    // 3 short FASTA records (< 32 KiB)
        (3, 150, true),     // 3 short FASTQ records (< 32 KiB)
    ];

    for (rec_count, rec_len, is_fastq) in cases {
        let bytes = if rec_count == 0 {
            Vec::new()
        } else if is_fastq {
            FastxGen::fastq(rec_count, rec_len)
                .with_seed(123)
                .generate_bytes()
        } else {
            FastxGen::fasta(rec_count, rec_len)
                .with_seed(123)
                .generate_bytes()
        };

        // Reference single-threaded counts
        let mut expected_records = 0;
        let mut expected_bases = 0;
        if !bytes.is_empty() {
            let mut single = FastxParser::<CONFIG_STR>::from_slice(&bytes).unwrap();
            while single.next().is_some() {
                expected_records += 1;
                expected_bases += single.get_dna_string().len();
            }
        }
        assert_eq!(expected_records, rec_count);

        let test_file = temp_dir.join(format!("test_fewer_{rec_count}_{rec_len}_{is_fastq}.tmp"));
        std::fs::write(&test_file, &bytes).unwrap();

        for thread_count in [2, 4, 8, 16] {
            // 1. from_slice
            let (slice_records, slice_bases) =
                ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
                    .threads(thread_count)
                    .reduce(
                        |mut parser| {
                            let mut recs = 0;
                            let mut bases = 0;
                            while parser.next().is_some() {
                                recs += 1;
                                bases += parser.get_dna_string().len();
                            }
                            (recs, bases)
                        },
                        |acc, val| (acc.0 + val.0, acc.1 + val.1),
                    )
                    .unwrap();
            assert_eq!(slice_records, expected_records);
            assert_eq!(slice_bases, expected_bases);

            // 2. from_file (multi-reader)
            let (reader_records, reader_bases) =
                ParallelFastxParser::<CONFIG_STR>::from_file(&test_file)
                    .unwrap()
                    .threads(thread_count)
                    .reduce(
                        |mut parser| {
                            let mut recs = 0;
                            let mut bases = 0;
                            while parser.next().is_some() {
                                recs += 1;
                                bases += parser.get_dna_string().len();
                            }
                            (recs, bases)
                        },
                        |acc, val| (acc.0 + val.0, acc.1 + val.1),
                    )
                    .unwrap();
            assert_eq!(reader_records, expected_records);
            assert_eq!(reader_bases, expected_bases);

            // 3. from_file_mmap
            if !bytes.is_empty() {
                let (mmap_records, mmap_bases) =
                    ParallelFastxParser::<CONFIG_STR>::from_file_mmap(&test_file)
                        .unwrap()
                        .threads(thread_count)
                        .reduce(
                            |mut parser| {
                                let mut recs = 0;
                                let mut bases = 0;
                                while parser.next().is_some() {
                                    recs += 1;
                                    bases += parser.get_dna_string().len();
                                }
                                (recs, bases)
                            },
                            |acc, val| (acc.0 + val.0, acc.1 + val.1),
                        )
                        .unwrap();
                assert_eq!(mmap_records, expected_records);
                assert_eq!(mmap_bases, expected_bases);
            }
        }

        let _ = std::fs::remove_file(test_file);
    }
}
