#![cfg(feature = "std")]

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || (0, 0),
                |(recs, len), mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *recs += 1;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_packed().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

        let par_col = ParallelFastxParser::<CONFIG_COL>::from_slice(&bytes)
            .threads(threads)
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_columnar().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || 0,
                |c, mut parser| {
                    for event in parser.by_ref() {
                        event?;
                        *c += 1;
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
            .for_each_with(
                || 0,
                |len, mut parser| {
                    while let Some(event) = parser.next() {
                        event?;
                        *len += parser.get_dna_string().len();
                    }
                    Ok(())
                },
            )
            .unwrap()
            .into_iter()
            .sum::<usize>();

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
        .for_each_with(
            || 0,
            |len, mut p| {
                while let Some(event) = p.next() {
                    event?;
                    *len += p.get_dna_string().len();
                }
                Ok(())
            },
        )
        .unwrap()
        .into_iter()
        .sum::<usize>();
    assert_eq!(slice_bp, 200 * 150);

    // 2. from_file (default multi-reader): reduce
    let file_bp = ParallelFastxParser::<CONFIG_STR>::from_file(&path)
        .unwrap()
        .threads(4)
        .for_each_with(
            || 0,
            |len, mut parser| {
                while let Some(event) = parser.next() {
                    event?;
                    *len += parser.get_dna_string().len();
                }
                Ok(())
            },
        )
        .unwrap()
        .into_iter()
        .sum::<usize>();
    assert_eq!(file_bp, slice_bp);

    let file_records = ParallelFastxParser::<CONFIG_STR>::from_file(&path)
        .unwrap()
        .threads(4)
        .for_each_with(
            || 0,
            |c, mut parser| {
                for event in parser.by_ref() {
                    event?;
                    *c += 1;
                }
                Ok(())
            },
        )
        .unwrap()
        .into_iter()
        .sum::<usize>();
    assert_eq!(file_records, 200);

    // 3. from_file_mmap (zero-copy memory mapping)
    let mmap_bp = ParallelFastxParser::<CONFIG_STR>::from_file_mmap(&path)
        .unwrap()
        .threads(4)
        .for_each_with(
            || 0,
            |len, mut parser| {
                while let Some(event) = parser.next() {
                    event?;
                    *len += parser.get_dna_string().len();
                }
                Ok(())
            },
        )
        .unwrap()
        .into_iter()
        .sum::<usize>();
    assert_eq!(mmap_bp, slice_bp);

    // 4. from_reader (fallback streaming reader)
    let file = std::fs::File::open(&path).unwrap();
    let stream_bp = ParallelFastxParser::<CONFIG_STR>::from_reader(file)
        .threads(4)
        .for_each_with(
            || 0,
            |len, mut parser| {
                while let Some(event) = parser.next() {
                    event?;
                    *len += parser.get_dna_string().len();
                }
                Ok(())
            },
        )
        .unwrap()
        .into_iter()
        .sum::<usize>();
    assert_eq!(stream_bp, slice_bp);

    // 5. for_each helper without reduction (reduce on unit)
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
                    .for_each_with(
                        || (0, 0),
                        |(recs, bases), mut parser| {
                            while let Some(event) = parser.next() {
                                event?;
                                *recs += 1;
                                *bases += parser.get_dna_string().len();
                            }
                            Ok(())
                        },
                    )
                    .unwrap()
                    .into_iter()
                    .fold((0, 0), |acc, val| (acc.0 + val.0, acc.1 + val.1));
            assert_eq!(slice_records, expected_records);
            assert_eq!(slice_bases, expected_bases);

            // 2. from_file (multi-reader)
            let (reader_records, reader_bases) =
                ParallelFastxParser::<CONFIG_STR>::from_file(&test_file)
                    .unwrap()
                    .threads(thread_count)
                    .for_each_with(
                        || (0, 0),
                        |(recs, bases), mut parser| {
                            while let Some(event) = parser.next() {
                                event?;
                                *recs += 1;
                                *bases += parser.get_dna_string().len();
                            }
                            Ok(())
                        },
                    )
                    .unwrap()
                    .into_iter()
                    .fold((0, 0), |acc, val| (acc.0 + val.0, acc.1 + val.1));
            assert_eq!(reader_records, expected_records);
            assert_eq!(reader_bases, expected_bases);

            // 3. from_file_mmap
            if !bytes.is_empty() {
                let (mmap_records, mmap_bases) =
                    ParallelFastxParser::<CONFIG_STR>::from_file_mmap(&test_file)
                        .unwrap()
                        .threads(thread_count)
                        .for_each_with(
                            || (0, 0),
                            |(recs, bases), mut parser| {
                                while let Some(event) = parser.next() {
                                    event?;
                                    *recs += 1;
                                    *bases += parser.get_dna_string().len();
                                }
                                Ok(())
                            },
                        )
                        .unwrap()
                        .into_iter()
                        .fold((0, 0), |acc, val| (acc.0 + val.0, acc.1 + val.1));
                assert_eq!(mmap_records, expected_records);
                assert_eq!(mmap_bases, expected_bases);
            }
        }

        let _ = std::fs::remove_file(test_file);
    }
}

#[cfg(feature = "paraseq")]
#[test]
fn test_process_parallel_batch_callbacks() {
    use paraseq::Record;
    use paraseq::parallel::ParallelProcessor;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[derive(Clone, Default)]
    struct Counter {
        records: usize,
        records_flushed: Arc<AtomicUsize>,
        batches: Arc<AtomicUsize>,
    }

    impl<R: Record> ParallelProcessor<R> for Counter {
        fn process_record(&mut self, _: R) -> Result<(), paraseq::Error> {
            self.records += 1;
            Ok(())
        }
        fn on_batch_complete(&mut self) -> Result<(), paraseq::Error> {
            self.records_flushed
                .fetch_add(self.records, Ordering::Relaxed);
            self.batches.fetch_add(1, Ordering::Relaxed);
            self.records = 0;
            Ok(())
        }
    }

    // Chunks are capped at 16 MiB, so ~34 MB of input is at least 3 batches, even on one thread
    let bytes = FastxGen::fastq(40_000, 400).with_seed(42).generate_bytes();
    for threads in [1, 4] {
        let mut counter = Counter::default();
        ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .process_parallel(&mut counter)
            .unwrap();
        assert_eq!(counter.records_flushed.load(Ordering::Relaxed), 40_000);
        assert!(counter.batches.load(Ordering::Relaxed) >= 3);
    }
}

#[cfg(feature = "paraseq")]
#[test]
fn test_process_parallel_thread_ids() {
    use paraseq::Record;
    use paraseq::parallel::ParallelProcessor;
    use std::sync::{Arc, Mutex};

    #[derive(Clone, Default)]
    struct IdCollector {
        id: usize,
        ids: Arc<Mutex<Vec<usize>>>,
    }

    impl<R: Record> ParallelProcessor<R> for IdCollector {
        fn process_record(&mut self, _: R) -> Result<(), paraseq::Error> {
            Ok(())
        }
        fn on_thread_complete(&mut self) -> Result<(), paraseq::Error> {
            self.ids.lock().unwrap().push(self.id);
            Ok(())
        }
        fn set_thread_id(&mut self, thread_id: usize) {
            self.id = thread_id;
        }
    }

    let bytes = FastxGen::fastq(5000, 100).with_seed(42).generate_bytes();
    for threads in [1, 4] {
        let mut collector = IdCollector::default();
        ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .process_parallel(&mut collector)
            .unwrap();
        let mut ids = collector.ids.lock().unwrap().clone();
        ids.sort_unstable();
        let expected: Vec<usize> = (0..threads).collect();
        assert_eq!(
            ids, expected,
            "threads should be called once per worker thread"
        );
    }
}

#[cfg(feature = "paraseq")]
#[test]
fn test_process_parallel_ordered() {
    use paraseq::Record;
    use paraseq::parallel::{Ordered, ParallelProcessor};
    use std::sync::{Arc, Mutex};

    #[derive(Clone, Default)]
    struct OrderedCollector {
        local: Vec<String>,
        output: Arc<Mutex<Vec<String>>>,
        ordered: bool,
    }

    impl<R: Record> ParallelProcessor<R> for OrderedCollector {
        fn process_record(&mut self, record: R) -> Result<(), paraseq::Error> {
            self.local.push(record.id_str().to_string());
            Ok(())
        }
        fn on_batch_complete(&mut self) -> Result<(), paraseq::Error> {
            self.output.lock().unwrap().append(&mut self.local);
            Ok(())
        }
        fn requires_ordering(&self) -> bool {
            self.ordered
        }
    }

    // Many short records, and records longer than the chunks, which leaves empty chunks
    for bytes in [
        FastxGen::fastq(25_000, 120).with_seed(42).generate_bytes(),
        FastxGen::fasta(12, 300_000).with_seed(42).generate_bytes(),
    ] {
        // Ground truth: single-threaded sequential order
        let expected: Vec<String> = {
            let mut coll = OrderedCollector::default();
            ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
                .threads(1)
                .process_parallel(&mut coll)
                .unwrap();
            std::mem::take(&mut *coll.output.lock().unwrap())
        };

        // 1. Test via requires_ordering override across thread counts
        for threads in [2, 4, 8] {
            let mut coll = OrderedCollector {
                local: Vec::new(),
                output: Arc::new(Mutex::new(Vec::new())),
                ordered: true,
            };
            ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
                .threads(threads)
                .process_parallel(&mut coll)
                .unwrap();
            let output = std::mem::take(&mut *coll.output.lock().unwrap());
            assert_eq!(
                output, expected,
                "mismatch in ordered mode with {threads} threads"
            );
        }

        // 2. Test via paraseq::parallel::Ordered wrapper
        for threads in [2, 4] {
            let coll = OrderedCollector {
                local: Vec::new(),
                output: Arc::new(Mutex::new(Vec::new())),
                ordered: false, // inner returns false, Ordered wrapper overrides to true
            };
            let out_ref = Arc::clone(&coll.output);
            let mut ordered_proc = Ordered(coll);
            ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
                .threads(threads)
                .process_parallel(&mut ordered_proc)
                .unwrap();
            let output = std::mem::take(&mut *out_ref.lock().unwrap());
            assert_eq!(
                output, expected,
                "mismatch in Ordered(wrapper) mode with {threads} threads"
            );
        }
    }
}

#[cfg(feature = "paraseq")]
#[test]
fn test_process_parallel_ordered_error_unwinds() {
    use paraseq::Record;
    use paraseq::parallel::ParallelProcessor;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[derive(Clone, Default)]
    struct FailingProcessor {
        processed: std::sync::Arc<AtomicUsize>,
    }

    impl<R: Record> ParallelProcessor<R> for FailingProcessor {
        fn process_record(&mut self, _: R) -> Result<(), paraseq::Error> {
            let count = self.processed.fetch_add(1, Ordering::Relaxed);
            if count == 500 {
                return Err(paraseq::Error::Io(std::io::Error::other(
                    "synthetic failure",
                )));
            }
            Ok(())
        }
        fn requires_ordering(&self) -> bool {
            true
        }
    }

    let bytes = FastxGen::fastq(5000, 100).with_seed(42).generate_bytes();
    let mut proc = FailingProcessor::default();
    let res = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
        .threads(4)
        .process_parallel(&mut proc);
    assert!(res.is_err(), "error should propagate and not deadlock");
}

#[cfg(feature = "paraseq")]
#[test]
fn test_process_parallel_ordered_panic_and_thread_error() {
    use paraseq::Record;
    use paraseq::parallel::ParallelProcessor;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::mpsc::{RecvTimeoutError, channel};
    use std::time::Duration;

    #[derive(Clone)]
    struct Failing {
        panic: bool,
        records: Arc<AtomicUsize>,
    }

    impl<R: Record> ParallelProcessor<R> for Failing {
        fn process_record(&mut self, _: R) -> Result<(), paraseq::Error> {
            if self.panic && self.records.fetch_add(1, Ordering::Relaxed) == 500 {
                panic!("synthetic panic");
            }
            Ok(())
        }
        fn on_thread_complete(&mut self) -> Result<(), paraseq::Error> {
            Err(paraseq::Error::Io(std::io::Error::other(
                "synthetic failure",
            )))
        }
        fn requires_ordering(&self) -> bool {
            true
        }
    }

    let bytes = Arc::new(FastxGen::fastq(25_000, 120).with_seed(42).generate_bytes());
    for panic in [true, false] {
        let (done, finished) = channel::<()>();
        let bytes = Arc::clone(&bytes);
        std::thread::spawn(move || {
            let _done = done;
            let mut proc = Failing {
                panic,
                records: Arc::default(),
            };
            let res = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
                .threads(4)
                .process_parallel(&mut proc);
            assert!(res.is_err());
        });
        // The sender is dropped when the call returns or panics, a hang times out
        assert_eq!(
            finished.recv_timeout(Duration::from_secs(60)),
            Err(RecvTimeoutError::Disconnected),
            "ordered run hangs (panic: {panic})"
        );
    }
}

#[test]
fn test_parallel_for_each_with() {
    let bytes = FastxGen::fastq(20_000, 150).with_seed(7).generate_bytes();
    let path = std::env::temp_dir().join("test_helicase_parallel_for_each_with.fq");
    std::fs::write(&path, &bytes).unwrap();

    for threads in [1, 4] {
        let sources = [
            ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes),
            ParallelFastxParser::from_file(&path).unwrap().chunks(64),
            ParallelFastxParser::from_file_mmap(&path).unwrap(),
            ParallelFastxParser::from_reader(bytes.as_slice()),
        ];
        for parser in sources {
            let states = parser
                .threads(threads)
                .for_each_with(
                    || (0, 0),
                    |(records, len), mut parser| {
                        while let Some(event) = parser.next() {
                            event?;
                            *records += 1;
                            *len += parser.get_dna_string().len();
                        }
                        Ok(())
                    },
                )
                .unwrap();
            assert!((1..=threads).contains(&states.len()), "{threads} threads");
            let records: usize = states.iter().map(|s| s.0).sum();
            let len: usize = states.iter().map(|s| s.1).sum();
            assert_eq!((records, len), (20_000, 20_000 * 150), "{threads} threads");
        }
    }
    let empty = ParallelFastxParser::<CONFIG_STR>::from_slice(b"").for_each_with(|| 1, |_, _| ());
    assert_eq!(empty.unwrap(), vec![1]);
    let _ = std::fs::remove_file(path);
}

#[test]
fn test_parallel_custom_chunks() {
    use std::sync::atomic::{AtomicUsize, Ordering};

    let bytes = FastxGen::fastq(20_000, 150).with_seed(3).generate_bytes();
    let path = std::env::temp_dir().join("test_helicase_parallel_chunks.fq");
    std::fs::write(&path, &bytes).unwrap();
    let max_chunks = bytes.len() / (16 * 1024);

    for (chunks, expected) in [
        (None, 8),
        (Some(3), 3),
        (Some(50), 50),
        (Some(100_000), max_chunks),
    ] {
        let sources = [
            ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes),
            ParallelFastxParser::from_file(&path).unwrap(),
            ParallelFastxParser::from_file_mmap(&path).unwrap(),
        ];
        for parser in sources {
            let parser = parser.threads(2);
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
    }
    let _ = std::fs::remove_file(path);
}

#[test]
#[allow(deprecated)]
fn test_parallel_reduce_deprecated() {
    let bytes = FastxGen::fastq(20_000, 150).with_seed(5).generate_bytes();
    for threads in [1, 4] {
        let lens: Vec<usize> = ParallelFastxParser::<CONFIG_STR>::from_slice(&bytes)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut lens = Vec::new();
                    while let Some(event) = parser.next() {
                        event?;
                        lens.push(parser.get_dna_len());
                    }
                    Ok(lens)
                },
                |mut a, b| {
                    a.extend(b);
                    a
                },
            )
            .unwrap();
        assert_eq!(lens, vec![150; 20_000]);
    }
    let empty: usize = ParallelFastxParser::<CONFIG_STR>::from_slice(b"")
        .reduce(|_| Ok(1), |a, b| a + b)
        .unwrap();
    assert_eq!(empty, 1);
}
