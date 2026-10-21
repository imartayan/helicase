use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use colored::Colorize;
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use paraseq::fastx;
use paraseq::parallel::{ParallelProcessor, ParallelReader};

const DNA_STRING: Config = ParserOptions::default()
    .ignore_headers()
    .dna_string()
    .config();

const DNA_PACKED: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .skip_non_actg()
    .config();

#[derive(Default, Clone)]
struct ParaseqCounter {
    len: Arc<AtomicUsize>,
    records: Arc<AtomicUsize>,
}

impl<R: paraseq::Record> ParallelProcessor<R> for ParaseqCounter {
    fn process_record(&mut self, record: R) -> Result<(), paraseq::Error> {
        self.len.fetch_add(record.seq().len(), Ordering::Relaxed);
        self.records.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }
}

pub fn bench_helicase_single(data: &[u8], repeat: usize, warmup: usize) -> (f64, usize) {
    for _ in 0..warmup {
        let mut parser = FastxParser::<DNA_STRING>::from_slice(data).unwrap();
        let mut total_len = 0;
        while parser.next().is_some() {
            total_len += parser.get_dna_string().len();
        }
        std::hint::black_box(total_len);
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let mut parser = FastxParser::<DNA_STRING>::from_slice(data).unwrap();
        let mut total_len = 0;
        while parser.next().is_some() {
            total_len += parser.get_dna_string().len();
        }
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = total_len;
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (data.len() as f64 / 1e9) / min_time;
    (gb_s, last_len)
}

pub fn bench_helicase_parallel(
    data: &[u8],
    threads: usize,
    repeat: usize,
    warmup: usize,
) -> (f64, usize) {
    for _ in 0..warmup {
        let total_len = ParallelFastxParser::<DNA_STRING>::from_slice(data)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_string().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        std::hint::black_box(total_len);
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let total_len = ParallelFastxParser::<DNA_STRING>::from_slice(data)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_string().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = total_len;
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (data.len() as f64 / 1e9) / min_time;
    (gb_s, last_len)
}

pub fn bench_helicase_packed(
    data: &[u8],
    threads: usize,
    repeat: usize,
    warmup: usize,
) -> (f64, usize) {
    for _ in 0..warmup {
        let total_len = ParallelFastxParser::<DNA_PACKED>::from_slice(data)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_packed().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        std::hint::black_box(total_len);
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let total_len = ParallelFastxParser::<DNA_PACKED>::from_slice(data)
            .threads(threads)
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_packed().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = total_len;
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (data.len() as f64 / 1e9) / min_time;
    (gb_s, last_len)
}

pub fn bench_paraseq(
    data: &[u8],
    threads: usize,
    repeat: usize,
    warmup: usize,
) -> Result<(f64, usize), String> {
    for _ in 0..warmup {
        let reader = fastx::Reader::new(data).map_err(|e| e.to_string())?;
        let mut counter = ParaseqCounter::default();
        reader
            .process_parallel(&mut counter, threads)
            .map_err(|e| e.to_string())?;
        std::hint::black_box(counter.len.load(Ordering::Relaxed));
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let reader = fastx::Reader::new(data).map_err(|e| e.to_string())?;
        let mut counter = ParaseqCounter::default();
        let t0 = Instant::now();
        reader
            .process_parallel(&mut counter, threads)
            .map_err(|e| e.to_string())?;
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = counter.len.load(Ordering::Relaxed);
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (data.len() as f64 / 1e9) / min_time;
    Ok((gb_s, last_len))
}

pub fn bench_helicase_file_parallel(
    path: &Path,
    size: u64,
    threads: usize,
    repeat: usize,
    warmup: usize,
    use_mmap: bool,
) -> (f64, usize) {
    for _ in 0..warmup {
        let builder = if use_mmap {
            ParallelFastxParser::<DNA_STRING>::from_file_mmap(path)
        } else {
            ParallelFastxParser::<DNA_STRING>::from_file(path)
        }
        .unwrap()
        .threads(threads);
        let total_len = builder
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_string().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        std::hint::black_box(total_len);
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let builder = if use_mmap {
            ParallelFastxParser::<DNA_STRING>::from_file_mmap(path)
        } else {
            ParallelFastxParser::<DNA_STRING>::from_file(path)
        }
        .unwrap()
        .threads(threads);
        let total_len = builder
            .reduce(
                |mut parser| {
                    let mut len = 0;
                    while let Some(event) = parser.next() {
                        event?;
                        len += parser.get_dna_string().len();
                    }
                    Ok(len)
                },
                |a, b| a + b,
            )
            .unwrap();
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = total_len;
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (size as f64 / 1e9) / min_time;
    (gb_s, last_len)
}

pub fn bench_paraseq_file(
    path: &Path,
    size: u64,
    threads: usize,
    repeat: usize,
    warmup: usize,
) -> Result<(f64, usize), String> {
    for _ in 0..warmup {
        let reader = paraseq::ReaderBuilder::path(path)
            .build()
            .map_err(|e| e.to_string())?;
        let mut counter = ParaseqCounter::default();
        reader
            .process_parallel(&mut counter, threads)
            .map_err(|e| e.to_string())?;
        std::hint::black_box(counter.len.load(Ordering::Relaxed));
    }

    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let reader = paraseq::ReaderBuilder::path(path)
            .build()
            .map_err(|e| e.to_string())?;
        let mut counter = ParaseqCounter::default();
        let t0 = Instant::now();
        reader
            .process_parallel(&mut counter, threads)
            .map_err(|e| e.to_string())?;
        let elapsed = t0.elapsed().as_secs_f64();
        times.push(elapsed);
        last_len = counter.len.load(Ordering::Relaxed);
    }

    let min_time = times.iter().copied().fold(f64::INFINITY, f64::min);
    let gb_s = (size as f64 / 1e9) / min_time;
    Ok((gb_s, last_len))
}

#[derive(Debug, Clone)]
pub struct ParallelRamResult {
    pub threads: usize,
    pub heli_str_gb_s: f64,
    pub speedup: f64,
    pub heli_pack_gb_s: f64,
    pub paraseq_str_gb_s: Option<f64>,
    pub delta_paraseq: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParallelRamSummary {
    pub name: String,
    pub size_mb: f64,
    pub single_gb_s: f64,
    pub results: Vec<ParallelRamResult>,
}

#[derive(Debug, Clone)]
pub struct ParallelFileResult {
    pub threads: usize,
    pub file_gb_s: f64,
    pub seq_gb_s: f64,
    pub speedup: f64,
    pub paraseq_file_gb_s: Option<f64>,
    pub delta_paraseq: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParallelIoSummary {
    pub name: String,
    pub size_mb: f64,
    pub backend: &'static str,
    pub single_file_gb_s: f64,
    pub single_seq_gb_s: f64,
    pub results: Vec<ParallelFileResult>,
}

pub fn run_parallel_benchmarks(
    data: &[u8],
    name: &str,
    thread_counts: &[usize],
    repeat: usize,
    warmup: usize,
    markdown: bool,
    no_baseline: bool,
) -> ParallelRamSummary {
    let size_mb = data.len() as f64 / (1024.0 * 1024.0);
    let (single_gb_s, single_len) = bench_helicase_single(data, repeat, warmup);

    if markdown {
        println!("\n### Parallel Scaling: {} ({:.1} MB)", name, size_mb);
        if no_baseline {
            println!("| Threads | Helicase String | Speedup | Helicase Packed |");
            println!("|---|---|---|---|");
        } else {
            println!(
                "| Threads | Helicase String | Speedup | Helicase Packed | Paraseq String | vs Paraseq |"
            );
            println!("|---|---|---|---|---|---|");
        }
    } else {
        println!(
            "\n{} {} ({:.1} MB) — 1T single-thread: {:6.2} GB/s",
            "==> Parallel Benchmark:".bold().magenta(),
            name.bold(),
            size_mb,
            single_gb_s
        );
        if no_baseline {
            println!(
                "  {:<8} | {:<16} | {:<10} | {:<16}",
                "Threads", "Heli Str (GB/s)", "Speedup", "Heli Pack(GB/s)"
            );
            println!("  {}", "-".repeat(58));
        } else {
            println!(
                "  {:<8} | {:<16} | {:<10} | {:<16} | {:<16} | {:<16}",
                "Threads",
                "Heli Str (GB/s)",
                "Speedup",
                "Heli Pack(GB/s)",
                "Paraseq (GB/s)",
                "Delta vs Paraseq"
            );
            println!("  {}", "-".repeat(95));
        }
    }

    let mut results = Vec::with_capacity(thread_counts.len());

    for &threads in thread_counts {
        let (heli_gb_s, heli_len) = bench_helicase_parallel(data, threads, repeat, warmup);
        assert_eq!(
            single_len, heli_len,
            "Helicase length mismatch on {} threads!",
            threads
        );
        let (heli_pack_gb_s, _heli_pack_len) = bench_helicase_packed(data, threads, repeat, warmup);
        let speedup = heli_gb_s / single_gb_s;

        let (par_str, delta_str, delta_num, paraseq_str_gb_s, delta_paraseq) = if no_baseline {
            (String::new(), String::new(), String::new(), None, None)
        } else {
            let paraseq_res = bench_paraseq(data, threads, repeat, warmup);
            match paraseq_res {
                Ok((par_gb_s, _)) => {
                    let delta = ((heli_gb_s - par_gb_s) / par_gb_s) * 100.0;
                    let delta_colored = if delta >= 0.0 {
                        format!("+{:5.1}%", delta).green()
                    } else {
                        format!("{:5.1}%", delta).red()
                    };
                    (
                        format!("{:6.2} GB/s", par_gb_s),
                        delta_colored.to_string(),
                        format!("+{:5.1}%", delta),
                        Some(par_gb_s),
                        Some(delta),
                    )
                }
                Err(err) => (
                    format!("ERR: {}", err).dimmed().to_string(),
                    "N/A".dimmed().to_string(),
                    "N/A".to_string(),
                    None,
                    None,
                ),
            }
        };

        if no_baseline {
            if markdown {
                println!(
                    "| {:<7} | **{:6.2} GB/s** | {:5.2}x | {:6.2} GB/s |",
                    format!("{}T", threads),
                    heli_gb_s,
                    speedup,
                    heli_pack_gb_s,
                );
            } else {
                println!(
                    "  {:<8} | {:6.2} GB/s       | {:5.2}x     | {:6.2} GB/s",
                    format!("{}T", threads),
                    heli_gb_s,
                    speedup,
                    heli_pack_gb_s,
                );
            }
        } else {
            if markdown {
                println!(
                    "| {:<7} | **{:6.2} GB/s** | {:5.2}x | {:6.2} GB/s | {:<14} | {:<10} |",
                    format!("{}T", threads),
                    heli_gb_s,
                    speedup,
                    heli_pack_gb_s,
                    par_str,
                    delta_num
                );
            } else {
                println!(
                    "  {:<8} | {:6.2} GB/s       | {:5.2}x     | {:6.2} GB/s      | {:<16} | {:<16}",
                    format!("{}T", threads),
                    heli_gb_s,
                    speedup,
                    heli_pack_gb_s,
                    par_str,
                    delta_str
                );
            }
        }

        results.push(ParallelRamResult {
            threads,
            heli_str_gb_s: heli_gb_s,
            speedup,
            heli_pack_gb_s,
            paraseq_str_gb_s,
            delta_paraseq,
        });
    }

    ParallelRamSummary {
        name: name.to_string(),
        size_mb,
        single_gb_s,
        results,
    }
}

pub fn run_parallel_file_benchmarks(
    path: &Path,
    name: &str,
    thread_counts: &[usize],
    repeat: usize,
    warmup: usize,
    markdown: bool,
    no_baseline: bool,
) -> ParallelIoSummary {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let size_mb = file_size as f64 / (1024.0 * 1024.0);

    // Baseline: 1T single-threaded from_file
    let mut baseline_times = Vec::with_capacity(repeat);
    let mut single_len = 0;
    for _ in 0..warmup {
        let mut parser = FastxParser::<DNA_STRING>::from_file(path).unwrap();
        let mut len = 0;
        while parser.next().is_some() {
            len += parser.get_dna_string().len();
        }
        std::hint::black_box(len);
    }
    for _ in 0..repeat {
        let t0 = Instant::now();
        let mut parser = FastxParser::<DNA_STRING>::from_file(path).unwrap();
        let mut len = 0;
        while parser.next().is_some() {
            len += parser.get_dna_string().len();
        }
        baseline_times.push(t0.elapsed().as_secs_f64());
        single_len = len;
    }
    let min_baseline = baseline_times.iter().copied().fold(f64::INFINITY, f64::min);
    let single_file_gb_s = (file_size as f64 / 1e9) / min_baseline;
    let single_dna_gb_s = (single_len as f64 / 1e9) / min_baseline;

    if markdown {
        println!(
            "\n### Parallel File Decompression & Parsing: {} ({:.1} MB file, {:.1} Mbp seq)",
            name,
            size_mb,
            single_len as f64 / 1e6
        );
        if no_baseline {
            println!("| Threads | File Read (GB/s) | Uncomp Seq (GB/s) | Speedup vs 1T |");
            println!("|---|---|---|---|");
            println!(
                "| 1T (base) | {:6.2} GB/s | {:6.2} GB/s | 1.00x |",
                single_file_gb_s, single_dna_gb_s
            );
        } else {
            println!(
                "| Threads | File Read (GB/s) | Uncomp Seq (GB/s) | Speedup vs 1T | Paraseq (GB/s) | vs Paraseq |"
            );
            println!("|---|---|---|---|---|---|");
            println!(
                "| 1T (base) | {:6.2} GB/s | {:6.2} GB/s | 1.00x | - | - |",
                single_file_gb_s, single_dna_gb_s
            );
        }
    } else {
        println!(
            "\n{} {} ({:.1} MB file, {:.1} Mbp seq) — 1T baseline: {:6.2} GB/s file ({:6.2} GB/s seq)",
            "==> Parallel File Benchmark:".bold().magenta(),
            name.bold(),
            size_mb,
            single_len as f64 / 1e6,
            single_file_gb_s,
            single_dna_gb_s
        );
        if no_baseline {
            println!(
                "  {:<8} | {:<16} | {:<18} | {:<10}",
                "Threads", "File (GB/s)", "Uncomp Seq (GB/s)", "Speedup"
            );
            println!("  {}", "-".repeat(60));
        } else {
            println!(
                "  {:<8} | {:<16} | {:<18} | {:<10} | {:<16} | {:<16}",
                "Threads",
                "File (GB/s)",
                "Uncomp Seq (GB/s)",
                "Speedup",
                "Paraseq (GB/s)",
                "Delta vs Paraseq"
            );
            println!("  {}", "-".repeat(95));
        }
    }

    let mut results = Vec::with_capacity(thread_counts.len());

    for &threads in thread_counts {
        let mut times = Vec::with_capacity(repeat);
        let mut last_len = 0;
        for _ in 0..warmup {
            let total_len = ParallelFastxParser::<DNA_STRING>::from_file(path)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut len = 0;
                        while let Some(event) = parser.next() {
                            event?;
                            len += parser.get_dna_string().len();
                        }
                        Ok(len)
                    },
                    |a, b| a + b,
                )
                .unwrap();
            std::hint::black_box(total_len);
        }
        for _ in 0..repeat {
            let t0 = Instant::now();
            let total_len = ParallelFastxParser::<DNA_STRING>::from_file(path)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut len = 0;
                        while let Some(event) = parser.next() {
                            event?;
                            len += parser.get_dna_string().len();
                        }
                        Ok(len)
                    },
                    |a, b| a + b,
                )
                .unwrap();
            times.push(t0.elapsed().as_secs_f64());
            last_len = total_len;
        }
        assert_eq!(
            single_len, last_len,
            "Length mismatch on {threads} threads!"
        );
        let min_t = times.iter().copied().fold(f64::INFINITY, f64::min);
        let file_gb_s = (file_size as f64 / 1e9) / min_t;
        let seq_gb_s = (last_len as f64 / 1e9) / min_t;
        let speedup = file_gb_s / single_file_gb_s;

        let (par_str, delta_str, delta_num, paraseq_file_gb_s, delta_paraseq) = if no_baseline {
            (String::new(), String::new(), String::new(), None, None)
        } else {
            let paraseq_res = bench_paraseq_file(path, file_size, threads, repeat, warmup);
            match paraseq_res {
                Ok((par_gb_s, _)) => {
                    let delta = ((file_gb_s - par_gb_s) / par_gb_s) * 100.0;
                    let delta_colored = if delta >= 0.0 {
                        format!("+{:5.1}%", delta).green()
                    } else {
                        format!("{:5.1}%", delta).red()
                    };
                    (
                        format!("{:6.2} GB/s", par_gb_s),
                        delta_colored.to_string(),
                        format!("+{:5.1}%", delta),
                        Some(par_gb_s),
                        Some(delta),
                    )
                }
                Err(err) => (
                    format!("ERR: {}", err).dimmed().to_string(),
                    "N/A".dimmed().to_string(),
                    "N/A".to_string(),
                    None,
                    None,
                ),
            }
        };

        if no_baseline {
            if markdown {
                println!(
                    "| {:<7} | **{:6.2} GB/s** | **{:6.2} GB/s** | **{:5.2}x** |",
                    format!("{}T", threads),
                    file_gb_s,
                    seq_gb_s,
                    speedup
                );
            } else {
                let speedup_colored = if speedup >= 1.05 {
                    format!("{:5.2}x", speedup).green()
                } else {
                    format!("{:5.2}x", speedup).normal()
                };
                println!(
                    "  {:<8} | {:6.2} GB/s       | {:6.2} GB/s         | {}",
                    format!("{}T", threads),
                    file_gb_s,
                    seq_gb_s,
                    speedup_colored
                );
            }
        } else {
            if markdown {
                println!(
                    "| {:<7} | **{:6.2} GB/s** | **{:6.2} GB/s** | **{:5.2}x** | {:<14} | {:<10} |",
                    format!("{}T", threads),
                    file_gb_s,
                    seq_gb_s,
                    speedup,
                    par_str,
                    delta_num
                );
            } else {
                let speedup_colored = if speedup >= 1.05 {
                    format!("{:5.2}x", speedup).green()
                } else {
                    format!("{:5.2}x", speedup).normal()
                };
                println!(
                    "  {:<8} | {:6.2} GB/s       | {:6.2} GB/s         | {:<10} | {:<16} | {:<16}",
                    format!("{}T", threads),
                    file_gb_s,
                    seq_gb_s,
                    speedup_colored,
                    par_str,
                    delta_str
                );
            }
        }

        results.push(ParallelFileResult {
            threads,
            file_gb_s,
            seq_gb_s,
            speedup,
            paraseq_file_gb_s,
            delta_paraseq,
        });
    }

    ParallelIoSummary {
        name: name.to_string(),
        size_mb,
        backend: helicase::input::CompressionFormat::from_path(path)
            .ok()
            .flatten()
            .map(|c| c.as_str())
            .unwrap_or("FileInput"),
        single_file_gb_s,
        single_seq_gb_s: single_dna_gb_s,
        results,
    }
}

pub fn run_parallel_mmap_benchmarks(
    path: &Path,
    name: &str,
    thread_counts: &[usize],
    repeat: usize,
    warmup: usize,
    markdown: bool,
    _no_baseline: bool,
) -> Option<ParallelIoSummary> {
    let is_compressed = crate::cli::is_compressed_file(path);
    if is_compressed {
        eprintln!(
            "{}",
            "Warning: Memory-mapping (-M) is not supported for compressed files; skipping."
                .yellow()
        );
        return None;
    }

    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let size_mb = file_size as f64 / (1024.0 * 1024.0);

    // Baseline: 1T single-threaded from_file_mmap
    let mut baseline_times = Vec::with_capacity(repeat);
    let mut single_len = 0;
    for _ in 0..warmup {
        let mut parser = match FastxParser::<DNA_STRING>::from_file_mmap(path) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("Failed to mmap file {}: {}", path.display(), e);
                return None;
            }
        };
        let mut len = 0;
        while parser.next().is_some() {
            len += parser.get_dna_string().len();
        }
        std::hint::black_box(len);
    }
    for _ in 0..repeat {
        let t0 = Instant::now();
        let mut parser = FastxParser::<DNA_STRING>::from_file_mmap(path).unwrap();
        let mut len = 0;
        while parser.next().is_some() {
            len += parser.get_dna_string().len();
        }
        baseline_times.push(t0.elapsed().as_secs_f64());
        single_len = len;
    }
    let min_baseline = baseline_times.iter().copied().fold(f64::INFINITY, f64::min);
    let single_mmap_gb_s = (file_size as f64 / 1e9) / min_baseline;
    let single_dna_gb_s = (single_len as f64 / 1e9) / min_baseline;

    if markdown {
        println!(
            "\n### Parallel Memory-Mapped Parsing: {} ({:.1} MB file, {:.1} Mbp seq)",
            name,
            size_mb,
            single_len as f64 / 1e6
        );
        println!("| Threads | Mmap Read (GB/s) | Uncomp Seq (GB/s) | Speedup vs 1T |");
        println!("|---|---|---|---|");
        println!(
            "| 1T (base) | {:6.2} GB/s | {:6.2} GB/s | 1.00x |",
            single_mmap_gb_s, single_dna_gb_s
        );
    } else {
        println!(
            "\n{} {} ({:.1} MB file, {:.1} Mbp seq) — 1T baseline: {:6.2} GB/s mmap ({:6.2} GB/s seq)",
            "==> Parallel Mmap Benchmark:".bold().cyan(),
            name.bold(),
            size_mb,
            single_len as f64 / 1e6,
            single_mmap_gb_s,
            single_dna_gb_s
        );
        println!(
            "  {:<8} | {:<16} | {:<18} | {:<10}",
            "Threads", "Mmap (GB/s)", "Uncomp Seq (GB/s)", "Speedup"
        );
        println!("  {}", "-".repeat(60));
    }

    let mut results = Vec::with_capacity(thread_counts.len());

    for &threads in thread_counts {
        let mut times = Vec::with_capacity(repeat);
        let mut last_len = 0;
        for _ in 0..warmup {
            let total_len = ParallelFastxParser::<DNA_STRING>::from_file_mmap(path)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut len = 0;
                        while let Some(event) = parser.next() {
                            event?;
                            len += parser.get_dna_string().len();
                        }
                        Ok(len)
                    },
                    |a, b| a + b,
                )
                .unwrap();
            std::hint::black_box(total_len);
        }
        for _ in 0..repeat {
            let t0 = Instant::now();
            let total_len = ParallelFastxParser::<DNA_STRING>::from_file_mmap(path)
                .unwrap()
                .threads(threads)
                .reduce(
                    |mut parser| {
                        let mut len = 0;
                        while let Some(event) = parser.next() {
                            event?;
                            len += parser.get_dna_string().len();
                        }
                        Ok(len)
                    },
                    |a, b| a + b,
                )
                .unwrap();
            times.push(t0.elapsed().as_secs_f64());
            last_len = total_len;
        }
        assert_eq!(
            single_len, last_len,
            "Length mismatch on {threads} threads!"
        );
        let min_t = times.iter().copied().fold(f64::INFINITY, f64::min);
        let mmap_gb_s = (file_size as f64 / 1e9) / min_t;
        let seq_gb_s = (last_len as f64 / 1e9) / min_t;
        let speedup = mmap_gb_s / single_mmap_gb_s;

        if markdown {
            println!(
                "| {:<7} | **{:6.2} GB/s** | **{:6.2} GB/s** | **{:5.2}x** |",
                format!("{}T", threads),
                mmap_gb_s,
                seq_gb_s,
                speedup
            );
        } else {
            let speedup_colored = if speedup >= 1.05 {
                format!("{:5.2}x", speedup).green()
            } else {
                format!("{:5.2}x", speedup).normal()
            };
            println!(
                "  {:<8} | {:6.2} GB/s       | {:6.2} GB/s         | {}",
                format!("{}T", threads),
                mmap_gb_s,
                seq_gb_s,
                speedup_colored
            );
        }

        results.push(ParallelFileResult {
            threads,
            file_gb_s: mmap_gb_s,
            seq_gb_s,
            speedup,
            paraseq_file_gb_s: None,
            delta_paraseq: None,
        });
    }

    Some(ParallelIoSummary {
        name: name.to_string(),
        size_mb,
        backend: "MmapInput",
        single_file_gb_s: single_mmap_gb_s,
        single_seq_gb_s: single_dna_gb_s,
        results,
    })
}

pub fn print_parallel_ram_summary_table(summaries: &[ParallelRamSummary], has_baseline: bool) {
    if summaries.is_empty() {
        return;
    }
    println!("\n{}", "=".repeat(113).dimmed());
    println!(
        "{}",
        "PARALLEL BENCHMARK SCOREBOARD (RAM Slice) (Throughput in GB/s)".bold()
    );
    println!("{}", "-".repeat(113).dimmed());

    if has_baseline {
        println!(
            "{:<24} {:>8} | {:<7} | {:>10} {:>12} {:>8} {:>11} | {:>10} | {:>10}",
            "Dataset",
            "Size",
            "Threads",
            "1T Single",
            "Heli Str",
            "Speedup",
            "Heli Pack",
            "Paraseq",
            "vs Paraseq"
        );
    } else {
        println!(
            "{:<24} {:>8} | {:<7} | {:>10} {:>12} {:>8} {:>11}",
            "Dataset", "Size", "Threads", "1T Single", "Heli Str", "Speedup", "Heli Pack"
        );
    }
    println!("{}", "-".repeat(113).dimmed());

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        for (i, res) in s.results.iter().enumerate() {
            let (name_col, size_col, single_col) = if i == 0 {
                (
                    s.name.as_str(),
                    size_str.as_str(),
                    format!("{:5.2} GB/s", s.single_gb_s),
                )
            } else {
                ("", "", String::new())
            };

            let threads_str = format!("{}T", res.threads);
            let h_str = format!("{:5.2} GB/s", res.heli_str_gb_s);
            let speedup_str = format!("{:4.2}x", res.speedup);
            let h_pack = format!("{:5.2} GB/s", res.heli_pack_gb_s);

            if has_baseline {
                let par_str = res
                    .paraseq_str_gb_s
                    .map(|p| format!("{:5.2} GB/s", p))
                    .unwrap_or_else(|| "N/A".to_string());
                let delta_str = match res.delta_paraseq {
                    Some(d) if d >= 0.0 => format!("+{:5.1}%", d).green().bold().to_string(),
                    Some(d) => format!("{:5.1}%", d).red().bold().to_string(),
                    None => "N/A".dimmed().to_string(),
                };

                println!(
                    "{:<24} {:>8} | {:<7} | {:>10} {:>12} {:>8} {:>11} | {:>10} | {:>10}",
                    name_col,
                    size_col,
                    threads_str.bold(),
                    single_col,
                    h_str.cyan().bold(),
                    speedup_str,
                    h_pack.bright_green(),
                    par_str,
                    delta_str
                );
            } else {
                println!(
                    "{:<24} {:>8} | {:<7} | {:>10} {:>12} {:>8} {:>11}",
                    name_col,
                    size_col,
                    threads_str.bold(),
                    single_col,
                    h_str.cyan().bold(),
                    speedup_str,
                    h_pack.bright_green()
                );
            }
        }
    }
    println!("{}\n", "=".repeat(113).dimmed());
}

pub fn print_parallel_ram_markdown_table(summaries: &[ParallelRamSummary], has_baseline: bool) {
    if summaries.is_empty() {
        return;
    }
    println!("\n### Parallel Scaling Scoreboard (RAM Slice)");
    if has_baseline {
        println!(
            "| Dataset | Size | Threads | 1T Single | Helicase String | Speedup | Helicase Packed | Paraseq String | vs Paraseq |"
        );
        println!("|---|---|---|---|---|---|---|---|---|");
    } else {
        println!(
            "| Dataset | Size | Threads | 1T Single | Helicase String | Speedup | Helicase Packed |"
        );
        println!("|---|---|---|---|---|---|---|");
    }

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        for res in &s.results {
            let single_str = format!("{:5.2} GB/s", s.single_gb_s);
            let h_str = format!("**{:5.2} GB/s**", res.heli_str_gb_s);
            let speedup_str = format!("{:4.2}x", res.speedup);
            let h_pack = format!("{:5.2} GB/s", res.heli_pack_gb_s);

            if has_baseline {
                let par_str = res
                    .paraseq_str_gb_s
                    .map(|p| format!("{:5.2} GB/s", p))
                    .unwrap_or_else(|| "N/A".to_string());
                let delta_str = match res.delta_paraseq {
                    Some(d) if d >= 0.0 => format!("+{:5.1}%", d),
                    Some(d) => format!("{:5.1}%", d),
                    None => "N/A".to_string(),
                };
                println!(
                    "| {} | {} | {}T | {} | {} | {} | {} | {} | {} |",
                    s.name,
                    size_str,
                    res.threads,
                    single_str,
                    h_str,
                    speedup_str,
                    h_pack,
                    par_str,
                    delta_str
                );
            } else {
                println!(
                    "| {} | {} | {}T | {} | {} | {} | {} |",
                    s.name, size_str, res.threads, single_str, h_str, speedup_str, h_pack
                );
            }
        }
    }
}

pub fn print_parallel_io_summary_table(summaries: &[ParallelIoSummary], has_baseline: bool) {
    if summaries.is_empty() {
        return;
    }
    println!("\n{}", "=".repeat(125).dimmed());
    println!(
        "{}",
        "PARALLEL I/O SCOREBOARD (Streaming File & Memory-Mapped) (Throughput in GB/s)".bold()
    );
    println!("{}", "-".repeat(125).dimmed());

    if has_baseline {
        println!(
            "{:<24} {:>8} | {:<10} {:<7} | {:>10} {:>12} {:>12} {:>8} | {:>10} | {:>10}",
            "Dataset",
            "Size",
            "Backend",
            "Threads",
            "1T Base",
            "File(GB/s)",
            "Seq(GB/s)",
            "Speedup",
            "Paraseq",
            "vs Paraseq"
        );
    } else {
        println!(
            "{:<24} {:>8} | {:<10} {:<7} | {:>10} {:>12} {:>12} {:>8}",
            "Dataset",
            "Size",
            "Backend",
            "Threads",
            "1T Base",
            "File(GB/s)",
            "Seq(GB/s)",
            "Speedup"
        );
    }
    println!("{}", "-".repeat(125).dimmed());

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        for (i, res) in s.results.iter().enumerate() {
            let (name_col, size_col, backend_col, base_col) = if i == 0 {
                (
                    s.name.as_str(),
                    size_str.as_str(),
                    s.backend,
                    format!("{:5.2} GB/s", s.single_file_gb_s),
                )
            } else {
                ("", "", "", String::new())
            };

            let threads_str = format!("{}T", res.threads);
            let file_str = format!("{:5.2} GB/s", res.file_gb_s);
            let seq_str = format!("{:5.2} GB/s", res.seq_gb_s);
            let speedup_str = format!("{:4.2}x", res.speedup);

            if has_baseline {
                let par_str = res
                    .paraseq_file_gb_s
                    .map(|p| format!("{:5.2} GB/s", p))
                    .unwrap_or_else(|| "N/A".to_string());
                let delta_str = match res.delta_paraseq {
                    Some(d) if d >= 0.0 => format!("+{:5.1}%", d).green().bold().to_string(),
                    Some(d) => format!("{:5.1}%", d).red().bold().to_string(),
                    None => "N/A".dimmed().to_string(),
                };

                println!(
                    "{:<24} {:>8} | {:<10} {:<7} | {:>10} {:>12} {:>12} {:>8} | {:>10} | {:>10}",
                    name_col,
                    size_col,
                    backend_col,
                    threads_str.bold(),
                    base_col,
                    file_str.cyan().bold(),
                    seq_str.bright_green(),
                    speedup_str,
                    par_str,
                    delta_str
                );
            } else {
                println!(
                    "{:<24} {:>8} | {:<10} {:<7} | {:>10} {:>12} {:>12} {:>8}",
                    name_col,
                    size_col,
                    backend_col,
                    threads_str.bold(),
                    base_col,
                    file_str.cyan().bold(),
                    seq_str.bright_green(),
                    speedup_str
                );
            }
        }
    }
    println!("{}\n", "=".repeat(125).dimmed());
}

pub fn print_parallel_io_markdown_table(summaries: &[ParallelIoSummary], has_baseline: bool) {
    if summaries.is_empty() {
        return;
    }
    println!("\n### Parallel I/O Scoreboard (Streaming File & Memory-Mapped)");
    if has_baseline {
        println!(
            "| Dataset | Size | Backend | Threads | 1T Base | File Read (GB/s) | Uncomp Seq (GB/s) | Speedup vs 1T | Paraseq (GB/s) | vs Paraseq |"
        );
        println!("|---|---|---|---|---|---|---|---|---|---|");
    } else {
        println!(
            "| Dataset | Size | Backend | Threads | 1T Base | File Read (GB/s) | Uncomp Seq (GB/s) | Speedup vs 1T |"
        );
        println!("|---|---|---|---|---|---|---|---|");
    }

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        for res in &s.results {
            let base_str = format!("{:5.2} GB/s", s.single_file_gb_s);
            let file_str = format!("**{:5.2} GB/s**", res.file_gb_s);
            let seq_str = format!("{:5.2} GB/s", res.seq_gb_s);
            let speedup_str = format!("{:4.2}x", res.speedup);

            if has_baseline {
                let par_str = res
                    .paraseq_file_gb_s
                    .map(|p| format!("{:5.2} GB/s", p))
                    .unwrap_or_else(|| "N/A".to_string());
                let delta_str = match res.delta_paraseq {
                    Some(d) if d >= 0.0 => format!("+{:5.1}%", d),
                    Some(d) => format!("{:5.1}%", d),
                    None => "N/A".to_string(),
                };
                println!(
                    "| {} | {} | {} | {}T | {} | {} | {} | {} | {} | {} |",
                    s.name,
                    size_str,
                    s.backend,
                    res.threads,
                    base_str,
                    file_str,
                    seq_str,
                    speedup_str,
                    par_str,
                    delta_str
                );
            } else {
                println!(
                    "| {} | {} | {} | {}T | {} | {} | {} | {} |",
                    s.name,
                    size_str,
                    s.backend,
                    res.threads,
                    base_str,
                    file_str,
                    seq_str,
                    speedup_str
                );
            }
        }
    }
}
