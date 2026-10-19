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

pub fn bench_helicase_single(data: &[u8], repeat: usize) -> (f64, usize) {
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

pub fn bench_helicase_parallel(data: &[u8], threads: usize, repeat: usize) -> (f64, usize) {
    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let total_len = ParallelFastxParser::<DNA_STRING>::from_slice(data)
            .threads(threads)
            .reduce(
                || 0,
                |acc, mut parser| {
                    while parser.next().is_some() {
                        *acc += parser.get_dna_string().len();
                    }
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

pub fn bench_helicase_packed(data: &[u8], threads: usize, repeat: usize) -> (f64, usize) {
    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let total_len = ParallelFastxParser::<DNA_PACKED>::from_slice(data)
            .threads(threads)
            .reduce(
                || 0,
                |acc, mut parser| {
                    while parser.next().is_some() {
                        *acc += parser.get_dna_packed().len();
                    }
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

pub fn bench_paraseq(data: &[u8], threads: usize, repeat: usize) -> Result<(f64, usize), String> {
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
    use_mmap: bool,
) -> (f64, usize) {
    let mut last_len = 0;
    let mut times = Vec::with_capacity(repeat);
    for _ in 0..repeat {
        let t0 = Instant::now();
        let builder = ParallelFastxParser::<DNA_STRING>::from_file(path)
            .unwrap()
            .threads(threads);
        let total_len = if use_mmap {
            builder.mmap().reduce(
                || 0,
                |acc, mut parser| {
                    while parser.next().is_some() {
                        *acc += parser.get_dna_string().len();
                    }
                },
                |a, b| a + b,
            )
        } else {
            builder.reduce(
                || 0,
                |acc, mut parser| {
                    while parser.next().is_some() {
                        *acc += parser.get_dna_string().len();
                    }
                },
                |a, b| a + b,
            )
        }
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
) -> Result<(f64, usize), String> {
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

pub fn run_parallel_benchmarks(
    data: &[u8],
    name: &str,
    thread_counts: &[usize],
    repeat: usize,
    markdown: bool,
    no_baseline: bool,
) {
    let size_mb = data.len() as f64 / (1024.0 * 1024.0);
    let (single_gb_s, single_len) = bench_helicase_single(data, repeat);

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

    for &threads in thread_counts {
        let (heli_gb_s, heli_len) = bench_helicase_parallel(data, threads, repeat);
        assert_eq!(
            single_len, heli_len,
            "Helicase length mismatch on {} threads!",
            threads
        );
        let (heli_pack_gb_s, _heli_pack_len) = bench_helicase_packed(data, threads, repeat);
        let speedup = heli_gb_s / single_gb_s;

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
            let paraseq_res = bench_paraseq(data, threads, repeat);
            let (par_str, delta_str, delta_num) = match paraseq_res {
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
                    )
                }
                Err(err) => (
                    format!("ERR: {}", err).dimmed().to_string(),
                    "N/A".dimmed().to_string(),
                    "N/A".to_string(),
                ),
            };

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
    }
}
