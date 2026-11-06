use std::hint::black_box;
use std::io::{BufReader, Cursor, Read};
use std::path::Path;
use std::time::Instant;

use helicase::HelicaseParser;
use helicase::config::advanced::DEFAULT_CONFIG;
use helicase::parser::FastxParser;

const ROUNDS: usize = 7;

fn benchmark_file(name: &str, path_str: &str) {
    let path = Path::new(path_str);
    if !path.exists() {
        println!("File {path_str} not found, skipping.");
        return;
    }
    let compressed = std::fs::read(path).unwrap();
    let comp_len = compressed.len();

    // Determine uncompressed length via zefiro or zstd
    let mut uncompressed = Vec::new();
    zefiro::Decoder::new()
        .decompress_into_vec(&compressed, &mut uncompressed)
        .unwrap();
    let uncomp_len = uncompressed.len();

    println!("==========================================================================");
    println!(
        "Benchmark: {name} (Compressed: {:.2} MB, Uncompressed: {:.2} MB, Ratio: {:.2}x)",
        comp_len as f64 / 1_000_000.0,
        uncomp_len as f64 / 1_000_000.0,
        uncomp_len as f64 / comp_len as f64
    );
    println!("--------------------------------------------------------------------------");

    // 1. One-shot In-Memory Decompression (Zstd C vs Zefiro)
    {
        let mut times_zstd = Vec::with_capacity(ROUNDS);
        let mut times_zefiro = Vec::with_capacity(ROUNDS);
        let mut out = vec![0u8; uncomp_len];

        for _ in 0..ROUNDS {
            // zstd one-shot
            let t0 = Instant::now();
            let n = zstd::bulk::Decompressor::new()
                .unwrap()
                .decompress_to_buffer(&compressed, &mut out)
                .unwrap();
            times_zstd.push(t0.elapsed().as_secs_f64());
            assert_eq!(n, uncomp_len);

            // zefiro one-shot (with decompress_into_vec + checksum bypass)
            let mut dec = zefiro::Decoder::new();
            dec.set_verify_checksum(false);
            let mut out_vec = Vec::with_capacity(uncomp_len);
            let t1 = Instant::now();
            let n = dec.decompress_into_vec(&compressed, &mut out_vec).unwrap();
            times_zefiro.push(t1.elapsed().as_secs_f64());
            assert_eq!(n, uncomp_len);
        }

        let best_zstd = times_zstd.iter().cloned().fold(f64::INFINITY, f64::min);
        let best_zefiro = times_zefiro.iter().cloned().fold(f64::INFINITY, f64::min);

        let gb_s_zstd = (uncomp_len as f64 / 1e9) / best_zstd;
        let gb_s_zefiro = (uncomp_len as f64 / 1e9) / best_zefiro;
        let delta = (gb_s_zefiro - gb_s_zstd) / gb_s_zstd * 100.0;

        println!(
            "  One-shot Decomp:  libzstd: {:>6.2} GB/s ({:.2} ms) | zefiro: {:>6.2} GB/s ({:.2} ms) | Delta: {:>+6.2}%",
            gb_s_zstd,
            best_zstd * 1000.0,
            gb_s_zefiro,
            best_zefiro * 1000.0,
            delta
        );
    }

    // 2. Streaming Decompression from RAM (Read in 256 KiB chunks)
    {
        let mut times_zstd = Vec::with_capacity(ROUNDS);
        let mut times_zefiro = Vec::with_capacity(ROUNDS);
        let mut buf = vec![0u8; 256 * 1024];

        for _ in 0..ROUNDS {
            // libzstd streaming
            let cursor = Cursor::new(&compressed);
            let mut reader = zstd::stream::read::Decoder::with_buffer(BufReader::with_capacity(
                64 * 1024,
                cursor,
            ))
            .unwrap();
            let t0 = Instant::now();
            let mut total = 0;
            loop {
                let n = reader.read(&mut buf).unwrap();
                if n == 0 {
                    break;
                }
                total += n;
            }
            times_zstd.push(t0.elapsed().as_secs_f64());
            assert_eq!(total, uncomp_len);

            // zefiro streaming
            let cursor = Cursor::new(&compressed);
            let mut dec = zefiro::Decoder::new();
            dec.set_verify_checksum(false);
            let mut reader =
                zefiro::Reader::with_decoder(dec, BufReader::with_capacity(64 * 1024, cursor));
            let t1 = Instant::now();
            let mut total = 0;
            loop {
                let n = reader.read(&mut buf).unwrap();
                if n == 0 {
                    break;
                }
                total += n;
            }
            times_zefiro.push(t1.elapsed().as_secs_f64());
            assert_eq!(total, uncomp_len);
        }

        let best_zstd = times_zstd.iter().cloned().fold(f64::INFINITY, f64::min);
        let best_zefiro = times_zefiro.iter().cloned().fold(f64::INFINITY, f64::min);

        let gb_s_zstd = (uncomp_len as f64 / 1e9) / best_zstd;
        let gb_s_zefiro = (uncomp_len as f64 / 1e9) / best_zefiro;
        let delta = (gb_s_zefiro - gb_s_zstd) / gb_s_zstd * 100.0;

        println!(
            "  Streaming Decomp: libzstd: {:>6.2} GB/s ({:.2} ms) | zefiro: {:>6.2} GB/s ({:.2} ms) | Delta: {:>+6.2}%",
            gb_s_zstd,
            best_zstd * 1000.0,
            gb_s_zefiro,
            best_zefiro * 1000.0,
            delta
        );
    }

    // 3. End-to-End Single-Threaded FastxParser (Decompression + SIMD parsing)
    {
        let mut times_zstd = Vec::with_capacity(ROUNDS);
        let mut times_zefiro = Vec::with_capacity(ROUNDS);

        for _ in 0..ROUNDS {
            // libzstd + Helicase SIMD
            let cursor = Cursor::new(&compressed);
            let reader = zstd::stream::read::Decoder::with_buffer(BufReader::with_capacity(
                64 * 1024,
                cursor,
            ))
            .unwrap();
            let mut parser = FastxParser::<DEFAULT_CONFIG>::from_reader(reader).unwrap();
            let t0 = Instant::now();
            let mut bases = 0;
            while let Some(e) = parser.next() {
                e.unwrap();
                bases += parser.get_dna_len();
            }
            times_zstd.push(t0.elapsed().as_secs_f64());
            black_box(bases);

            // zefiro + Helicase SIMD
            let cursor = Cursor::new(&compressed);
            let mut dec = zefiro::Decoder::new();
            dec.set_verify_checksum(false);
            let reader =
                zefiro::Reader::with_decoder(dec, BufReader::with_capacity(64 * 1024, cursor));
            let mut parser = FastxParser::<DEFAULT_CONFIG>::from_reader(reader).unwrap();
            let t1 = Instant::now();
            let mut bases = 0;
            while let Some(e) = parser.next() {
                e.unwrap();
                bases += parser.get_dna_len();
            }
            times_zefiro.push(t1.elapsed().as_secs_f64());
            black_box(bases);
        }

        let best_zstd = times_zstd.iter().cloned().fold(f64::INFINITY, f64::min);
        let best_zefiro = times_zefiro.iter().cloned().fold(f64::INFINITY, f64::min);

        let gb_s_zstd = (uncomp_len as f64 / 1e9) / best_zstd;
        let gb_s_zefiro = (uncomp_len as f64 / 1e9) / best_zefiro;
        let delta = (gb_s_zefiro - gb_s_zstd) / gb_s_zstd * 100.0;

        println!(
            "  End-to-End SIMD:  libzstd: {:>6.2} GB/s ({:.2} ms) | zefiro: {:>6.2} GB/s ({:.2} ms) | Delta: {:>+6.2}%",
            gb_s_zstd,
            best_zstd * 1000.0,
            gb_s_zefiro,
            best_zefiro * 1000.0,
            delta
        );
    }
}

fn benchmark_parallel_file(name: &str, path_str: &str) {
    let path = Path::new(path_str);
    if !path.exists() {
        println!("File {path_str} not found, skipping.");
        return;
    }
    let compressed = std::fs::read(path).unwrap();
    let file_len = compressed.len();
    let mut uncompressed = Vec::new();
    zefiro::Decoder::new()
        .decompress_into_vec(&compressed, &mut uncompressed)
        .unwrap();
    let uncomp_len = uncompressed.len();

    println!("==========================================================================");
    println!(
        "Parallel Benchmark: {name} (File: {:.2} MB, Uncompressed: {:.2} MB)",
        file_len as f64 / 1_000_000.0,
        uncomp_len as f64 / 1_000_000.0,
    );
    println!("--------------------------------------------------------------------------");

    for threads in [1, 2, 4, 8] {
        let mut times = Vec::with_capacity(ROUNDS);
        for _ in 0..ROUNDS {
            let t0 = Instant::now();
            let total_bases = helicase::ParallelFastxParser::<DEFAULT_CONFIG>::from_file(path)
                .unwrap()
                .threads(threads)
                .for_each_with(
                    || 0usize,
                    |bases, mut parser| {
                        while let Some(e) = parser.next() {
                            e.unwrap();
                            *bases += parser.get_dna_len();
                        }
                        Ok(())
                    },
                )
                .unwrap()
                .into_iter()
                .sum::<usize>();
            times.push(t0.elapsed().as_secs_f64());
            black_box(total_bases);
        }
        let best_time = times.iter().cloned().fold(f64::INFINITY, f64::min);
        let uncomp_gb_s = (uncomp_len as f64 / 1e9) / best_time;
        let file_gb_s = (file_len as f64 / 1e9) / best_time;
        println!(
            "  ParallelFastxParser ({threads:>2} threads): {:>6.2} GB/s uncomp seq ({:>5.2} GB/s file) in {:>6.2} ms",
            uncomp_gb_s,
            file_gb_s,
            best_time * 1000.0
        );
    }
}

fn benchmark_pzstd() {
    let path_str = "../zefiro/sample_long_reads.fa.zst";
    if !Path::new(path_str).exists() {
        return;
    }
    let compressed = std::fs::read(path_str).unwrap();
    let mut uncompressed = Vec::new();
    zefiro::Decoder::new()
        .decompress_into_vec(&compressed, &mut uncompressed)
        .unwrap();
    let uncomp_len = uncompressed.len();

    // Create a multi-frame pzstd file (1 MB frames)
    let pzstd_path = "/tmp/benchmark_test.pzst";
    let chunk_size = 1024 * 1024;
    let mut pzstd_data = Vec::new();
    for chunk in uncompressed.chunks(chunk_size) {
        let comp = zstd::bulk::compress(chunk, 1).unwrap();
        pzstd_data.extend_from_slice(&comp);
    }
    std::fs::write(pzstd_path, &pzstd_data).unwrap();
    let pzstd_file_len = pzstd_data.len();

    println!("==========================================================================");
    println!(
        "Benchmark: Parallel pzstd Multi-Frame Decompression (File: {:.2} MB, Uncompressed: {:.2} MB, Frames: {})",
        pzstd_file_len as f64 / 1_000_000.0,
        uncomp_len as f64 / 1_000_000.0,
        uncompressed.chunks(chunk_size).count()
    );
    println!("--------------------------------------------------------------------------");

    for threads in [1, 2, 4, 8] {
        let mut times = Vec::with_capacity(ROUNDS);
        for _ in 0..ROUNDS {
            let t0 = Instant::now();
            let total_bases =
                helicase::ParallelFastxParser::<DEFAULT_CONFIG>::from_file(pzstd_path)
                    .unwrap()
                    .threads(threads)
                    .for_each_with(
                        || 0usize,
                        |bases, mut parser| {
                            while let Some(e) = parser.next() {
                                e.unwrap();
                                *bases += parser.get_dna_len();
                            }
                            Ok(())
                        },
                    )
                    .unwrap()
                    .into_iter()
                    .sum::<usize>();
            times.push(t0.elapsed().as_secs_f64());
            black_box(total_bases);
        }
        let best_time = times.iter().cloned().fold(f64::INFINITY, f64::min);
        let uncomp_gb_s = (uncomp_len as f64 / 1e9) / best_time;
        let file_gb_s = (pzstd_file_len as f64 / 1e9) / best_time;
        println!(
            "  ParallelFastxParser ({threads:>2} threads): {:>6.2} GB/s uncomp seq ({:>5.2} GB/s file) in {:>6.2} ms",
            uncomp_gb_s,
            file_gb_s,
            best_time * 1000.0
        );
    }
    let _ = std::fs::remove_file(pzstd_path);
}

fn main() {
    println!(
        "Comparing libzstd (v1.5.7 C) vs zefiro (pure-Rust) - Best of {ROUNDS} interleaved rounds"
    );
    benchmark_file("Short Reads FASTQ", "../zefiro/sample_short_reads.fq.zst");
    benchmark_file("Short Reads FASTA", "../zefiro/sample_short_reads.fa.zst");
    benchmark_file("Long Reads FASTQ", "../zefiro/sample_long_reads.fq.zst");
    benchmark_file("Long Reads FASTA", "../zefiro/sample_long_reads.fa.zst");

    println!("\n=== Parallel Single-Frame zstd Decompression ===");
    benchmark_parallel_file(
        "Short Reads FASTQ (single-frame .zst)",
        "../zefiro/sample_short_reads.fq.zst",
    );
    benchmark_parallel_file(
        "Long Reads FASTA (single-frame .zst)",
        "../zefiro/sample_long_reads.fa.zst",
    );

    println!("\n=== Parallel Multi-Frame pzstd Decompression (17 frames, 8 MB chunks) ===");
    benchmark_parallel_file(
        "Short Reads FASTQ (pzstd 17 frames)",
        "/tmp/sample_short_reads.fq.pzst",
    );
    benchmark_parallel_file(
        "Long Reads FASTA (pzstd 17 frames)",
        "/tmp/sample_long_reads.fa.pzst",
    );

    println!(
        "\n=== Parallel Synthetic pzstd Multi-Frame Decompression (134 frames, 1 MB chunks) ==="
    );
    benchmark_pzstd();
}
