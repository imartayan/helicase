use std::io::Read;
use std::path::{Path, PathBuf};

use crate::legacy::{LegacyArgs, run_legacy_bench};
use crate::measurement::{BaseTime, Measurement};
use crate::parallel::{
    print_parallel_io_markdown_table, print_parallel_io_summary_table,
    print_parallel_ram_markdown_table, print_parallel_ram_summary_table, run_parallel_benchmarks,
    run_parallel_file_benchmarks, run_parallel_mmap_benchmarks,
};
use crate::synthetic::{
    BenchOptions, OutputFormat, generate_datasets, run_custom_file_bench, run_synthetic_bench,
};
use clap::Parser as ClapParser;
use colored::Colorize;

#[derive(ClapParser, Debug)]
#[command(
    name = "bench",
    author,
    version,
    about = "Helicase Benchmark Suite",
    long_about = "High-performance FASTX parsing and bitpacking benchmark suite comparing Helicase against Needletail and Paraseq across memory, streaming I/O, mmap, and parallel backends."
)]
pub struct Args {
    /// Input file (FASTA/FASTQ) to benchmark. If omitted (or if -s/--synth is passed), runs synthetic datasets.
    #[arg(value_name = "FILE")]
    pub input: Option<PathBuf>,

    /// Force running synthetic benchmarks even if an input file is provided
    #[arg(short = 's', long, hide = true)]
    pub synth: bool,

    /// Run legacy v0.2 file benchmark mode (use `bench legacy --help`)
    #[arg(short = 'l', long, hide = true)]
    pub legacy: bool,

    /// Number of benchmark repetitions per test
    #[arg(short = 'r', long, default_value_t = 5)]
    pub repeat: u64,

    /// Number of warmup iterations before benchmarking
    #[arg(short = 'w', long, default_value_t = 1)]
    pub warmup: u64,

    /// Benchmark in-memory RAM slices/readers (enabled by default)
    #[arg(short = 'R', long, hide = true)]
    pub ram: bool,

    /// Add streaming FileInput benchmark (with 128 KiB buffer)
    #[arg(short = 'F', long)]
    pub file: bool,

    /// Add memory-mapped MmapInput benchmark
    #[arg(short = 'M', long)]
    pub mmap: bool,

    /// Disable in-memory RAM benchmarks (only run streaming File or Mmap)
    #[arg(long, hide = true)]
    pub no_ram: bool,

    /// Run multi-threaded parallel benchmarks (comma-separated thread counts, e.g. "1,2,4,8" or "4")
    #[arg(
        short = 'p',
        long = "parallel",
        num_args = 0..=1,
        default_missing_value = "4",
        value_name = "THREADS"
    )]
    pub parallel: Option<String>,

    /// Detailed bitpacking comparison (Needletail+PackedSeqVec vs Helicase str+PackedSeqVec vs native 2-bit)
    #[arg(short = 'b', long, alias = "bitpack")]
    pub pack_compare: bool,

    /// Output format: text (colored scoreboard), markdown (for docs/AGENTS.md), csv, json
    #[arg(short = 'f', long, default_value = "text")]
    pub format: String,

    /// Output in CSV format (shortcut for -f csv, hidden by default)
    #[arg(short = 'c', long, hide = true)]
    pub csv: bool,

    /// Suppress CSV header in legacy mode
    #[arg(short = 'H', long, hide = true)]
    pub no_csv_header: bool,

    /// Disable baseline parsers (Needletail & Paraseq)
    #[arg(short = 'B', long)]
    pub no_baseline: bool,

    /// Disable slice parsing in legacy mode
    #[arg(short = 'S', long, hide = true)]
    pub no_slice: bool,

    /// Disable Linux hardware perf metrics in legacy mode
    #[arg(short = 'P', long, hide = true)]
    pub no_perf: bool,

    /// Show result values (length, record count, verification)
    #[arg(short = 'v', long, hide = true)]
    pub show_val: bool,

    /// Compatibility flag passed automatically by cargo bench
    #[arg(long, hide = true)]
    pub bench: bool,

    /// Compatibility flag passed automatically by cargo test
    #[arg(long, hide = true)]
    pub test: bool,
}

fn parse_thread_counts(s: &str) -> Vec<usize> {
    s.split(',')
        .filter_map(|part| part.trim().parse::<usize>().ok())
        .filter(|&t| t > 0)
        .collect()
}

/// Detects whether a file is compressed by peeking its first 4 magic bytes
/// (gzip/bgzf, zstd, bz2, xz) and falling back to its file extension.
pub fn is_compressed_file(path: impl AsRef<Path>) -> bool {
    let path = path.as_ref();
    if let Ok(mut file) = std::fs::File::open(path) {
        let mut magic = [0u8; 4];
        if let Ok(n) = file.read(&mut magic) {
            // Gzip / BGZF: 0x1F, 0x8B
            if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
                return true;
            }
            // Zstandard: 0xFD2FB528 (little-endian: 0x28, 0xb5, 0x2f, 0xfd)
            if n >= 4 && magic == [0x28, 0xb5, 0x2f, 0xfd] {
                return true;
            }
            // Bzip2: 'B', 'Z'
            if n >= 2 && magic[0] == b'B' && magic[1] == b'Z' {
                return true;
            }
            // Xz: 0xFD, '7', 'z'
            if n >= 3 && magic[0] == 0xfd && magic[1] == b'7' && magic[2] == b'z' {
                return true;
            }
        }
    }
    path.extension().is_some_and(|ext| {
        ext == "gz"
            || ext == "bgz"
            || ext == "bgzf"
            || ext == "gzip"
            || ext == "zst"
            || ext == "zstd"
            || ext == "bz2"
            || ext == "xz"
    })
}

fn run_with_measurement<M: Measurement>(args: &Args) {
    if args.test {
        return;
    }

    let output_format = if args.csv {
        OutputFormat::Csv
    } else {
        OutputFormat::from_str_lossy(&args.format)
    };

    // 1. Legacy mode (v0.2.0 compatibility)
    if args.legacy {
        let Some(p) = &args.input else {
            eprintln!(
                "{} Legacy mode requires an input file path: cargo bench -- -l <path>",
                "Error:".red().bold()
            );
            std::process::exit(1);
        };
        let path = p.to_str().expect("Invalid file path");

        let legacy_args = LegacyArgs {
            repeat: args.repeat,
            warmup: args.warmup,
            csv: args.csv,
            no_csv_header: args.no_csv_header,
            no_baseline: args.no_baseline,
            no_slice: args.no_slice,
            file: args.file,
            mmap: args.mmap,
            show_val: args.show_val,
        };

        run_legacy_bench::<M>(&legacy_args, path);
        return;
    }

    // 2. Parallel scaling mode (Helicase vs Paraseq)
    if let Some(thread_str) = &args.parallel {
        let threads = parse_thread_counts(thread_str);
        let threads = if threads.is_empty() { vec![4] } else { threads };
        let is_markdown = output_format == OutputFormat::Markdown;

        let bench_file = args.file;
        let bench_mmap = args.mmap;

        let mut ram_summaries = Vec::new();
        let mut io_summaries = Vec::new();

        if let Some(path) = &args.input
            && !args.synth
        {
            let name = path.file_name().unwrap_or_default().to_string_lossy();
            let is_compressed = is_compressed_file(path);

            // If neither -F nor -M is specified:
            // - compressed files default to FileInput (-F)
            // - uncompressed files default to in-memory slice (RAM)
            let (run_file, run_mmap, run_ram) = if !bench_file && !bench_mmap {
                if is_compressed {
                    (true, false, false)
                } else {
                    (false, false, true)
                }
            } else {
                (bench_file, bench_mmap, false)
            };

            if run_file {
                let sum = run_parallel_file_benchmarks(
                    path,
                    &name,
                    &threads,
                    args.repeat as usize,
                    args.warmup as usize,
                    is_markdown,
                    args.no_baseline,
                );
                io_summaries.push(sum);
            }

            if run_mmap
                && let Some(sum) = run_parallel_mmap_benchmarks(
                    path,
                    &name,
                    &threads,
                    args.repeat as usize,
                    args.warmup as usize,
                    is_markdown,
                    args.no_baseline,
                )
            {
                io_summaries.push(sum);
            }

            if run_ram {
                let data = match std::fs::read(path) {
                    Ok(d) => d,
                    Err(e) => {
                        eprintln!("Error reading {}: {}", path.display(), e);
                        std::process::exit(1);
                    }
                };
                let sum = run_parallel_benchmarks(
                    &data,
                    &name,
                    &threads,
                    args.repeat as usize,
                    args.warmup as usize,
                    is_markdown,
                    args.no_baseline,
                );
                ram_summaries.push(sum);
            }
        } else {
            // Synthetic datasets mode with -p
            let run_file = bench_file;
            let run_mmap = bench_mmap;
            let run_ram = !args.no_ram && (args.ram || (!run_file && !run_mmap));

            let datasets = generate_datasets();

            if run_ram {
                for ds in &datasets {
                    let sum = run_parallel_benchmarks(
                        &ds.data,
                        ds.name,
                        &threads,
                        args.repeat as usize,
                        args.warmup as usize,
                        is_markdown,
                        args.no_baseline,
                    );
                    ram_summaries.push(sum);
                }
            }

            if run_file || run_mmap {
                use std::io::Write;
                let temp_dir = std::env::temp_dir().join("helicase_parallel_bench_temp");
                let _ = std::fs::create_dir_all(&temp_dir);

                for (idx, ds) in datasets.iter().enumerate() {
                    let file_ext = if ds.name.starts_with("FASTQ") {
                        "fastq"
                    } else {
                        "fasta"
                    };
                    let temp_file_path = temp_dir.join(format!("ds_{idx}.{file_ext}"));
                    if let Ok(mut f) = std::fs::File::create(&temp_file_path) {
                        let _ = f.write_all(&ds.data);
                    }

                    if run_file {
                        let sum = run_parallel_file_benchmarks(
                            &temp_file_path,
                            ds.name,
                            &threads,
                            args.repeat as usize,
                            args.warmup as usize,
                            is_markdown,
                            args.no_baseline,
                        );
                        io_summaries.push(sum);
                    }

                    if run_mmap
                        && let Some(sum) = run_parallel_mmap_benchmarks(
                            &temp_file_path,
                            ds.name,
                            &threads,
                            args.repeat as usize,
                            args.warmup as usize,
                            is_markdown,
                            args.no_baseline,
                        )
                    {
                        io_summaries.push(sum);
                    }

                    let _ = std::fs::remove_file(&temp_file_path);
                }

                let _ = std::fs::remove_dir_all(&temp_dir);
            }
        }

        if is_markdown {
            if !ram_summaries.is_empty() {
                print_parallel_ram_markdown_table(&ram_summaries, !args.no_baseline);
            }
            if !io_summaries.is_empty() {
                print_parallel_io_markdown_table(&io_summaries, !args.no_baseline);
            }
        } else {
            if !ram_summaries.is_empty() {
                print_parallel_ram_summary_table(&ram_summaries, !args.no_baseline);
            }
            if !io_summaries.is_empty() {
                print_parallel_io_summary_table(&io_summaries, !args.no_baseline);
            }
        }
        return;
    }

    // Additive I/O backends: RAM is default; -F and -M add their backends
    let mut bench_ram = !args.no_ram;
    let mut bench_file = args.file;
    let bench_mmap = args.mmap;

    if let Some(path) = &args.input
        && !args.synth
    {
        let is_compressed = is_compressed_file(path);
        if is_compressed {
            // Compressed files cannot be sliced in RAM directly; disable RAM and default to streaming FileInput
            bench_ram = false;
            if !bench_file && !bench_mmap {
                bench_file = true;
            }
        }
    }

    let opts = BenchOptions {
        repeat: args.repeat,
        warmup: args.warmup,
        no_baseline: args.no_baseline,
        bench_ram,
        bench_file,
        bench_mmap,
        pack_compare: args.pack_compare,
        show_val: args.show_val,
        format: output_format,
    };

    // 3. User-specified custom file
    if let Some(path) = &args.input
        && !args.synth
    {
        run_custom_file_bench::<M>(path, &opts);
        return;
    }

    // 4. Default: Standard synthetic suite
    run_synthetic_bench::<M>(&opts);
}

fn print_legacy_help() {
    println!(
        "Helicase Benchmark Suite — Legacy v0.2.0 Mode\n\n\
         Usage: bench -l [OPTIONS] <FILE>\n       \
                bench legacy [OPTIONS] <FILE>\n\n\
         Arguments:\n  \
           <FILE>                 Input FASTA or FASTQ file to benchmark\n\n\
         Options:\n  \
           -r, --repeat <REPEAT>  Number of benchmark repetitions per test [default: 5]\n  \
           -w, --warmup <WARMUP>  Number of warmup iterations before benchmarking [default: 1]\n  \
           -c, --csv              Output in CSV format\n  \
           -H, --no-csv-header    Suppress CSV header\n  \
           -B, --no-baseline      Disable baseline parsers (Needletail & Paraseq)\n  \
           -S, --no-slice         Disable slice (in-memory) parsing\n  \
           -F, --file             Enable streaming FileInput benchmark\n  \
           -M, --mmap             Enable memory-mapped MmapInput benchmark\n  \
           -P, --no-perf          Disable Linux hardware perf metrics (instructions, branches)\n  \
           -v, --show-val         Show result values (length, record count, verification)\n  \
           -h, --help             Print help"
    );
}

fn is_thread_spec(s: &str) -> bool {
    !s.is_empty()
        && s.split(',')
            .all(|part| part.trim().parse::<usize>().is_ok_and(|n| n > 0))
}

pub fn normalize_args(raw_args: &[String]) -> Vec<String> {
    let mut normalized = Vec::with_capacity(raw_args.len());
    let mut i = 0;
    while i < raw_args.len() {
        let arg = &raw_args[i];

        // 1. Program name (index 0)
        if i == 0 {
            normalized.push(arg.clone());
            i += 1;
            continue;
        }

        // 2. Legacy alias: `bench legacy <FILE>` -> `bench -l <FILE>`
        if arg == "legacy" {
            normalized.push("-l".to_string());
            i += 1;
            continue;
        }

        // 3. Flags with equals: e.g. -p=8, --parallel=1,2,4,8
        if let Some(val) = arg
            .strip_prefix("-p=")
            .or_else(|| arg.strip_prefix("--parallel="))
        {
            normalized.push(format!("--parallel={val}"));
            i += 1;
            continue;
        }

        // 4. Standalone -p / --parallel
        if arg == "-p" || arg == "--parallel" {
            let next_is_thread = raw_args.get(i + 1).is_some_and(|next| is_thread_spec(next));

            if next_is_thread {
                let threads = &raw_args[i + 1];
                normalized.push(format!("--parallel={threads}"));
                i += 2;
            } else {
                normalized.push("--parallel=4".to_string());
                i += 1;
            }
            continue;
        }

        // 5. Short flag with thread spec directly attached: e.g. -p8, -p2,4,8
        if let Some(rem) = arg.strip_prefix("-p")
            && is_thread_spec(rem)
        {
            normalized.push(format!("--parallel={rem}"));
            i += 1;
            continue;
        }

        // 6. Clustered short flags containing 'p': e.g. -pFM, -pF, -pM, -FMp
        if arg.starts_with('-')
            && !arg.starts_with("--")
            && arg.len() > 2
            && arg.contains('p')
            && arg[1..].chars().all(|c| c.is_alphabetic())
        {
            let next_is_thread = raw_args.get(i + 1).is_some_and(|next| is_thread_spec(next));
            if next_is_thread {
                let threads = &raw_args[i + 1];
                normalized.push(format!("--parallel={threads}"));
                i += 2;
            } else {
                normalized.push("--parallel=4".to_string());
                i += 1;
            }
            for c in arg[1..].chars() {
                if c != 'p' {
                    normalized.push(format!("-{c}"));
                }
            }
            continue;
        }

        normalized.push(arg.clone());
        i += 1;
    }
    normalized
}

pub fn run() {
    let raw_args: Vec<String> = std::env::args().collect();

    // Check if legacy help was requested: `bench -l --help`, `bench --legacy -h`, or `bench legacy --help`
    let is_help = raw_args.iter().any(|a| a == "-h" || a == "--help");
    let is_legacy = raw_args
        .iter()
        .any(|a| a == "-l" || a == "--legacy" || a == "legacy");

    if is_help && is_legacy {
        print_legacy_help();
        return;
    }

    let normalized = normalize_args(&raw_args);
    let mut args = Args::parse_from(&normalized);

    // If an input path was provided as a relative path and cargo bench / cargo test
    // changed the current working directory to CARGO_MANIFEST_DIR (bench/), resolve
    // it against the user's shell working directory ($PWD).
    if let Some(ref mut path) = args.input
        && path.is_relative()
        && let Some(pwd) = std::env::var_os("PWD")
    {
        let candidate = PathBuf::from(pwd).join(&path);
        if candidate.exists() || !path.exists() {
            *path = candidate;
        }
    }

    #[cfg(target_os = "linux")]
    {
        if !args.no_perf && args.legacy {
            use crate::linux_perf::PerfMeasurement;
            run_with_measurement::<PerfMeasurement>(&args);
        } else {
            run_with_measurement::<BaseTime>(&args);
        }
    }

    #[cfg(not(target_os = "linux"))]
    {
        run_with_measurement::<BaseTime>(&args);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_string_vec(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_normalize_args_parallel() {
        // Standalone -p defaults to 4
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p"])),
            to_string_vec(&["bench", "--parallel=4"])
        );

        // -p with thread count
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p", "8"])),
            to_string_vec(&["bench", "--parallel=8"])
        );

        // -p with thread list
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p", "2,4,8"])),
            to_string_vec(&["bench", "--parallel=2,4,8"])
        );

        // -p with equals
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p=2,4,8"])),
            to_string_vec(&["bench", "--parallel=2,4,8"])
        );

        // --parallel with equals
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "--parallel=16"])),
            to_string_vec(&["bench", "--parallel=16"])
        );

        // -p followed by file path (should NOT consume file as thread count)
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p", "reads.fastq"])),
            to_string_vec(&["bench", "--parallel=4", "reads.fastq"])
        );

        // -p followed by threads and file path
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p", "2,4,8", "reads.fastq"])),
            to_string_vec(&["bench", "--parallel=2,4,8", "reads.fastq"])
        );

        // File path followed by -p
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "reads.fastq", "-p"])),
            to_string_vec(&["bench", "reads.fastq", "--parallel=4"])
        );

        // File path followed by -p 8
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "reads.fastq", "-p", "8"])),
            to_string_vec(&["bench", "reads.fastq", "--parallel=8"])
        );

        // Legacy mode rewrite
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "legacy", "reads.fastq"])),
            to_string_vec(&["bench", "-l", "reads.fastq"])
        );

        // Combined short flags with -p: -pFM
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-pFM"])),
            to_string_vec(&["bench", "--parallel=4", "-F", "-M"])
        );

        // Combined short flags with -p and thread count: -pFM 8
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-pFM", "8"])),
            to_string_vec(&["bench", "--parallel=8", "-F", "-M"])
        );

        // -p with attached thread count: -p8
        assert_eq!(
            normalize_args(&to_string_vec(&["bench", "-p8"])),
            to_string_vec(&["bench", "--parallel=8"])
        );
    }

    #[test]
    fn test_args_parsing_variants() {
        let parsed = Args::parse_from(normalize_args(&to_string_vec(&["bench"])));
        assert_eq!(parsed.parallel, None);
        assert_eq!(parsed.input, None);

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&["bench", "-p"])));
        assert_eq!(parsed.parallel.as_deref(), Some("4"));
        assert_eq!(parsed.input, None);

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&["bench", "-p", "2,4,8"])));
        assert_eq!(parsed.parallel.as_deref(), Some("2,4,8"));
        assert_eq!(parsed.input, None);

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&[
            "bench",
            "-p",
            "sample.fa",
        ])));
        assert_eq!(parsed.parallel.as_deref(), Some("4"));
        assert_eq!(
            parsed.input.as_ref().map(|p| p.to_str().unwrap()),
            Some("sample.fa")
        );

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&[
            "bench",
            "-p",
            "8",
            "sample.fa",
        ])));
        assert_eq!(parsed.parallel.as_deref(), Some("8"));
        assert_eq!(
            parsed.input.as_ref().map(|p| p.to_str().unwrap()),
            Some("sample.fa")
        );
    }

    #[test]
    fn test_warmup_arg() {
        let parsed = Args::parse_from(normalize_args(&to_string_vec(&["bench"])));
        assert_eq!(parsed.warmup, 1);
        assert_eq!(parsed.repeat, 5);

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&[
            "bench", "-w", "3", "-r", "10",
        ])));
        assert_eq!(parsed.warmup, 3);
        assert_eq!(parsed.repeat, 10);

        let parsed = Args::parse_from(normalize_args(&to_string_vec(&["bench", "--warmup=0"])));
        assert_eq!(parsed.warmup, 0);
    }
}
