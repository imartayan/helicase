use std::fs::File;
use std::io::Write;
use std::path::Path;

use colored::Colorize;
use fastx_gen::FastxGen;
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use needletail::{parse_fastx_file, parse_fastx_reader};
use packed_seq::{PackedSeqVec, Seq, SeqVec};

use crate::measurement::Measurement;

const DNA_STRING: Config = ParserOptions::default()
    .ignore_headers()
    .dna_string()
    .config();
const DNA_COLUMNAR: Config = ParserOptions::default()
    .ignore_headers()
    .dna_columnar()
    .skip_non_actg()
    .config();
const DNA_PACKED: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .skip_non_actg()
    .config();
const DNA_PACKED_KEEP: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .keep_non_actg()
    .config();

pub struct SynthDataset {
    pub name: &'static str,
    pub description: &'static str,
    pub data: Vec<u8>,
}

pub fn generate_datasets() -> Vec<SynthDataset> {
    vec![
        SynthDataset {
            name: "FASTQ 150bp (Illumina)",
            description: "40,000 short reads (150bp fixed)",
            data: FastxGen::illumina(40_000).with_seed(42).generate_bytes(),
        },
        SynthDataset {
            name: "FASTQ 15kb (PacBio)",
            description: "800 long reads (~15kb, 10% std_dev)",
            data: FastxGen::pacbio(800).with_seed(42).generate_bytes(),
        },
        SynthDataset {
            name: "FASTA wrapped 60-col",
            description: "100 contigs x 50kb (wrapped at 60 cols)",
            data: FastxGen::wrapped_fasta(100, 50_000, 60)
                .with_seed(42)
                .generate_bytes(),
        },
        SynthDataset {
            name: "FASTA wrapped 50-col",
            description: "100 contigs x 50kb (wrapped at 50 cols)",
            data: FastxGen::wrapped_fasta(100, 50_000, 50)
                .with_seed(42)
                .generate_bytes(),
        },
        SynthDataset {
            name: "FASTA unwrapped 250bp",
            description: "40,000 unwrapped contigs (250bp fixed)",
            data: FastxGen::fasta(40_000, 250).with_seed(42).generate_bytes(),
        },
        SynthDataset {
            name: "FASTA unwrapped 500kbp",
            description: "50 large contigs (500,000bp fixed)",
            data: FastxGen::fasta(50, 500_000).with_seed(42).generate_bytes(),
        },
    ]
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Text,
    Markdown,
    Csv,
    Json,
}

impl OutputFormat {
    pub fn from_str_lossy(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "markdown" | "md" => Self::Markdown,
            "csv" => Self::Csv,
            "json" => Self::Json,
            _ => Self::Text,
        }
    }
}

#[derive(Debug, Clone)]
pub struct BenchOptions {
    pub repeat: u64,
    pub no_baseline: bool,
    pub bench_ram: bool,
    pub bench_file: bool,
    pub bench_mmap: bool,
    pub pack_compare: bool,
    pub show_val: bool,
    pub format: OutputFormat,
}

impl Default for BenchOptions {
    fn default() -> Self {
        Self {
            repeat: 5,
            no_baseline: false,
            bench_ram: true,
            bench_file: false,
            bench_mmap: false,
            pack_compare: false,
            show_val: false,
            format: OutputFormat::Text,
        }
    }
}

#[derive(Default, Clone)]
pub struct DatasetSummary {
    pub name: String,
    pub size_mb: f64,
    pub needletail_str: Option<(f64, f64)>,
    pub needletail_packed: Option<(f64, f64)>,
    pub needletail_file: Option<(f64, f64)>,
    pub helicase_str: Option<(f64, f64)>,
    pub helicase_str_reader: Option<(f64, f64)>,
    pub helicase_str_packed: Option<(f64, f64)>,
    pub helicase_str_file: Option<(f64, f64)>,
    pub helicase_str_mmap: Option<(f64, f64)>,
    pub helicase_packed: Option<(f64, f64)>,
    pub helicase_packed_keep: Option<(f64, f64)>,
    pub helicase_packed_reader: Option<(f64, f64)>,
    pub helicase_packed_file: Option<(f64, f64)>,
    pub helicase_packed_mmap: Option<(f64, f64)>,
    pub helicase_col: Option<(f64, f64)>,
}

pub fn run_synthetic_bench<M: Measurement>(opts: &BenchOptions) {
    let is_human_text = opts.format == OutputFormat::Text;
    if is_human_text {
        println!(
            "{}",
            "Generating synthetic benchmark datasets in RAM..."
                .cyan()
                .bold()
        );
    }
    let datasets = generate_datasets();
    let mut summaries = Vec::new();

    let temp_dir = std::env::temp_dir().join("helicase_bench_temp");
    if opts.bench_file || opts.bench_mmap {
        let _ = std::fs::create_dir_all(&temp_dir);
    }

    for (idx, ds) in datasets.iter().enumerate() {
        let size = ds.data.len() as u64;
        let size_mb = size as f64 / 1_000_000.0;
        let slice = ds.data.as_slice();

        if is_human_text {
            println!(
                "\n{} {} ({:.1} MB, {})",
                "==>".bold().blue(),
                ds.name.bold(),
                size_mb,
                ds.description.dimmed()
            );
        }

        let mut summary = DatasetSummary {
            name: ds.name.to_string(),
            size_mb,
            ..Default::default()
        };

        // 1 warmup iteration to prime instruction/data caches
        {
            let mut parser = FastxParser::<DNA_STRING>::from_slice(slice).unwrap();
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_string().len();
            }
            std::hint::black_box(len);
        }

        let is_csv = opts.format == OutputFormat::Csv;
        let show_ticks = is_human_text || is_csv;
        let mut m = M::new(ds.name);
        let show_tick = |m: &M, label: &str, dna_len: usize| {
            if show_ticks {
                m.show(
                    label,
                    size,
                    if opts.show_val { Some(dna_len) } else { None },
                    is_csv,
                );
            }
        };

        // 1. Needletail baseline (slice reader)
        if !opts.no_baseline && opts.bench_ram {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut reader = parse_fastx_reader(slice).expect("Failed to parse slice");
                dna_len = 0;
                while let Some(r) = reader.next() {
                    let record = r.expect("Invalid record");
                    let seq = record.seq();
                    dna_len += seq.len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Needletail string \t(reader)", dna_len);
            summary.needletail_str = m.throughput(size);
        }

        // 2. Needletail string + PackedSeqVec (reused buffer)
        if !opts.no_baseline && (opts.bench_ram || opts.pack_compare) {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut reader = parse_fastx_reader(slice).expect("Failed to parse slice");
                let mut packed = PackedSeqVec::default();
                dna_len = 0;
                while let Some(r) = reader.next() {
                    let record = r.expect("Invalid record");
                    packed.clear();
                    packed.push_ascii(&record.seq());
                    dna_len += packed.as_slice().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Needletail+PackedSeqVec", dna_len);
            summary.needletail_packed = m.throughput(size);
        }

        // 3. Helicase string (slice)
        if opts.bench_ram {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_STRING>::from_slice(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_string().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase string \t(slice)", dna_len);
            summary.helicase_str = m.throughput(size);
        }

        // 4. Helicase string (reader)
        if opts.bench_ram {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_STRING>::from_reader(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_string().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase string \t(reader)", dna_len);
            summary.helicase_str_reader = m.throughput(size);
        }

        // 5. Helicase string + PackedSeqVec (reused buffer)
        if opts.bench_ram || opts.pack_compare {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_STRING>::from_slice(slice).unwrap();
                let mut packed = PackedSeqVec::default();
                dna_len = 0;
                while parser.next().is_some() {
                    packed.clear();
                    packed.push_ascii(parser.get_dna_string());
                    dna_len += packed.as_slice().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase str+PackedSeqVec", dna_len);
            summary.helicase_str_packed = m.throughput(size);
        }

        // 6. Helicase packed (keep_non_actg) (slice)
        if opts.bench_ram || opts.pack_compare {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_PACKED_KEEP>::from_slice(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_packed().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase packed (keep)\t(slice)", dna_len);
            summary.helicase_packed_keep = m.throughput(size);
        }

        // 7. Helicase packed (keep_non_actg) (reader)
        if opts.bench_ram {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_PACKED_KEEP>::from_reader(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_packed().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase packed (keep)\t(reader)", dna_len);
            summary.helicase_packed_reader = m.throughput(size);
        }

        // 8. Helicase packed (skip_non_actg)
        if opts.bench_ram || opts.pack_compare {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_PACKED>::from_slice(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_packed().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase packed (skip)\t(slice)", dna_len);
            summary.helicase_packed = m.throughput(size);
        }

        // 9. Helicase columnar
        if opts.bench_ram || opts.pack_compare {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut parser = FastxParser::<DNA_COLUMNAR>::from_slice(slice).unwrap();
                dna_len = 0;
                while parser.next().is_some() {
                    dna_len += parser.get_dna_columnar().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Helicase columnar \t(slice)", dna_len);
            summary.helicase_col = m.throughput(size);
        }

        // 10. File & Mmap benchmarks on synthetic dataset
        if opts.bench_file || opts.bench_mmap {
            let file_ext = if ds.name.starts_with("FASTQ") {
                "fastq"
            } else {
                "fasta"
            };
            let temp_file_path = temp_dir.join(format!("ds_{idx}.{file_ext}"));
            if let Ok(mut f) = File::create(&temp_file_path) {
                let _ = f.write_all(slice);
            }

            if opts.bench_file {
                // Needletail file
                if !opts.no_baseline {
                    m.start();
                    let mut dna_len = 0;
                    for _ in 0..opts.repeat {
                        let mut reader = parse_fastx_file(&temp_file_path).unwrap();
                        dna_len = 0;
                        while let Some(r) = reader.next() {
                            dna_len += r.unwrap().seq().len();
                        }
                        m.tick();
                    }
                    show_tick(&m, "Needletail string \t(file)", dna_len);
                    summary.needletail_file = m.throughput(size);
                }

                // Helicase string file
                m.start();
                let mut dna_len = 0;
                for _ in 0..opts.repeat {
                    let mut parser = FastxParser::<DNA_STRING>::from_file(&temp_file_path).unwrap();
                    dna_len = 0;
                    while parser.next().is_some() {
                        dna_len += parser.get_dna_string().len();
                    }
                    m.tick();
                }
                show_tick(&m, "Helicase string \t(file)", dna_len);
                summary.helicase_str_file = m.throughput(size);

                // Helicase packed file
                m.start();
                let mut dna_len = 0;
                for _ in 0..opts.repeat {
                    let mut parser =
                        FastxParser::<DNA_PACKED_KEEP>::from_file(&temp_file_path).unwrap();
                    dna_len = 0;
                    while parser.next().is_some() {
                        dna_len += parser.get_dna_packed().len();
                    }
                    m.tick();
                }
                show_tick(&m, "Helicase packed \t(file)", dna_len);
                summary.helicase_packed_file = m.throughput(size);
            }

            if opts.bench_mmap {
                // Helicase string mmap
                m.start();
                let mut dna_len = 0;
                for _ in 0..opts.repeat {
                    let mut parser =
                        FastxParser::<DNA_STRING>::from_file_mmap(&temp_file_path).unwrap();
                    dna_len = 0;
                    while parser.next().is_some() {
                        dna_len += parser.get_dna_string().len();
                    }
                    m.tick();
                }
                show_tick(&m, "Helicase string \t(mmap)", dna_len);
                summary.helicase_str_mmap = m.throughput(size);

                // Helicase packed mmap
                m.start();
                let mut dna_len = 0;
                for _ in 0..opts.repeat {
                    let mut parser =
                        FastxParser::<DNA_PACKED_KEEP>::from_file_mmap(&temp_file_path).unwrap();
                    dna_len = 0;
                    while parser.next().is_some() {
                        dna_len += parser.get_dna_packed().len();
                    }
                    m.tick();
                }
                show_tick(&m, "Helicase packed \t(mmap)", dna_len);
                summary.helicase_packed_mmap = m.throughput(size);
            }

            let _ = std::fs::remove_file(&temp_file_path);
        }

        summaries.push(summary);
    }

    if opts.bench_file || opts.bench_mmap {
        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    // Output formatting
    match opts.format {
        OutputFormat::Text => {
            print_summary_table(&summaries, !opts.no_baseline);
            if opts.pack_compare {
                print_pack_comparison_table(&summaries, !opts.no_baseline);
            }
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_table(&summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Markdown => {
            print_markdown_table(&summaries, !opts.no_baseline);
            if opts.pack_compare {
                print_pack_comparison_markdown(&summaries, !opts.no_baseline);
            }
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_markdown(&summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Csv => {}
        OutputFormat::Json => {
            print_json_summary(&summaries);
        }
    }
}

pub fn run_custom_file_bench<M: Measurement>(path: &Path, opts: &BenchOptions) {
    let size = match std::fs::metadata(path) {
        Ok(m) => m.len(),
        Err(e) => {
            eprintln!("Error reading file {}: {}", path.display(), e);
            return;
        }
    };
    let size_mb = size as f64 / 1_000_000.0;
    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "custom_file".to_string());

    if opts.format == OutputFormat::Text {
        println!(
            "\n{} {} ({:.1} MB)",
            "==> Benchmarking Custom File:".bold().green(),
            name.bold(),
            size_mb
        );
    }

    let mut summary = DatasetSummary {
        name: name.clone(),
        size_mb,
        ..Default::default()
    };

    let is_csv = opts.format == OutputFormat::Csv;
    let is_human_text = opts.format == OutputFormat::Text;
    let show_ticks = is_human_text || is_csv;
    let mut m = M::new(&name);
    let show_tick = |m: &M, label: &str, dna_len: usize| {
        if show_ticks {
            m.show(
                label,
                size,
                if opts.show_val { Some(dna_len) } else { None },
                is_csv,
            );
        }
    };

    // File input (streaming)
    if opts.bench_file {
        if !opts.no_baseline {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut reader = parse_fastx_file(path).expect("Failed to parse file");
                dna_len = 0;
                while let Some(r) = reader.next() {
                    dna_len += r.expect("Invalid record").seq().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Needletail string \t(file)", dna_len);
            summary.needletail_file = m.throughput(size);
        }

        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser = FastxParser::<DNA_STRING>::from_file(path).expect("Cannot open file");
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_string().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase string \t(file)", dna_len);
        summary.helicase_str_file = m.throughput(size);

        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser =
                FastxParser::<DNA_PACKED_KEEP>::from_file(path).expect("Cannot open file");
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_packed().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase packed \t(file)", dna_len);
        summary.helicase_packed_file = m.throughput(size);
    }

    // Mmap input
    if opts.bench_mmap {
        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser =
                FastxParser::<DNA_STRING>::from_file_mmap(path).expect("Cannot open file");
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_string().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase string \t(mmap)", dna_len);
        summary.helicase_str_mmap = m.throughput(size);

        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser =
                FastxParser::<DNA_PACKED_KEEP>::from_file_mmap(path).expect("Cannot open file");
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_packed().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase packed \t(mmap)", dna_len);
        summary.helicase_packed_mmap = m.throughput(size);
    }

    // RAM slice input (reads whole file into RAM)
    if opts.bench_ram
        && let Ok(data) = std::fs::read(path)
    {
        let slice = data.as_slice();

        if !opts.no_baseline {
            m.start();
            let mut dna_len = 0;
            for _ in 0..opts.repeat {
                let mut reader = parse_fastx_reader(slice).expect("Failed to parse slice");
                dna_len = 0;
                while let Some(r) = reader.next() {
                    dna_len += r.expect("Invalid record").seq().len();
                }
                m.tick();
            }
            std::hint::black_box(dna_len);
            show_tick(&m, "Needletail string \t(reader)", dna_len);
            summary.needletail_str = m.throughput(size);
        }

        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser = FastxParser::<DNA_STRING>::from_slice(slice).unwrap();
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_string().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase string \t(slice)", dna_len);
        summary.helicase_str = m.throughput(size);

        m.start();
        let mut dna_len = 0;
        for _ in 0..opts.repeat {
            let mut parser = FastxParser::<DNA_PACKED_KEEP>::from_slice(slice).unwrap();
            dna_len = 0;
            while parser.next().is_some() {
                dna_len += parser.get_dna_packed().len();
            }
            m.tick();
        }
        std::hint::black_box(dna_len);
        show_tick(&m, "Helicase packed \t(slice)", dna_len);
        summary.helicase_packed_keep = m.throughput(size);
    }

    let summaries = vec![summary];
    match opts.format {
        OutputFormat::Text => {
            print_summary_table(&summaries, !opts.no_baseline);
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_table(&summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Markdown => {
            print_markdown_table(&summaries, !opts.no_baseline);
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_markdown(&summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Csv => {}
        OutputFormat::Json => {
            print_json_summary(&summaries);
        }
    }
}

// ============================================================================
// Scoreboard Printing Routines
// ============================================================================

pub fn print_summary_table(summaries: &[DatasetSummary], has_baseline: bool) {
    println!("\n{}", "=".repeat(126).dimmed());
    println!(
        "{}",
        "SYNTHETIC BENCHMARK SCOREBOARD (Throughput in GB/s)".bold()
    );
    println!("{}", "-".repeat(126).dimmed());

    if has_baseline {
        println!(
            "{:<22} {:>7} | {:>10} {:>10} {:>10} {:>8} | {:>10} {:>10} {:>10} {:>9}",
            "Dataset",
            "Size",
            "Needletail",
            "Heli(slc)",
            "Heli(rdr)",
            "vs Base",
            "Needle+pk",
            "H-pk(slc)",
            "H-pk(rdr)",
            "vs N+pk"
        );
    } else {
        println!(
            "{:<24} {:>8} {:>14} {:>14} {:>14} {:>14} {:>14}",
            "Dataset",
            "Size",
            "Heli str(slc)",
            "Heli str(rdr)",
            "Heli pk(slc)",
            "Heli pk(rdr)",
            "Heli col(slc)"
        );
    }
    println!("{}", "-".repeat(126).dimmed());

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        let h_str = s
            .helicase_str
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_default();
        let h_str_rdr = s
            .helicase_str_reader
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_default();
        let h_pack_val = s.helicase_packed_keep.or(s.helicase_packed);
        let h_pack = h_pack_val
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_default();
        let h_pack_rdr = s
            .helicase_packed_reader
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_default();

        if has_baseline {
            let n_str = s
                .needletail_str
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_default();
            let diff_str = match (s.needletail_str, s.helicase_str) {
                (Some((n, _)), Some((h, _))) => {
                    let diff_pct = ((h - n) / n) * 100.0;
                    if diff_pct >= 0.0 {
                        format!("+{:5.1}%", diff_pct).green().bold().to_string()
                    } else {
                        format!("{:5.1}%", diff_pct).red().bold().to_string()
                    }
                }
                _ => String::new(),
            };

            let np_str = s
                .needletail_packed
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_default();
            let diff_pack_str = match (s.needletail_packed, h_pack_val) {
                (Some((np, _)), Some((hp, _))) => {
                    let diff_pct = ((hp - np) / np) * 100.0;
                    if diff_pct >= 0.0 {
                        format!("+{:5.1}%", diff_pct).green().bold().to_string()
                    } else {
                        format!("{:5.1}%", diff_pct).red().bold().to_string()
                    }
                }
                _ => String::new(),
            };

            println!(
                "{:<22} {:>7} | {:>10} {:>10} {:>10} {:>8} | {:>10} {:>10} {:>10} {:>9}",
                s.name,
                size_str,
                n_str,
                h_str.bright_green(),
                h_str_rdr.cyan(),
                diff_str,
                np_str,
                h_pack.bright_green(),
                h_pack_rdr.cyan(),
                diff_pack_str
            );
        } else {
            let h_col = s
                .helicase_col
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_default();
            println!(
                "{:<24} {:>8} {:>14} {:>14} {:>14} {:>14} {:>14}",
                s.name,
                size_str,
                h_str.bright_green(),
                h_str_rdr.cyan(),
                h_pack,
                h_pack_rdr.cyan(),
                h_col
            );
        }
    }
    println!("{}\n", "=".repeat(126).dimmed());
}

pub fn print_markdown_table(summaries: &[DatasetSummary], has_baseline: bool) {
    println!(
        "\n| Benchmark | Contig / Read Size | Helicase (str) | Needletail (str) | Delta vs Needletail | Helicase (pack) | Helicase (col) |"
    );
    println!("|---|---|---|---|---|---|---|");

    for s in summaries {
        let h_str = s
            .helicase_str
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());
        let n_str = if has_baseline {
            s.needletail_str
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_else(|| "N/A".to_string())
        } else {
            "-".to_string()
        };
        let delta = match (s.needletail_str, s.helicase_str) {
            (Some((n, _)), Some((h, _))) if has_baseline => {
                let diff_pct = ((h - n) / n) * 100.0;
                format!("**+{:4.1}%**", diff_pct)
            }
            _ => "-".to_string(),
        };
        let h_pack_val = s.helicase_packed_keep.or(s.helicase_packed);
        let h_pack = h_pack_val
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());
        let h_col = s
            .helicase_col
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());

        println!(
            "| **{}** | {:.1} MB | {} | {} | {} | {} | {} |",
            s.name, s.size_mb, h_str, n_str, delta, h_pack, h_col
        );
    }
}

pub fn print_pack_comparison_table(summaries: &[DatasetSummary], has_baseline: bool) {
    println!("{}", "=".repeat(118).dimmed());
    println!(
        "{}",
        "BITPACKING BREAKDOWN SCOREBOARD (Throughput in GB/s & Delta vs ASCII conversion)"
            .bold()
            .yellow()
    );
    println!("{}", "-".repeat(118).dimmed());
    if has_baseline {
        println!(
            "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12} | {:>14} {:>14}",
            "Dataset",
            "Size",
            "Needle+pkSeq",
            "Heli+pkSeq",
            "Heli(pack)",
            "Heli(col)",
            "vs Needle+pk",
            "vs Heli+pkSeq"
        );
    } else {
        println!(
            "{:<22} {:>7} | {:>12} {:>12} {:>12} | {:>14}",
            "Dataset", "Size", "Heli+pkSeq", "Heli(pack)", "Heli(col)", "vs Heli+pkSeq"
        );
    }
    println!("{}", "-".repeat(118).dimmed());

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        let np_str = s
            .needletail_packed
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let hp_str = s
            .helicase_str_packed
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let h_pack_val = s.helicase_packed_keep.or(s.helicase_packed);
        let pack_str = h_pack_val
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let col_str = s
            .helicase_col
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());

        let vs_needle = match (s.needletail_packed, h_pack_val) {
            (Some((n, _)), Some((h, _))) => {
                let diff = ((h - n) / n) * 100.0;
                if diff >= 0.0 {
                    format!("+{:5.1}%", diff).green().to_string()
                } else {
                    format!("{:6.1}%", diff).red().to_string()
                }
            }
            _ => "-".to_string(),
        };

        let vs_heli_str_pack = match (s.helicase_str_packed, h_pack_val) {
            (Some((h_sp, _)), Some((h, _))) => {
                let diff = ((h - h_sp) / h_sp) * 100.0;
                if diff >= 0.0 {
                    format!("+{:5.1}%", diff).cyan().to_string()
                } else {
                    format!("{:6.1}%", diff).red().to_string()
                }
            }
            _ => "-".to_string(),
        };

        if has_baseline {
            println!(
                "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12} | {:>14} {:>14}",
                s.name,
                size_str,
                np_str,
                hp_str,
                pack_str.bright_green(),
                col_str,
                vs_needle,
                vs_heli_str_pack
            );
        } else {
            println!(
                "{:<22} {:>7} | {:>12} {:>12} {:>12} | {:>14}",
                s.name,
                size_str,
                hp_str,
                pack_str.bright_green(),
                col_str,
                vs_heli_str_pack
            );
        }
    }
    println!("{}\n", "=".repeat(118).dimmed());
}

pub fn print_pack_comparison_markdown(summaries: &[DatasetSummary], has_baseline: bool) {
    println!("\n### 2-Bit Bitpacking Comparison");
    if has_baseline {
        println!(
            "| Dataset | Size | Needletail + pkSeq | Heli str + pkSeq | Heli native pack | Heli native col | Delta vs Needle+pk | Delta vs Heli+pk |"
        );
        println!("|---|---|---|---|---|---|---|---|");
    } else {
        println!(
            "| Dataset | Size | Heli str + pkSeq | Heli native pack | Heli native col | Delta vs Heli+pk |"
        );
        println!("|---|---|---|---|---|---|");
    }

    for s in summaries {
        let np = s
            .needletail_packed
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let hp = s
            .helicase_str_packed
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let h_pack_val = s.helicase_packed_keep.or(s.helicase_packed);
        let pack = h_pack_val
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "-".to_string());
        let col = s
            .helicase_col
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let vs_needle = match (s.needletail_packed, h_pack_val) {
            (Some((n, _)), Some((h, _))) => {
                let diff = ((h - n) / n) * 100.0;
                if diff >= 0.0 {
                    format!("+{:4.1}%", diff)
                } else {
                    format!("{:5.1}%", diff)
                }
            }
            _ => "-".to_string(),
        };
        let vs_heli = match (s.helicase_str_packed, h_pack_val) {
            (Some((sp, _)), Some((h, _))) => {
                let diff = ((h - sp) / sp) * 100.0;
                if diff >= 0.0 {
                    format!("+{:4.1}%", diff)
                } else {
                    format!("{:5.1}%", diff)
                }
            }
            _ => "-".to_string(),
        };

        if has_baseline {
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} | {} | {} |",
                s.name, s.size_mb, np, hp, pack, col, vs_needle, vs_heli
            );
        } else {
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} |",
                s.name, s.size_mb, hp, pack, col, vs_heli
            );
        }
    }
}

pub fn print_io_comparison_table(summaries: &[DatasetSummary], has_baseline: bool) {
    println!("{}", "=".repeat(118).dimmed());
    println!(
        "{}",
        "I/O & BUFFERING SCOREBOARD (RAM vs 128KB FileInput vs MmapInput)"
            .bold()
            .blue()
    );
    println!("{}", "-".repeat(118).dimmed());
    if has_baseline {
        println!(
            "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12} | {:>12} {:>12}",
            "Dataset",
            "Size",
            "RAM (slice)",
            "RAM (reader)",
            "FileInput",
            "MmapInput",
            "Needle(file)",
            "vs Needle"
        );
    } else {
        println!(
            "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12}",
            "Dataset", "Size", "RAM (slice)", "RAM (reader)", "FileInput", "MmapInput"
        );
    }
    println!("{}", "-".repeat(118).dimmed());

    for s in summaries {
        let size_str = format!("{:.1} MB", s.size_mb);
        let ram_slc = s
            .helicase_str
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let ram_rdr = s
            .helicase_str_reader
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let file_str = s
            .helicase_str_file
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let mmap_str = s
            .helicase_str_mmap
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());

        if has_baseline {
            let needle_f = s
                .needletail_file
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_else(|| "-".to_string());
            let vs_needle = match (s.needletail_file, s.helicase_str_file) {
                (Some((n, _)), Some((h, _))) => {
                    let diff = ((h - n) / n) * 100.0;
                    format!("+{:5.1}%", diff).green().to_string()
                }
                _ => "-".to_string(),
            };
            println!(
                "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12} | {:>12} {:>12}",
                s.name,
                size_str,
                ram_slc,
                ram_rdr,
                file_str.bright_green(),
                mmap_str.cyan(),
                needle_f,
                vs_needle
            );
        } else {
            println!(
                "{:<22} {:>7} | {:>12} {:>12} {:>12} {:>12}",
                s.name,
                size_str,
                ram_slc,
                ram_rdr,
                file_str.bright_green(),
                mmap_str.cyan(),
            );
        }
    }
    println!("{}\n", "=".repeat(118).dimmed());
}

pub fn print_io_comparison_markdown(summaries: &[DatasetSummary], has_baseline: bool) {
    println!("\n### I/O & Streaming Buffering Comparison");
    if has_baseline {
        println!(
            "| Dataset | Size | RAM Slice | RAM Reader | FileInput (128KB) | MmapInput | Needletail File | vs Needletail |"
        );
        println!("|---|---|---|---|---|---|---|---|");
    } else {
        println!("| Dataset | Size | RAM Slice | RAM Reader | FileInput (128KB) | MmapInput |");
        println!("|---|---|---|---|---|---|");
    }

    for s in summaries {
        let ram_slc = s
            .helicase_str
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let ram_rdr = s
            .helicase_str_reader
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());
        let file_str = s
            .helicase_str_file
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "-".to_string());
        let mmap_str = s
            .helicase_str_mmap
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "-".to_string());

        if has_baseline {
            let needle_f = s
                .needletail_file
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_else(|| "-".to_string());
            let vs_needle = match (s.needletail_file, s.helicase_str_file) {
                (Some((n, _)), Some((h, _))) => {
                    format!("+{:4.1}%", ((h - n) / n) * 100.0)
                }
                _ => "-".to_string(),
            };
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} | {} | {} |",
                s.name, s.size_mb, ram_slc, ram_rdr, file_str, mmap_str, needle_f, vs_needle
            );
        } else {
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} |",
                s.name, s.size_mb, ram_slc, ram_rdr, file_str, mmap_str
            );
        }
    }
}

pub fn print_json_summary(summaries: &[DatasetSummary]) {
    println!("[");
    for (i, s) in summaries.iter().enumerate() {
        let comma = if i + 1 < summaries.len() { "," } else { "" };
        println!(
            "  {{\n    \"name\": \"{}\",\n    \"size_mb\": {:.2},\n    \"helicase_str\": {},\n    \"needletail_str\": {},\n    \"helicase_pack\": {},\n    \"helicase_col\": {}\n  }}{}",
            s.name,
            s.size_mb,
            s.helicase_str
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.needletail_str
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.helicase_packed_keep
                .or(s.helicase_packed)
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.helicase_col
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            comma
        );
    }
    println!("]");
}
