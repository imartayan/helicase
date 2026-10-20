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
    pub warmup: u64,
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
            warmup: 1,
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

        if opts.bench_ram || opts.pack_compare {
            bench_ram_dataset::<M>(ds.name, slice, size, opts, &mut summary);
        }

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
            bench_file_and_mmap::<M>(ds.name, &temp_file_path, size, opts, &mut summary);
            let _ = std::fs::remove_file(&temp_file_path);
        }

        summaries.push(summary);
    }

    if opts.bench_file || opts.bench_mmap {
        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    print_results(&summaries, opts);
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

    let is_compressed = crate::cli::is_compressed_file(path);

    if !is_compressed && (opts.bench_ram || opts.pack_compare) {
        match std::fs::read(path) {
            Ok(data) => {
                bench_ram_dataset::<M>(&name, &data, size, opts, &mut summary);
            }
            Err(e) => {
                eprintln!("Error reading file {} into RAM: {}", path.display(), e);
            }
        }
    }

    if opts.bench_file || opts.bench_mmap {
        bench_file_and_mmap::<M>(&name, path, size, opts, &mut summary);
    }

    let summaries = vec![summary];
    print_results(&summaries, opts);
}

macro_rules! bench_run {
    ($m:ident, $opts:expr, $code:block) => {{
        for _ in 0..$opts.warmup {
            let len: usize = $code;
            std::hint::black_box(len);
        }
        $m.start();
        let mut last_len: usize = 0;
        for _ in 0..$opts.repeat {
            last_len = $code;
            $m.tick();
        }
        std::hint::black_box(last_len);
        last_len
    }};
}

fn bench_ram_dataset<M: Measurement>(
    name: &str,
    slice: &[u8],
    size: u64,
    opts: &BenchOptions,
    summary: &mut DatasetSummary,
) {
    if slice.is_empty() {
        return;
    }

    let is_csv = opts.format == OutputFormat::Csv;
    let is_human_text = opts.format == OutputFormat::Text;
    let show_ticks = is_human_text || is_csv;
    let mut m = M::new(name);
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
        let dna_len = bench_run!(m, opts, {
            let mut reader = parse_fastx_reader(slice).expect("Failed to parse slice");
            let mut len = 0;
            while let Some(r) = reader.next() {
                let record = r.expect("Invalid record");
                let seq = record.seq();
                len += seq.len();
            }
            len
        });
        show_tick(&m, "Needletail string \t(reader)", dna_len);
        summary.needletail_str = m.throughput(size);
    }

    // 2. Needletail string + PackedSeqVec (reused buffer)
    if !opts.no_baseline && (opts.bench_ram || opts.pack_compare) {
        let mut packed = PackedSeqVec::default();
        let dna_len = bench_run!(m, opts, {
            let mut reader = parse_fastx_reader(slice).expect("Failed to parse slice");
            let mut len = 0;
            while let Some(r) = reader.next() {
                let record = r.expect("Invalid record");
                packed.clear();
                packed.push_ascii(&record.seq());
                len += packed.as_slice().len();
            }
            len
        });
        show_tick(&m, "Needletail+PackedSeqVec", dna_len);
        summary.needletail_packed = m.throughput(size);
    }

    // 3. Helicase string (reader)
    if opts.bench_ram {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_STRING>::from_reader(slice).expect("Failed to parse reader");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_string().len();
            }
            len
        });
        show_tick(&m, "Helicase string \t(reader)", dna_len);
        summary.helicase_str_reader = m.throughput(size);
    }

    // 4. Helicase string (slice)
    if opts.bench_ram {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_STRING>::from_slice(slice).expect("Failed to parse slice");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_string().len();
            }
            len
        });
        show_tick(&m, "Helicase string \t(slice)", dna_len);
        summary.helicase_str = m.throughput(size);
    }

    // 5. Helicase string + PackedSeqVec (reused buffer)
    if opts.bench_ram || opts.pack_compare {
        let mut packed = PackedSeqVec::default();
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_STRING>::from_slice(slice).expect("Failed to parse slice");
            let mut len = 0;
            while parser.next().is_some() {
                packed.clear();
                packed.push_ascii(parser.get_dna_string());
                len += packed.as_slice().len();
            }
            len
        });
        show_tick(&m, "Helicase str+PackedSeqVec", dna_len);
        summary.helicase_str_packed = m.throughput(size);
    }

    // 6. Helicase packed (keep_non_actg) (reader)
    if opts.bench_ram {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_PACKED_KEEP>::from_reader(slice).expect("Failed to parse reader");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_packed().len();
            }
            len
        });
        show_tick(&m, "Helicase packed (keep)\t(reader)", dna_len);
        summary.helicase_packed_reader = m.throughput(size);
    }

    // 7. Helicase packed (keep_non_actg) (slice)
    if opts.bench_ram || opts.pack_compare {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_PACKED_KEEP>::from_slice(slice).expect("Failed to parse slice");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_packed().len();
            }
            len
        });
        show_tick(&m, "Helicase packed (keep)\t(slice)", dna_len);
        summary.helicase_packed_keep = m.throughput(size);
    }

    // 8. Helicase packed (skip_non_actg)
    if opts.bench_ram || opts.pack_compare {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_PACKED>::from_slice(slice).expect("Failed to parse slice");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_packed().len();
            }
            len
        });
        show_tick(&m, "Helicase packed (skip)\t(slice)", dna_len);
        summary.helicase_packed = m.throughput(size);
    }

    // 9. Helicase columnar
    if opts.bench_ram || opts.pack_compare {
        let dna_len = bench_run!(m, opts, {
            let mut parser =
                FastxParser::<DNA_COLUMNAR>::from_slice(slice).expect("Failed to parse slice");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_columnar().len();
            }
            len
        });
        show_tick(&m, "Helicase columnar \t(slice)", dna_len);
        summary.helicase_col = m.throughput(size);
    }
}

pub fn open_needletail_reader(
    path: &Path,
) -> Result<Box<dyn needletail::FastxReader>, Box<dyn std::error::Error + Send + Sync>> {
    if let Ok(reader) = parse_fastx_file(path) {
        return Ok(reader);
    }

    // Needletail only checks standard Zstd magic [0x28, 0xB5] and fails to recognize
    // multi-frame Zstd skippable frames (0x184D2A50..0x184D2A5F) emitted by pzstd.
    // Wrap with zstd decoder if the file is Zstd/Pzstd.
    let format = helicase::input::CompressionFormat::from_path(path)
        .ok()
        .flatten();
    if matches!(
        format,
        Some(helicase::input::CompressionFormat::Zstd | helicase::input::CompressionFormat::Pzstd)
    ) {
        let file = File::open(path)?;
        let decoder = zstd::stream::read::Decoder::new(file)?;
        return Ok(parse_fastx_reader(decoder)?);
    }

    Ok(parse_fastx_file(path)?)
}

fn bench_file_and_mmap<M: Measurement>(
    name: &str,
    path: &Path,
    size: u64,
    opts: &BenchOptions,
    summary: &mut DatasetSummary,
) {
    let is_csv = opts.format == OutputFormat::Csv;
    let is_human_text = opts.format == OutputFormat::Text;
    let show_ticks = is_human_text || is_csv;
    let mut m = M::new(name);
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

    if opts.bench_file {
        // Needletail file
        if !opts.no_baseline {
            match open_needletail_reader(path) {
                Ok(_) => {
                    let dna_len = bench_run!(m, opts, {
                        let mut reader = open_needletail_reader(path)
                            .expect("Failed to parse file with Needletail");
                        let mut len = 0;
                        while let Some(r) = reader.next() {
                            len += r.expect("Invalid record").seq().len();
                        }
                        len
                    });
                    show_tick(&m, "Needletail string \t(file)", dna_len);
                    summary.needletail_file = m.throughput(size);
                }
                Err(e) => {
                    eprintln!(
                        "Warning: Needletail does not support {}: {e}; skipping baseline.",
                        path.display()
                    );
                }
            }
        }

        // Helicase string file
        let dna_len = bench_run!(m, opts, {
            let mut parser = FastxParser::<DNA_STRING>::from_file(path)
                .expect("Failed to open file for Helicase string");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_string().len();
            }
            len
        });
        show_tick(&m, "Helicase string \t(file)", dna_len);
        summary.helicase_str_file = m.throughput(size);

        // Helicase packed file
        let dna_len = bench_run!(m, opts, {
            let mut parser = FastxParser::<DNA_PACKED_KEEP>::from_file(path)
                .expect("Failed to open file for Helicase packed");
            let mut len = 0;
            while parser.next().is_some() {
                len += parser.get_dna_packed().len();
            }
            len
        });
        show_tick(&m, "Helicase packed \t(file)", dna_len);
        summary.helicase_packed_file = m.throughput(size);
    }

    let is_compressed = crate::cli::is_compressed_file(path);

    if opts.bench_mmap {
        if is_compressed {
            eprintln!(
                "{}",
                "Warning: Memory-mapping (-M) is not supported for compressed files in single-thread mode; skipping."
                    .yellow()
            );
        } else {
            // Helicase string mmap
            let dna_len = bench_run!(m, opts, {
                let mut parser = FastxParser::<DNA_STRING>::from_file_mmap(path)
                    .expect("Failed to mmap file for Helicase string");
                let mut len = 0;
                while parser.next().is_some() {
                    len += parser.get_dna_string().len();
                }
                len
            });
            show_tick(&m, "Helicase string \t(mmap)", dna_len);
            summary.helicase_str_mmap = m.throughput(size);

            // Helicase packed mmap
            let dna_len = bench_run!(m, opts, {
                let mut parser = FastxParser::<DNA_PACKED_KEEP>::from_file_mmap(path)
                    .expect("Failed to mmap file for Helicase packed");
                let mut len = 0;
                while parser.next().is_some() {
                    len += parser.get_dna_packed().len();
                }
                len
            });
            show_tick(&m, "Helicase packed \t(mmap)", dna_len);
            summary.helicase_packed_mmap = m.throughput(size);
        }
    }
}

pub fn print_results(summaries: &[DatasetSummary], opts: &BenchOptions) {
    let has_ram = summaries.iter().any(|s| {
        s.helicase_str.is_some()
            || s.helicase_str_reader.is_some()
            || s.needletail_str.is_some()
            || s.helicase_packed.is_some()
            || s.helicase_col.is_some()
    });

    match opts.format {
        OutputFormat::Text => {
            if has_ram {
                print_summary_table(summaries, !opts.no_baseline);
            }
            if opts.pack_compare {
                print_pack_comparison_table(summaries, !opts.no_baseline);
            }
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_table(summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Markdown => {
            if has_ram {
                print_markdown_table(summaries, !opts.no_baseline);
            }
            if opts.pack_compare {
                print_pack_comparison_markdown(summaries, !opts.no_baseline);
            }
            if opts.bench_file || opts.bench_mmap {
                print_io_comparison_markdown(summaries, !opts.no_baseline);
            }
        }
        OutputFormat::Csv => {}
        OutputFormat::Json => {
            print_json_summary(summaries);
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
            "{:<22} {:>7} | {:>10} {:>10} {:>8} {:>10} | {:>10} {:>10} {:>9} {:>10}",
            "Dataset",
            "Size",
            "Needletail",
            "Heli(rdr)",
            "vs Base",
            "Heli(slc)",
            "Needle+pk",
            "H-pk(rdr)",
            "vs N+pk",
            "H-pk(slc)"
        );
    } else {
        println!(
            "{:<24} {:>8} {:>14} {:>14} {:>14} {:>14} {:>14}",
            "Dataset",
            "Size",
            "Heli str(rdr)",
            "Heli str(slc)",
            "Heli pk(rdr)",
            "Heli pk(slc)",
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
            let h_str_comp = s.helicase_str_reader.or(s.helicase_str);
            let diff_str = match (s.needletail_str, h_str_comp) {
                (Some((n, _)), Some((h, _))) => {
                    let diff_pct = ((h - n) / n) * 100.0;
                    if diff_pct >= 0.0 {
                        format!("+{:.1}%", diff_pct).green().bold().to_string()
                    } else {
                        format!("{:.1}%", diff_pct).red().bold().to_string()
                    }
                }
                _ => String::new(),
            };

            let np_str = s
                .needletail_packed
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_default();
            let h_pack_comp = s.helicase_packed_reader.or(h_pack_val);
            let diff_pack_str = match (s.needletail_packed, h_pack_comp) {
                (Some((np, _)), Some((hp, _))) => {
                    let diff_pct = ((hp - np) / np) * 100.0;
                    if diff_pct >= 0.0 {
                        format!("+{:.1}%", diff_pct).green().bold().to_string()
                    } else {
                        format!("{:.1}%", diff_pct).red().bold().to_string()
                    }
                }
                _ => String::new(),
            };

            println!(
                "{:<22} {:>7} | {:>10} {:>10} {:>8} {:>10} | {:>10} {:>10} {:>9} {:>10}",
                s.name,
                size_str,
                n_str,
                h_str_rdr.bright_green(),
                diff_str,
                h_str.cyan(),
                np_str,
                h_pack_rdr.bright_green(),
                diff_pack_str,
                h_pack.cyan()
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
                h_str_rdr.bright_green(),
                h_str.cyan(),
                h_pack_rdr.bright_green(),
                h_pack.cyan(),
                h_col
            );
        }
    }
    println!("{}\n", "=".repeat(126).dimmed());
}

pub fn print_markdown_table(summaries: &[DatasetSummary], has_baseline: bool) {
    if has_baseline {
        println!(
            "\n| Benchmark | Contig / Read Size | Helicase (rdr) | Needletail (str) | Delta vs Needletail | Helicase (slc) | Helicase (pack) | Helicase (col) |"
        );
        println!("|---|---|---|---|---|---|---|---|");
    } else {
        println!(
            "\n| Benchmark | Contig / Read Size | Helicase (rdr) | Helicase (slc) | Helicase (pack) | Helicase (col) |"
        );
        println!("|---|---|---|---|---|---|");
    }

    for s in summaries {
        let h_rdr = s
            .helicase_str_reader
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());
        let h_slc = s
            .helicase_str
            .map(|(m, _)| format!("{:5.2} GB/s", m))
            .unwrap_or_else(|| "N/A".to_string());
        let n_str = if has_baseline {
            s.needletail_str
                .map(|(m, _)| format!("{:5.2} GB/s", m))
                .unwrap_or_else(|| "N/A".to_string())
        } else {
            "-".to_string()
        };
        let delta = match (s.needletail_str, s.helicase_str_reader.or(s.helicase_str)) {
            (Some((n, _)), Some((h, _))) if has_baseline => {
                let diff_pct = ((h - n) / n) * 100.0;
                if diff_pct >= 0.0 {
                    format!("**+{:4.1}%**", diff_pct)
                } else {
                    format!("**{:5.1}%**", diff_pct)
                }
            }
            _ => "-".to_string(),
        };
        let h_pack_val = s
            .helicase_packed_reader
            .or(s.helicase_packed_keep)
            .or(s.helicase_packed);
        let h_pack = h_pack_val
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());
        let h_col = s
            .helicase_col
            .map(|(m, _)| format!("**{:5.2} GB/s**", m))
            .unwrap_or_else(|| "N/A".to_string());

        if has_baseline {
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} | {} | {} |",
                s.name, s.size_mb, h_rdr, n_str, delta, h_slc, h_pack, h_col
            );
        } else {
            println!(
                "| **{}** | {:.1} MB | {} | {} | {} | {} |",
                s.name, s.size_mb, h_rdr, h_slc, h_pack, h_col
            );
        }
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
                    if diff >= 0.0 {
                        format!("+{:5.1}%", diff).green().to_string()
                    } else {
                        format!("{:6.1}%", diff).red().to_string()
                    }
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
                    let diff = ((h - n) / n) * 100.0;
                    if diff >= 0.0 {
                        format!("+{:4.1}%", diff)
                    } else {
                        format!("{:5.1}%", diff)
                    }
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
            "  {{\n    \"name\": \"{}\",\n    \"size_mb\": {:.2},\n    \"helicase_str_rdr\": {},\n    \"helicase_str_slc\": {},\n    \"needletail_str\": {},\n    \"helicase_pack_rdr\": {},\n    \"helicase_pack_slc\": {},\n    \"helicase_col\": {}\n  }}{}",
            s.name,
            s.size_mb,
            s.helicase_str_reader
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.helicase_str
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.needletail_str
                .map(|(m, _)| format!("{:.2}", m))
                .unwrap_or_else(|| "null".to_string()),
            s.helicase_packed_reader
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
