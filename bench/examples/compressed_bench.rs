use fastx_gen::FastxGen;
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use needletail::parse_fastx_file;
use niffler::compression::{Format as NifflerFormat, Level as NifflerLevel};
use paraseq::{ReaderBuilder, Record};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use std::time::Instant;

const CONFIG_STR: Config = ParserOptions::default().dna_string().config();
const CONFIG_PACKED: Config = ParserOptions::default().dna_packed().config();
const CONFIG_COL: Config = ParserOptions::default().dna_columnar().config();

struct BenchmarkResult {
    name: &'static str,
    raw_size_mb: f64,
    compressed_size_mb: f64,
    comp_ratio: f64,
    decomp_only_gb_s: f64,
    needletail_gb_s: f64,
    paraseq_gb_s: f64,
    helicase_str_gb_s: f64,
    helicase_pack_gb_s: f64,
    helicase_col_gb_s: f64,
}

fn measure_decomp_only(path: &Path, raw_len: usize) -> f64 {
    let mut min_time = f64::MAX;
    let mut buf = vec![0u8; 128 * 1024];

    for _ in 0..4 {
        let file = File::open(path).unwrap();
        let (mut reader, _) = niffler::get_reader(Box::new(file)).unwrap();
        let start = Instant::now();
        let mut total = 0;
        loop {
            let n = reader.read(&mut buf).unwrap();
            if n == 0 {
                break;
            }
            total += n;
        }
        let elapsed = start.elapsed().as_secs_f64();
        assert_eq!(total, raw_len);
        if elapsed < min_time {
            min_time = elapsed;
        }
    }
    (raw_len as f64 / 1_000_000_000.0) / min_time
}

fn measure_needletail(path: &Path, raw_len: usize) -> f64 {
    let mut min_time = f64::MAX;
    for _ in 0..4 {
        let start = Instant::now();
        let mut reader = parse_fastx_file(path).unwrap();
        let mut count = 0;
        let mut bases = 0;
        while let Some(r) = reader.next() {
            let rec = r.unwrap();
            count += 1;
            bases += rec.seq().len();
        }
        let elapsed = start.elapsed().as_secs_f64();
        assert!(count > 0 && bases > 0);
        if elapsed < min_time {
            min_time = elapsed;
        }
    }
    (raw_len as f64 / 1_000_000_000.0) / min_time
}

fn measure_paraseq(path: &Path, raw_len: usize) -> f64 {
    let mut min_time = f64::MAX;
    for _ in 0..4 {
        let start = Instant::now();
        let mut reader = ReaderBuilder::path(path).batch_size(1).build().unwrap();
        let mut record_set = reader.new_record_set();
        let mut bases = 0;
        while record_set.fill(&mut reader).unwrap() {
            for r in record_set.iter() {
                let rec = r.unwrap();
                bases += rec.seq().len();
            }
        }
        let elapsed = start.elapsed().as_secs_f64();
        assert!(bases > 0);
        if elapsed < min_time {
            min_time = elapsed;
        }
    }
    (raw_len as f64 / 1_000_000_000.0) / min_time
}

fn measure_helicase<const CONFIG: Config>(path: &Path, raw_len: usize, mode: &str) -> f64 {
    let mut min_time = f64::MAX;
    for _ in 0..4 {
        let start = Instant::now();
        let mut parser = FastxParser::<CONFIG>::from_file(path).unwrap();
        let mut count = 0;
        let mut bases = 0;
        while parser.next().is_some() {
            count += 1;
            match mode {
                "str" => bases += parser.get_dna_string().len(),
                "pack" => bases += parser.get_dna_packed().len(),
                "col" => bases += parser.get_dna_columnar().len(),
                _ => unreachable!(),
            }
        }
        let elapsed = start.elapsed().as_secs_f64();
        assert!(count > 0 && bases > 0);
        if elapsed < min_time {
            min_time = elapsed;
        }
    }
    (raw_len as f64 / 1_000_000_000.0) / min_time
}

fn run_compressed_suite(
    dataset_name: &'static str,
    raw_bytes: &[u8],
    format: NifflerFormat,
    ext: &str,
    temp_dir: &Path,
) -> BenchmarkResult {
    let file_path = temp_dir.join(format!("{dataset_name}.{ext}"));
    {
        let file = File::create(&file_path).unwrap();
        let level = match format {
            NifflerFormat::Gzip => NifflerLevel::Five,
            NifflerFormat::Zstd => NifflerLevel::Three,
            _ => unreachable!(),
        };
        let mut encoder = niffler::get_writer(Box::new(file), format, level).unwrap();
        encoder.write_all(raw_bytes).unwrap();
        encoder.flush().unwrap();
    }

    let compressed_len = std::fs::metadata(&file_path).unwrap().len() as usize;
    let raw_len = raw_bytes.len();
    let raw_mb = raw_len as f64 / (1024.0 * 1024.0);
    let comp_mb = compressed_len as f64 / (1024.0 * 1024.0);
    let ratio = raw_mb / comp_mb;

    let decomp_only = measure_decomp_only(&file_path, raw_len);
    let needletail_spd = measure_needletail(&file_path, raw_len);
    let paraseq_spd = measure_paraseq(&file_path, raw_len);
    let heli_str = measure_helicase::<CONFIG_STR>(&file_path, raw_len, "str");
    let heli_pack = measure_helicase::<CONFIG_PACKED>(&file_path, raw_len, "pack");
    let heli_col = measure_helicase::<CONFIG_COL>(&file_path, raw_len, "col");

    let _ = std::fs::remove_file(&file_path);

    BenchmarkResult {
        name: dataset_name,
        raw_size_mb: raw_mb,
        compressed_size_mb: comp_mb,
        comp_ratio: ratio,
        decomp_only_gb_s: decomp_only,
        needletail_gb_s: needletail_spd,
        paraseq_gb_s: paraseq_spd,
        helicase_str_gb_s: heli_str,
        helicase_pack_gb_s: heli_pack,
        helicase_col_gb_s: heli_col,
    }
}

fn main() {
    let temp_dir = std::env::temp_dir();

    println!("Generating benchmark datasets in RAM...");
    let datasets: Vec<(&'static str, Vec<u8>)> = vec![
        (
            "FASTQ_150bp_Illumina",
            FastxGen::illumina(50_000).with_seed(101).generate_bytes(),
        ),
        (
            "FASTQ_15kb_PacBio",
            FastxGen::pacbio(1_000).with_seed(202).generate_bytes(),
        ),
        (
            "FASTA_wrapped_60col",
            FastxGen::wrapped_fasta(150, 40_000, 60)
                .with_seed(303)
                .generate_bytes(),
        ),
        (
            "FASTA_unwrapped_500K",
            FastxGen::fasta(40, 500_000).with_seed(404).generate_bytes(),
        ),
    ];

    for format_type in [NifflerFormat::Gzip, NifflerFormat::Zstd] {
        let (fmt_label, ext) = match format_type {
            NifflerFormat::Gzip => ("GZIP (gz, flate2/zlib-rs level 5)", "gz"),
            NifflerFormat::Zstd => ("Zstandard (zstd level 3)", "zst"),
            _ => unreachable!(),
        };

        println!(
            "\n===================================================================================================================="
        );
        println!("BENCHMARK: {fmt_label}");
        println!(
            "===================================================================================================================="
        );
        println!(
            "{:<22} | {:>7} -> {:>7} ({:>4}x) | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10}",
            "Dataset",
            "Raw",
            "Comp",
            "Ratio",
            "DecompOnly",
            "Needletail",
            "Paraseq",
            "Heli (str)",
            "Heli (pk)",
            "Heli (col)"
        );
        println!(
            "--------------------------------------------------------------------------------------------------------------------"
        );

        for (name, raw_bytes) in &datasets {
            let res = run_compressed_suite(name, raw_bytes, format_type, ext, &temp_dir);
            println!(
                "{:<22} | {:>5.1} MB -> {:>5.1} MB ({:>4.1}x) | {:>8.2} GB/s | {:>8.2} GB/s | {:>8.2} GB/s | {:>8.2} GB/s | {:>8.2} GB/s | {:>8.2} GB/s",
                res.name,
                res.raw_size_mb,
                res.compressed_size_mb,
                res.comp_ratio,
                res.decomp_only_gb_s,
                res.needletail_gb_s,
                res.paraseq_gb_s,
                res.helicase_str_gb_s,
                res.helicase_pack_gb_s,
                res.helicase_col_gb_s,
            );
        }
    }
}
