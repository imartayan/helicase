use clap::Parser;
use eyre::{Result, eyre};
use helicase::*;
use memmap2::Mmap;
use std::fs::read;
use tracing::{info, info_span};
use tracing_subscriber::fmt::time::UtcTime;
use tracing_subscriber::{EnvFilter, fmt};

fn init_tracing(log_level: &str) {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(log_level));

    fmt()
        .with_env_filter(filter)
        .with_timer(UtcTime::rfc_3339())
        .with_target(false)
        .init();
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA/FASTQ)
    input: String,
    /// kmer size (must be <=32)
    #[arg(short, long)]
    k: usize,
    /// log level (warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
    /// Prevent usage of mmap
    #[arg(long)]
    no_mmap: bool,
    /// Exact with hash map
    #[arg(long)]
    exact: bool,
    /// Do not normalize kmers
    #[arg(long)]
    no_normalize: bool,
}

fn kmer_approximate_count<const K: usize, T: BitStorage>(
    data: &[u8],
    args: &Args,
) -> Result<usize> {
    use cardinality_estimator::CardinalityEstimator;
    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .split_non_actg()
        .return_dna_chunk(true)
        .return_record(false)
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;

    let mut estimator = CardinalityEstimator::<u64>::default();
    while parser.next().is_some() {
        let seq = parser.get_dna_columnar();
        for mut val in Mer::<K, T>::from_dna_columns(seq) {
            if !args.no_normalize {
                val = val.normalize()
            }
            estimator.insert(&val.hash());
        }
    }
    Ok(estimator.estimate())
}

fn kmer_exact_count<const K: usize, T: BitStorage>(data: &[u8], args: &Args) -> Result<usize> {
    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .split_non_actg()
        .return_dna_chunk(true)
        .return_record(false)
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;

    let mut estimator = std::collections::HashSet::<Mer<K, T>>::new();
    while parser.next().is_some() {
        let seq = parser.get_dna_columnar();
        for mut val in Mer::<K, T>::from_dna_columns(seq) {
            if !args.no_normalize {
                val = val.normalize()
            }
            estimator.insert(val);
        }
    }
    Ok(estimator.len())
}

fn main_dispatched<const K: usize, T: BitStorage>(args: &Args) -> Result<()> {
    let path = &args.input;
    let _span = info_span!("main");
    let exec = |data| {
        if args.exact {
            kmer_exact_count::<K, T>(data, args)
        } else {
            kmer_approximate_count::<K, T>(data, args)
        }
    };
    let approx = if args.no_mmap {
        info!("Loading data in RAM");
        let data = read(path)?;
        exec(&data)?
    } else {
        info!("Access to data with mmap");
        let file = std::fs::File::open(path)?;
        let data = unsafe { Mmap::map(&file)? };
        exec(&data)?
    };
    info!("uniq kmer:{approx}");
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.k > 32 || args.k < 7 {
        return Err(eyre!("Support only k<=32, got '{}'", args.k));
    }
    init_tracing(&args.log_level);
    dispatch_k!(args.k, |K, T| main_dispatched::<K, T>(&args))
}
