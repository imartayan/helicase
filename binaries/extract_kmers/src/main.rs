mod timer;
use timer::TraceTimer;
mod arrow_types;
use arrow_types::ArrowDispatch;

mod arrow_builders;
mod arrow_writer;
use arrow_writer::*;

use clap::Parser as ClapParser;
use clap::ValueEnum;
use eyre::{Result, eyre};
use helicase::*;
use memmap2::Mmap;
use std::fmt::Display;
use std::fs::read;
use std::path::Path;
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

#[derive(ClapParser, Debug, Clone, ValueEnum)]
enum TypeLayout {
    U8,
    U16,
    U32,
    U64,
    U128,
}

#[derive(ClapParser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA/FASTQ)
    input: String,
    /// kmer size
    #[arg(short, long)]
    k: usize,
    /// without deduplication
    #[arg(long, default_value_t = false)]
    no_dedup: bool,
    /// log level (warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
    /// will create a directory with arrow file fragment there
    #[arg(short = 'o', long)]
    output: String,
    /// log bucket nb
    #[arg(short = 'b', long, default_value_t = 4)]
    bucket_log_nb: usize,
    /// bucket starting capacity
    #[arg(short='c', long, default_value_t=1<<23)]
    capacity: usize,
    /// numbers of threads. Carefull, high parallelism
    /// will be detrimental as it will staturate IO.
    /// Also, RAM behaviour are per threads. Hence
    /// the number of threads is multiplicative with
    /// respect to RAM usage.
    #[arg(long, default_value_t = 4)]
    max_num_threads: usize,
    /// keep the directory and prevent the final concatenation of buckets
    #[arg(long, default_value_t = false)]
    keep_bucket: bool,
    /// The Hash-Sort aggregation needs space. This number of power of two we use as a multiplying
    /// factor. For bucket of size n, the allocated space will be
    /// 2^(log2(n) + hashsort_log_factor)
    #[arg(long, default_value_t = 1)]
    hashsort_log_factor: u32,
    /// Prevent usage of mmap
    #[arg(long)]
    no_mmap: bool,
}

fn main_dispatched<const K: usize, T: ArrowDispatch + Send + 'static + Display + Sync>(
    args: &Args,
) -> Result<()> {
    let path = &args.input;
    let _span = info_span!("main");
    let exec = |data| {
        fastx_slice_to_arrow::<K, T>(
            data,
            &args.output,
            args.bucket_log_nb,
            args.capacity,
            args.hashsort_log_factor,
            !args.no_dedup,
            !args.keep_bucket,
        )
    };
    if args.no_mmap {
        info!("Loading data in RAM");
        let data = read(path)?;
        exec(&data)
    } else {
        info!("Access to data with mmap");
        let file = std::fs::File::open(path)?;
        let data = unsafe { Mmap::map(&file)? };
        exec(&data)
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.max_num_threads)
        .build_global()?;

    if Path::new(&args.output).exists() {
        return Err(eyre!("Output file '{}' already exists", &args.output));
    }
    init_tracing(&args.log_level);
    dispatch_k!(args.k, |K, T| main_dispatched::<K, T>(&args))
}
