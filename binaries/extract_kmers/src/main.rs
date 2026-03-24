mod timer;
use timer::TraceTimer;
mod arrow_types;
use arrow_types::ArrowDispatch;

mod arrow_writer;
use arrow_writer::*;

mod kernels;

use clap::Parser as ClapParser;
use clap::ValueEnum;
use helicase::*;
use std::fmt::Display;
use std::fs::read;
use tracing::info_span;
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
    /// without sorting
    #[arg(long, default_value_t = false)]
    no_sort: bool,
    /// log level (warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
    /// will create a directory with arrow file fragment there
    #[arg(short = 'o', long)]
    output: String,
    /// log bucket nb
    #[arg(short = 'b', long, default_value_t = 7)]
    bucket_log_nb: usize,
    /// bucket starting capacity
    #[arg(short='c', long, default_value_t=1<<10)]
    capacity: usize,
    /// keep the directory and prevent the final concatenation of buckets
    #[arg(long, default_value_t = false)]
    keep_bucket: bool,
}

fn main_dispatched<const K: usize, T: ArrowDispatch + Send + 'static + Display + Sync>(
    args: &Args,
) {
    let path = &args.input;
    let data = read(path).expect("Cannot open file");
    let _span = info_span!("main");
    let _t = TraceTimer::new("main");
    fastx_slice_to_arrow::<K, T>(
        &data,
        &args.output,
        args.bucket_log_nb,
        args.capacity,
        !args.no_sort,
        !args.no_dedup,
        !args.keep_bucket,
    )
    .unwrap();
}

fn main() {
    let args = Args::parse();
    init_tracing(&args.log_level);
    dispatch_k!(args.k, |K, T| main_dispatched::<K, T>(&args))
}
