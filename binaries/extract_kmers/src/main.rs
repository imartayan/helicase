mod parquet_writer;
use parquet_writer::*;

use clap::Parser as ClapParser;
use helicase::*;
use std::fs::read;
use tracing::*;
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

#[derive(ClapParser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA/FASTQ)
    input: String,
    /// kmer size
    #[arg(short, long)]
    k: usize,
    /// without deduplication
    #[arg(short = 'w', long)]
    dedup: bool,
    /// sort instead of using a hashset
    #[arg(short = 'h', long, default_value_t = false)]
    sort: bool,
    /// slow printing kmers
    #[arg(short = 'p', long, default_value_t = false)]
    print_kmer: bool,
    /// log level (warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
    /// dump to parquet
    #[arg(long)]
    parquet_path: Option<String>,
}

fn main_dispatched<const K: usize, T: ArrowDispatch>(args: &Args) {
    let path = &args.input;
    let data = read(path).expect("Cannot open file");
    info!("Loading the data");
    if args.dedup {
        let mut kmers = kmer_from_fastx_slice::<K, T>(&data).expect("invalid data");
        info!("Computing {} kmers", kmers.len());
        if args.sort {
            info!("Sorting kmers");
            kmers.sort(true);
            info!("Got {} uniq kmers", kmers.len());
        } else {
            use std::collections::HashSet;
            use std::hash::RandomState;
            info!("Building hset of kmers");
            let hash_set: HashSet<_, RandomState> = HashSet::from_iter(kmers.iter());
            info!("Got {} uniq kmers", hash_set.len());
        }
        if let Some(path) = &args.parquet_path {
            mer_slice_to_parquet::<K, T>(path, kmers.as_slice()).unwrap();
        }
        if args.print_kmer {
            for kmer in kmers.iter() {
                println!("{}", kmer.to_string());
            }
        }
    } else {
        if let Some(path) = &args.parquet_path {
            fastx_slice_to_parquet::<K, T>(path, &data).unwrap();
        }
        if args.print_kmer {
            chunk_process_from_fastx_slice::<K, T>(&data, |mer_slice| {
                for kmer in mer_slice.iter() {
                    println!("{}", kmer.to_string());
                }
                Ok(())
            })
            .unwrap();
        }
    }
}

fn main() {
    let args = Args::parse();
    init_tracing(&args.log_level);
    dispatch_k!(args.k, |K, T| main_dispatched::<K, T>(&args))
}
