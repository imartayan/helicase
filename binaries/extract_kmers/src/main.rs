mod timer;
use timer::TraceTimer;
mod arrow_types;
use arrow_types::ArrowDispatch;
mod mmap_arrow_reader;

mod arrow_writer;
use arrow_writer::*;

use clap::Parser as ClapParser;
use clap::ValueEnum;
use helicase::*;
use std::fmt::Display;
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
    #[arg(short = 'w', long)]
    dedup: bool,
    /// sort instead of using a hashset
    #[arg(short = 'h', long, default_value_t = false)]
    sort: bool,
    /// par_sort instead of using a hashset
    #[arg(short = 'p', long, default_value_t = false)]
    par_sort: bool,
    /// threads number
    #[arg(long, default_value_t = 4)]
    thread_number: usize,
    /// slow printing kmers
    #[arg(long, default_value_t = false)]
    print_kmer: bool,
    /// log level (warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
    /// will create a directory with arrow file fragment there
    #[arg(long)]
    arrow_path: Option<String>,
    /// log bucket nb
    #[arg(short = 'b', long, default_value_t = 8)]
    bucket_log_nb: usize,
    /// bucket capacity
    #[arg(short='c', long, default_value_t=1<<20)]
    capacity: usize,
}

fn main_dispatched<const K: usize, T: ArrowDispatch + Send + 'static + Display + Sync>(
    args: &Args,
) {
    let path = &args.input;
    let data = read(path).expect("Cannot open file");
    let _span = info_span!("main");
    let _t = TraceTimer::new("main");
    if let Some(path) = &args.arrow_path {
        fastx_slice_to_arrow::<K, T>(&data, path, args.bucket_log_nb, args.capacity).unwrap();
    }
    //if args.dedup | args.sort | args.par_sort {
    //    let mut kmers = MerChunk::<K, T>::new();
    //    if args.sort | args.par_sort {
    //        kmers = kmer_from_fastx_slice::<K, T>(&data).expect("invalid data");
    //        info!("Computing {} kmers", kmers.len());
    //        info!("Sorting kmers");
    //        if args.par_sort {
    //            kmers = kmers.par_sort(8, true, args.thread_number);
    //        } else {
    //            kmers.sort(true);
    //        }
    //        info!("Got {} uniq kmers", kmers.len());
    //    } else {
    //        use std::collections::HashSet;
    //        use std::hash::RandomState;
    //        info!("Building hset of kmers");
    //        let mut hash_set: HashSet<_, RandomState> = HashSet::from_iter(kmers.iter());
    //        chunk_process_from_fastx_slice::<K, T>(&data, |mer_chunk| {
    //            hash_set.extend(mer_chunk.iter());
    //            Ok(())
    //        })
    //        .unwrap();
    //        info!("Got {} uniq kmers", hash_set.len());
    //        kmers.extend(hash_set);
    //        info!("Building the KmerChunk from the hashset");
    //    }
    //    if let Some(path) = &args.parquet_path {
    //        mer_chunk_to_parquet::<K, T>(path, &mut kmers).unwrap();
    //    }
    //    if let Some(path) = &args.arrow_path {
    //        mer_chunk_to_arrow::<K, T>(path, &mut kmers).unwrap();
    //    }
    //    if args.print_kmer {
    //        for kmer in kmers.iter() {
    //            println!("{}", kmer.to_string());
    //        }
    //    }
    //} else {
    //    if let Some(path) = &args.parquet_path {
    //        todo!();
    //        //fastx_slice_to_parquet::<K, T>(path, &data).unwrap();
    //    }
    //    if let Some(path) = &args.arrow_path {
    //        fastx_slice_to_arrow::<K, T>(path, &data).unwrap();
    //    }
    //    if args.print_kmer {
    //        chunk_process_from_fastx_slice::<K, T>(&data, |mer_slice| {
    //            for kmer in mer_slice.iter() {
    //                println!("{}", kmer.to_string());
    //            }
    //            Ok(())
    //        })
    //        .unwrap();
    //    }
    //}
    info!("Done!");
}

fn main() {
    let args = Args::parse();
    init_tracing(&args.log_level);
    dispatch_k!(args.k, |K, T| main_dispatched::<K, T>(&args))
}
