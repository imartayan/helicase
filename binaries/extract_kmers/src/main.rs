use clap::{Parser, ValueEnum};
use helicase::*;
use std::fs::read;
use std::time::Instant;

#[derive(Debug, Clone, ValueEnum)]
enum TypeLayout {
    U8,
    U16,
    U32,
    U64,
    U128,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA/FASTQ)
    input: String,
    /// kmer size
    #[arg(short, long)]
    k: usize,
    /// type layout
    #[arg(short = 'c', long)]
    t: TypeLayout,
}
fn main() {
    let args = Args::parse();
    let path = &args.input;
    let start = Instant::now();
    let data = read(path).expect("Cannot open file");
    let elapsed = start.elapsed().as_secs_f64();
    println!("Loading the data in {}s", elapsed);
    let start = Instant::now();
    let kmers = kmer_from_fastx_slice::<63, u64>(&data).expect("invalid data");
    let elapsed = start.elapsed().as_secs_f64();
    println!("Computing kmers in {}s", elapsed);
    println!("{}", kmers.len());
}
