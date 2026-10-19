use clap::Parser;
use fastx_gen::{Alphabet, FastxGen, Format};
use std::io::{self, BufWriter};

#[derive(Parser, Debug)]
#[command(author, version, about = "Synthetic FASTA/FASTQ generator")]
struct Cli {
    /// Output format (fasta or fastq)
    #[arg(short, long, default_value = "fastq")]
    format: String,

    /// Number of records to generate
    #[arg(short = 'n', long, default_value_t = 100)]
    records: usize,

    /// Mean sequence length in base pairs
    #[arg(short = 'l', long, default_value_t = 150)]
    len: usize,

    /// Relative standard deviation of sequence length (e.g. 0.1 for 10% variance, 0.0 for fixed)
    #[arg(short = 's', long, default_value_t = 0.0)]
    std_dev: f64,

    /// Random seed
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Wrap width for sequences (0 for unwrapped single-line)
    #[arg(short = 'w', long, default_value_t = 0)]
    wrap: usize,

    /// Alphabet to use: dna, dna-n, iupac, protein, mixed
    #[arg(short = 'a', long, default_value = "dna")]
    alphabet: String,

    /// Record header prefix
    #[arg(long, default_value = "record_")]
    prefix: String,

    /// Frequency / proportion of non-ACTG bases (N) to inject (0.0 to 1.0)
    #[arg(short = 'N', long, default_value_t = 0.0)]
    non_actg_prob: f64,

    /// Emit Windows CRLF (\r\n) line endings
    #[arg(long)]
    crlf: bool,
}

fn main() -> io::Result<()> {
    let args = Cli::parse();

    let format = match args.format.to_lowercase().as_str() {
        "fa" | "fasta" => Format::Fasta,
        "fq" | "fastq" => Format::Fastq,
        other => {
            eprintln!("Unknown format '{}', expected 'fasta' or 'fastq'", other);
            std::process::exit(1);
        }
    };

    let alphabet = match args.alphabet.to_lowercase().as_str() {
        "dna" => Alphabet::Dna,
        "dna-n" | "dnan" => Alphabet::DnaWithN,
        "iupac" => Alphabet::DnaIupac,
        "protein" | "prot" | "aa" => Alphabet::Protein,
        "mixed" => Alphabet::DnaMixed,
        other => {
            eprintln!(
                "Unknown alphabet '{}', expected dna, dna-n, iupac, protein, mixed",
                other
            );
            std::process::exit(1);
        }
    };

    let mut generator = match format {
        Format::Fasta => FastxGen::fasta(args.records, args.len),
        Format::Fastq => FastxGen::fastq(args.records, args.len),
    }
    .with_seed(args.seed)
    .with_std_dev(args.std_dev)
    .with_alphabet(alphabet);

    if args.non_actg_prob > 0.0 {
        generator = generator.with_non_actg_prob(args.non_actg_prob);
    }
    if args.wrap > 0 {
        generator = generator.with_wrap(args.wrap);
    }
    if args.crlf {
        generator.crlf = true;
    }
    generator.header_prefix = args.prefix;

    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    generator.write_to(&mut writer)?;

    Ok(())
}
