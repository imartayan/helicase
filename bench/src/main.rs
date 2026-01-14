#[cfg(target_os = "linux")]
mod linux_perf;
mod measurement;

use helicase::config::advanced::*;
use helicase::config::*;
use helicase::input::*;
use helicase::*;
use measurement::{BaseTime, Measurement};

use clap::Parser as ClapParser;
use needletail::{parse_fastx_file, parse_fastx_reader};
use paraseq::{Record, fastx};
use std::fs::{metadata, read};

#[derive(ClapParser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA/FASTQ)
    input: String,
    /// Number of repetitions
    #[arg(short, long, default_value_t = 3)]
    repeat: u64,
    /// Disable perf metrics (Linux only)
    #[arg(short = 'P', long)]
    no_perf: bool,
    /// Disable baseline (needletail & paraseq)
    #[arg(short = 'B', long)]
    no_baseline: bool,
    /// Disable slice bench
    #[arg(short = 'S', long)]
    no_slice: bool,
    /// Enable (compressed) file bench
    #[arg(short, long, default_value_t = false)]
    file: bool,
    /// Enable mmap bench
    #[arg(short, long, default_value_t = false)]
    mmap: bool,
    /// Show DNA length
    #[arg(short = 'l', long, default_value_t = false)]
    show_len: bool,
}

const MINIMAL: Config = ParserOptions::default()
    .ignore_headers()
    .ignore_dna()
    .config();
const DNA_LEN: Config = COMPUTE_DNA_LEN | SPLIT_NON_ACTG | MERGE_DNA_CHUNKS | MERGE_RECORDS;
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

fn run_bench<M: Measurement>(args: &Args) {
    let path = &args.input;
    let size = metadata(path).expect("Cannot get file metadata").len();
    let mut input_file = FileInput::new(path).expect("Cannot open file");
    let compressed = input_file.is_compressed().unwrap();
    let data = if !args.no_baseline {
        read(path).expect("Cannot open file")
    } else {
        Vec::new()
    };
    let data = data.as_slice();
    let mut record_len = 0;
    let mut dna_len = 0;
    let mut m = M::new();

    if !args.no_baseline {
        if !args.no_slice {
            m.start();
            for _ in 0..args.repeat {
                let mut reader = parse_fastx_reader(data).expect("Failed to parse slice");
                dna_len = 0;
                while let Some(r) = reader.next() {
                    let record = r.expect("Invalid record");
                    let clean_seq = record.seq();
                    dna_len += clean_seq.len();
                }
            }
            let result = if args.show_len { Some(dna_len) } else { None };
            m.show("Needletail string \t(reader)", size, args.repeat, result);

            if !compressed {
                m.start();
                for _ in 0..args.repeat {
                    // let mut reader = fastx::Reader::new(data).expect("Failed to parse slice"); // crashes on human genome
                    let mut reader =
                        fastx::Reader::new_with_batch_size(data, 1).expect("Failed to parse slice");
                    let mut record_set = reader.new_record_set();
                    dna_len = 0;
                    while record_set.fill(&mut reader).unwrap() {
                        for r in record_set.iter() {
                            let record = r.expect("Invalid record");
                            let clean_seq = record.seq();
                            dna_len += clean_seq.len();
                        }
                    }
                }
                let result = if args.show_len { Some(dna_len) } else { None };
                m.show("Paraseq string  \t(reader)", size, args.repeat, result);
            }
        }
        if args.file {
            m.start();
            for _ in 0..args.repeat {
                let mut reader = parse_fastx_file(path).expect("Failed to parse file");
                dna_len = 0;
                while let Some(r) = reader.next() {
                    let record = r.expect("Invalid record");
                    let clean_seq = record.seq();
                    dna_len += clean_seq.len();
                }
            }
            let result = if args.show_len { Some(dna_len) } else { None };
            m.show("Needletail string \t(file)", size, args.repeat, result);

            m.start();
            for _ in 0..args.repeat {
                // let mut reader = fastx::Reader::from_path_with_batch_size(path).expect("Failed to parse file"); // crashes on human genome
                let mut reader = fastx::Reader::from_path_with_batch_size(path, 1)
                    .expect("Failed to parse file");
                let mut record_set = reader.new_record_set();
                dna_len = 0;
                while record_set.fill(&mut reader).unwrap() {
                    for r in record_set.iter() {
                        let record = r.expect("Invalid record");
                        let clean_seq = record.seq();
                        dna_len += clean_seq.len();
                    }
                }
            }
            let result = if args.show_len { Some(dna_len) } else { None };
            m.show("Paraseq string  \t(file)", size, args.repeat, result);
        }
    }

    if !args.no_slice && !compressed {
        m.start();
        for _ in 0..args.repeat {
            let parser = FastxParser::<MINIMAL>::from_slice(data);
            record_len = 0;
            for _ in parser {
                record_len += 1;
            }
        }
        let result = if args.show_len {
            Some(record_len)
        } else {
            None
        };
        m.show("Helicase #records \t(slice)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_LEN>::from_slice(data);
            dna_len = 0;
            parser.next();
            dna_len += parser.get_dna_len();
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase #bases \t(slice)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_STRING>::from_slice(data);
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_string().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase string \t(slice)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_PACKED>::from_slice(data);
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_packed().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase packed \t(slice)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_COLUMNAR>::from_slice(data);
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_columnar().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase columnar \t(slice)", size, args.repeat, result);
    }
    if args.file {
        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_STRING>::from_file(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_string().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase string \t(file)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser = FastxParser::<DNA_PACKED>::from_file(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_packed().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase packed \t(file)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser =
                FastxParser::<DNA_COLUMNAR>::from_file(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_columnar().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase columnar \t(file)", size, args.repeat, result);
    }
    if args.mmap && !compressed {
        m.start();
        for _ in 0..args.repeat {
            let mut parser =
                FastxParser::<DNA_STRING>::from_file_mmap(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_string().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase string \t(mmap)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser =
                FastxParser::<DNA_PACKED>::from_file_mmap(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_packed().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase packed \t(mmap)", size, args.repeat, result);

        m.start();
        for _ in 0..args.repeat {
            let mut parser =
                FastxParser::<DNA_COLUMNAR>::from_file_mmap(path).expect("Cannot open file");
            dna_len = 0;
            while let Some(_) = parser.next() {
                dna_len += parser.get_dna_columnar().len();
            }
        }
        let result = if args.show_len { Some(dna_len) } else { None };
        m.show("Helicase columnar \t(mmap)", size, args.repeat, result);
    }
}

fn main() {
    let args = Args::parse();
    #[cfg(target_os = "linux")]
    {
        if args.no_perf {
            use linux_perf::PerfMeasurement;
            run_bench::<PerfMeasurement>(&args);
        } else {
            run_bench::<BaseTime>(&args);
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        run_bench::<BaseTime>(&args);
    }
}
