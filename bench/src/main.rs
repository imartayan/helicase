mod measurement;
use measurement::{BaseTime, Measurement};

#[cfg(target_os = "linux")]
mod linux_perf;

use helicase::config::*;
use helicase::input::*;
use helicase::*;

use needletail::parse_fastx_reader;
use paraseq::{Record, fastx};

use std::env::args;
use std::fs::read;
use std::path::Path;

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

#[allow(unused)]
struct Setup<'a, P: AsRef<Path>> {
    path: P,
    data: &'a [u8],
    size: u64,
    rep: u64,
    compressed: bool,
}

fn measurment_variant<M: Measurement, P: AsRef<Path>>(s: Setup<P>) {
    let mut m = M::new();
    let mut dna_len = 0usize;

    #[cfg(feature = "baselines")]
    {
        // Needletail
        m.start();
        for _ in 0..s.rep {
            let mut reader = parse_fastx_reader(s.data).expect("invalid reader");
            dna_len = 0usize;
            while let Some(r) = reader.next() {
                let record = r.expect("invalid record");
                let clean_seq = record.seq();
                dna_len += clean_seq.len();
            }
        }
        m.show::<_>("Needletail", s.size, s.rep, dna_len);

        m.start();
        for _ in 0..s.rep {
            // let mut reader = fastx::Reader::new(data).expect("invalid reader"); // crashes on human genome
            let mut reader = fastx::Reader::new_with_batch_size(s.data, 1).expect("invalid reader");
            let mut record_set = reader.new_record_set();
            dna_len = 0;
            while record_set.fill(&mut reader).unwrap() {
                for r in record_set.iter() {
                    let record = r.expect("invalid record");
                    let clean_seq = record.seq();
                    dna_len += clean_seq.len();
                }
            }
        }
        m.show("Paraseq", s.size, s.rep, dna_len);
    }

    m.start();
    for _ in 0..s.rep {
        let mut parser = FastxParser::<DNA_STRING>::from_slice(s.data);
        dna_len = 0;
        while let Some(_) = parser.next() {
            dna_len += parser.get_dna_string().len();
        }
    }
    m.show("helicase (DNA string)", s.size, s.rep, dna_len);

    m.start();
    for _ in 0..s.rep {
        let mut parser = FastxParser::<DNA_PACKED>::from_slice(s.data);
        dna_len = 0;
        while let Some(_) = parser.next() {
            dna_len += parser.get_dna_packed().len();
        }
    }
    m.show("helicase (DNA packed)", s.size, s.rep, dna_len);

    m.start();
    for _ in 0..s.rep {
        let mut parser = FastxParser::<DNA_COLUMNAR>::from_slice(s.data);
        dna_len = 0;
        while let Some(_) = parser.next() {
            dna_len += parser.get_dna_columnar().len();
        }
    }
    m.show("helicase (DNA columnar)", s.size, s.rep, dna_len);
}

fn main() {
    let path = args().nth(1).expect("No input file given");
    let content = read(&path).expect("Cannot open file");
    let data = content.as_slice();
    let size = data.len() as u64;
    let mut input_file = FileInput::new(&path).expect("Cannot open file");
    let compressed = input_file.is_compressed().unwrap();
    let rep = 3;

    let s = Setup {
        path: &path,
        data,
        size,
        compressed,
        rep,
    };
    #[cfg(target_os = "linux")]
    {
        use linux_perf::PerfMeasurement;
        measurment_variant::<PerfMeasurement, _>(s);
    }
    #[cfg(not(target_os = "linux"))]
    {
        measurment_variant::<BaseTime, _>(s);
    }
}
