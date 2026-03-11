use helicase::input::*;
use helicase::*;
use needletail::parse_fastx_file;

use simd_minimizers::packed_seq::{PackedSeqVec, Seq, SeqVec};

const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .dna_string()
    .and_dna_packed()
    .and_dna_columnar()
    .keep_non_actg()
    .config();

fn main() {
    let path = std::env::args().nth(1).expect("No input file given");
    let mut reader = parse_fastx_file(&path).expect("Cannot open file");
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
    while let Some(r) = reader.next() {
        let record = r.expect("Invalid record");
        let _event = parser.next();
        let line = record.start_line_number();

        let packed_from_ascii = PackedSeqVec::from_ascii(&record.seq());
        let packed_native = parser.get_packed_seq();

        let (left, right) = (packed_from_ascii.as_slice(), packed_native);
        let mut eq = true;
        if left.len() != right.len() {
            eq = false;
            eprintln!(
                "Len mismatch line {line}: N={}, H={}, NP={}, HP={}, HDP={}, HDS={}",
                record.seq().len(),
                parser.get_dna_string().len(),
                left.len(),
                right.len(),
                parser.get_dna_packed().len(),
                parser.get_dna_columnar().len(),
            );
        }
        let len = left.len().min(right.len());
        for i in (0..len).step_by(29) {
            let len = (len - i).min(29);
            let this = left.slice(i..i + len).as_u64();
            let that = right.slice(i..i + len).as_u64();
            if this != that {
                eq = false;
                eprintln!("Bases mismatch line {line} pos {i}:");
                eprintln!("NP: {:058b}", this);
                eprintln!("HP: {:058b}", that);
                break;
            }
        }
        if len < right.len() {
            eprintln!("Extra bases line {line}:");
            for i in (len..right.len()).step_by(29) {
                let len = (right.len() - i).min(29);
                let that = right.slice(i..i + len).as_u64();
                eprintln!("{:058b}", that);
            }
        }
        if !eq {
            return;
        }
    }
    eprintln!("Identical results ✓")
}
