use helicase::input::*;
use helicase::*;
use needletail::parse_fastx_file;
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::fs::metadata;
use std::time::Instant;

const CONFIG_STRING: Config = ParserOptions::default().ignore_headers().config();

const CONFIG_PACKED: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .keep_non_actg()
    .config();

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).expect("No input file given");
    let size = metadata(&path).expect("Cannot get file metadata").len() as usize;
    let _ = std::fs::read(&path)?.last();
    let reader = parse_fastx_file(&path);
    let parser_string = FastxParser::<CONFIG_STRING>::from_file(&path);
    let parser_packed = FastxParser::<CONFIG_PACKED>::from_file(&path);
    let parser_packed_par = ParallelFastxParser::<CONFIG_PACKED>::from_file(&path);

    let k = 21;
    let w = 11;
    let builder = minimizers(k, w);

    let buffer_size = size * 2 / (w + 1) * 11 / 10;
    let mut min_pos = Vec::with_capacity(buffer_size);

    if let Ok(mut reader) = reader {
        let now = Instant::now();
        while let Some(r) = reader.next() {
            let record = r.expect("Invalid record");
            let packed_seq = PackedSeqVec::from_ascii(&record.seq());
            let seq = packed_seq.as_slice();
            builder.run(seq, &mut min_pos);
        }
        eprintln!(
            "Needletail & compute minimizers:\t {:6.3} GB/s (single thread)",
            size as f64 / 1e9 / now.elapsed().as_secs_f64()
        );
        min_pos.clear();
    }

    if let Ok(mut parser_string) = parser_string {
        let now = Instant::now();
        while let Some(event) = parser_string.next() {
            event?;
            let packed_seq = PackedSeqVec::from_ascii(parser_string.get_dna_string());
            let seq = packed_seq.as_slice();
            builder.run(seq, &mut min_pos);
        }
        eprintln!(
            "Parse string & compute minimizers:\t {:6.3} GB/s (single thread)",
            size as f64 / 1e9 / now.elapsed().as_secs_f64()
        );
        min_pos.clear();
    }

    if let Ok(mut parser_packed) = parser_packed {
        let now = Instant::now();
        while let Some(event) = parser_packed.next() {
            event?;
            let seq = parser_packed.get_packed_seq();
            builder.run(seq, &mut min_pos);
        }
        eprintln!(
            "Parse packed & compute minimizers:\t {:6.3} GB/s (single thread)",
            size as f64 / 1e9 / now.elapsed().as_secs_f64()
        );
        min_pos.clear();
    }

    if let Ok(parser_packed_par) = parser_packed_par {
        let threads = 4;
        let now = Instant::now();
        parser_packed_par.threads(threads).for_each(|mut parser| {
            let mut min_pos = Vec::with_capacity(buffer_size / threads);
            while let Some(event) = parser.next() {
                event?;
                let seq = parser.get_packed_seq();
                builder.run(seq, &mut min_pos);
            }
            Ok(())
        })?;
        eprintln!(
            "Par parse packed & compute minimizers:\t {:6.3} GB/s ({threads} threads)",
            size as f64 / 1e9 / now.elapsed().as_secs_f64()
        );
    }

    Ok(())
}
