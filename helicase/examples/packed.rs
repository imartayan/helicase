#![allow(unused)]

use helicase::input::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default()
    .dna_packed()
    // don't stop the iterator at the end of a record
    .return_record(false)
    .config();

fn main() {
    let path = std::env::args().nth(1).expect("No input file given");

    // create a parser with the desired options
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    // iterate over each chunk of ACTG
    while let Some(_event) = parser.next() {
        // we still have access to the header
        let header = parser.get_header();

        // get a reference to the packed sequence
        let seq = parser.get_dna_packed();

        // or directly get a PackedSeq (requires the packed-seq feature)
        let packed_seq = parser.get_packed_seq();

        // ...
    }
}
