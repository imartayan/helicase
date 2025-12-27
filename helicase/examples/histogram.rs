#![allow(unused)]

use helicase::input::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .dna_columnar()
    // don't stop the iterator at the end of a record
    .return_record(false)
    .config();

fn main() {
    let path = std::env::args().nth(1).expect("No input file given");

    // create a parser with the desired options
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
    let mut num_a = 0;
    let mut num_c = 0;
    let mut num_t = 0;
    let mut num_g = 0;

    while let Some(_event) = parser.next() {
        let cdna = parser.get_dna_columnar();
        let (high_bits, hi) = cdna.high_bits();
        let (low_bits, lo) = cdna.low_bits();
        let rem = cdna.len() % 64;
        for (&hi, &lo) in high_bits.into_iter().zip(low_bits) {
            num_a += (!hi & !lo).count_ones() as usize;
            num_c += (!hi & lo).count_ones() as usize;
            num_t += (hi & !lo).count_ones() as usize;
            num_g += (hi & lo).count_ones() as usize;
        }
        if rem > 0 {
            num_a += (!hi & !lo).count_ones() as usize - (64 - rem);
            num_c += (!hi & lo).count_ones() as usize;
            num_t += (hi & !lo).count_ones() as usize;
            num_g += (hi & lo).count_ones() as usize;
        }
    }

    println!("a: {num_a}");
    println!("c: {num_c}");
    println!("t: {num_t}");
    println!("g: {num_g}");
}
