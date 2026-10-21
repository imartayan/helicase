#![allow(unused)]

use helicase::input::*;
use helicase::*;

// set the options of the parser (at compile-time)
const CONFIG: Config = ParserOptions::default().config();

fn main() -> std::io::Result<()> {
    let path = std::env::args().nth(1).expect("No input file given");

    // create a parser with the desired options
    let mut parser = FastxParser::<CONFIG>::from_file(&path)?;

    // iterate over records
    while let Some(event) = parser.next() {
        event?;
        // get a reference to the header
        let header = parser.get_header();

        // get a reference to the sequence (without newlines)
        let seq = parser.get_dna_string();

        // ...
    }
    Ok(())
}
