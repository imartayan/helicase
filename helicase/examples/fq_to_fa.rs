use helicase::input::*;
use helicase::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const CONFIG: Config = ParserOptions::default().config();

fn main() {
    let path_str = std::env::args().nth(1).expect("No input file given");
    let path = Path::new(&path_str);
    let out_path = path.with_extension("fa");
    assert_ne!(path, out_path, "Input and output paths must be different!");

    let mut parser = FastxParser::<CONFIG>::from_file(path).expect("Cannot open file");
    let out_file = File::create(out_path).expect("Cannot create output file");
    let mut writer = BufWriter::new(out_file);

    while let Some(_event) = parser.next() {
        writer.write_all(b">").unwrap();
        writer.write_all(parser.get_header()).unwrap();
        writer.write_all(b"\n").unwrap();
        writer.write_all(parser.get_dna_string()).unwrap();
        writer.write_all(b"\n").unwrap();
    }
}
