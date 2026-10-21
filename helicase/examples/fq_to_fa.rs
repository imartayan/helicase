use helicase::input::*;
use helicase::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const CONFIG: Config = ParserOptions::default().config();

fn main() -> std::io::Result<()> {
    let path_str = std::env::args().nth(1).expect("No input file given");
    let path = Path::new(&path_str);
    let out_path = path.with_extension("fa");
    assert_ne!(path, out_path, "Input and output paths must be different!");

    let mut parser = FastxParser::<CONFIG>::from_file(path)?;
    let out_file = File::create(out_path)?;
    let mut writer = BufWriter::new(out_file);

    while let Some(event) = parser.next() {
        event?;
        writer.write_all(b">")?;
        writer.write_all(parser.get_header())?;
        writer.write_all(b"\n")?;
        writer.write_all(parser.get_dna_string())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
