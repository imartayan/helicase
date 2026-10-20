use crate::config::*;
use crate::input::*;
use crate::parser::FastxParser;
use std::io;
use std::path::Path;

/// Default configuration for Needletail compatibility (headers, DNA strings, and quality enabled).
pub const NEEDLETAIL_CONFIG: Config = ParserOptions::default().compute_quality().config();

/// Needletail compatibility alias: parse a FASTA/FASTQ file with default options.
pub fn parse_fastx_file<P: AsRef<Path>>(
    path: P,
) -> io::Result<FastxParser<'static, NEEDLETAIL_CONFIG>> {
    FastxParser::from_file(path)
}

/// Needletail compatibility alias: parse from an arbitrary reader with default options.
pub fn parse_fastx_reader<R: io::Read + Send + 'static>(
    reader: R,
) -> io::Result<FastxParser<'static, NEEDLETAIL_CONFIG>> {
    FastxParser::from_reader(reader)
}

/// Needletail compatibility alias: parse from an in-memory byte slice with default options.
pub fn parse_fastx_slice<'a>(slice: &'a [u8]) -> io::Result<FastxParser<'a, NEEDLETAIL_CONFIG>> {
    FastxParser::from_slice(slice)
}
