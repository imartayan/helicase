use super::{FastxParser, Format, HelicaseParser, Result};
use crate::config::*;
use crate::input::*;
use std::borrow::Cow;
use std::io;
use std::path::Path;

/// Default configuration for Needletail compatibility (headers, DNA strings, and quality enabled).
pub const NEEDLETAIL_CONFIG: Config = ParserOptions::default().compute_quality().config();

/// Needletail compatibility reader: [`next`](NeedletailReader::next) yields a [`SequenceRecord`].
pub struct NeedletailReader<'a>(FastxParser<'a, NEEDLETAIL_CONFIG>);

/// Needletail compatibility record, a view of the parser's current record.
pub struct SequenceRecord<'a>(&'a FastxParser<'a, NEEDLETAIL_CONFIG>);

impl<'a> NeedletailReader<'a> {
    /// Advance to the next record, like needletail's `FastxReader::next`.
    #[allow(clippy::should_implement_trait)]
    #[inline(always)]
    pub fn next(&mut self) -> Option<Result<SequenceRecord<'_>>> {
        match self.0.next()? {
            Ok(_) => Some(Ok(SequenceRecord(&self.0))),
            Err(e) => Some(Err(e)),
        }
    }

    /// Get the underlying [`FastxParser`].
    pub fn into_inner(self) -> FastxParser<'a, NEEDLETAIL_CONFIG> {
        self.0
    }
}

impl<'a> SequenceRecord<'a> {
    #[inline(always)]
    pub fn id(&self) -> &'a [u8] {
        self.0.get_header()
    }

    #[inline(always)]
    pub fn seq(&self) -> Cow<'a, [u8]> {
        Cow::Borrowed(self.0.get_dna_string())
    }

    /// Unlike needletail, wrapped FASTA lines are already joined.
    #[inline(always)]
    pub fn raw_seq(&self) -> &'a [u8] {
        self.0.get_dna_string()
    }

    #[inline(always)]
    pub fn qual(&self) -> Option<&'a [u8]> {
        self.0.get_quality()
    }

    #[inline(always)]
    pub fn num_bases(&self) -> usize {
        self.0.get_dna_len()
    }

    #[inline(always)]
    pub fn format(&self) -> Format {
        self.0.format()
    }
}

/// Needletail compatibility alias: parse a FASTA/FASTQ file with default options.
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<NeedletailReader<'static>> {
    FastxParser::from_file(path).map(NeedletailReader)
}

/// Needletail compatibility alias: parse from an arbitrary reader with default options.
pub fn parse_fastx_reader<R: io::Read + Send + 'static>(
    reader: R,
) -> Result<NeedletailReader<'static>> {
    FastxParser::from_reader(reader).map(NeedletailReader)
}

/// Needletail compatibility alias: parse from an in-memory byte slice with default options.
pub fn parse_fastx_slice(slice: &[u8]) -> Result<NeedletailReader<'_>> {
    FastxParser::from_slice(slice).map(NeedletailReader)
}
