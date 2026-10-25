//! Parser types and traits.

#[macro_use]
mod macros;
mod bitpack;
mod error;
mod fasta;
mod fastq;
mod fastx;
mod line_index;
#[cfg(feature = "needletail")]
pub(crate) mod needletail_helpers;
#[cfg(feature = "paraseq")]
mod paraseq_traits;
mod traits;

use alloc::vec::Vec;
use core::ops::Range;

pub(crate) use bitpack::*;
pub use error::{Error, ErrorKind};
pub use fasta::*;
pub use fastq::*;
pub use fastx::*;
use line_index::{LineIndex, scan_block};
#[cfg(feature = "needletail")]
pub use needletail_helpers::{NeedletailReader, SequenceRecord};
pub use traits::*;

/// A specialized [`Result`] type for Helicase parser operations.
pub type Result<T> = core::result::Result<T, Error>;

/// Sequence file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Records starting with `>`.
    Fasta,
    /// Records starting with `@`.
    Fastq,
}

impl Format {
    /// Identifies the format from the initial record character (`>` for FASTA, `@` for FASTQ).
    #[inline]
    pub const fn from_byte(byte: u8) -> Option<Self> {
        match byte {
            b'>' => Some(Self::Fasta),
            b'@' => Some(Self::Fastq),
            _ => None,
        }
    }

    /// Identifies the format from `byte`, returning an [`Error`] if neither FASTA nor FASTQ.
    #[inline]
    pub fn try_from_byte(byte: u8) -> Result<Self> {
        Self::from_byte(byte).ok_or_else(|| Error::invalid_record_start("'>' or '@'", byte))
    }

    /// Returns the index of the first record starting at or after `pos` (always after a newline),
    /// or `None` if there is none or a FASTQ header is cut by the end of `data`.
    pub fn find_record_start(self, data: &[u8], pos: usize) -> Option<usize> {
        let mut start = pos.saturating_sub(1);
        while let Some(newline) = memchr::memchr(b'\n', data.get(start..)?) {
            let cand = start + newline + 1;
            let is_valid = match self {
                Self::Fasta => data.get(cand) == Some(&b'>'),
                Self::Fastq => data.get(cand) == Some(&b'@') && is_fastq_header(&data[cand..])?,
            };
            if is_valid {
                return Some(cand);
            }
            start = cand;
        }
        None
    }
}

/// Item yielded by the parsers, see [`ParserOptions`](crate::ParserOptions) to choose which ones.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    /// End of a record, with its byte range in the input.
    Record(Range<usize>),
    /// Run of ACTG bases, with its byte range in the input.
    DnaChunk(Range<usize>),
}

// A line starting with `@` is a header if the next line does not start with `@` (quality otherwise).
// `None` if `line` ends before the first byte of the next line.
fn is_fastq_header(line: &[u8]) -> Option<bool> {
    let next = memchr::memchr(b'\n', line)? + 1;
    line.get(next).map(|&b| b != b'@')
}
