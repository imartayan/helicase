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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
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
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    Record(Range<usize>),
    DnaChunk(Range<usize>),
}
