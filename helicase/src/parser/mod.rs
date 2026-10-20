//! Parser types and traits.

mod bitpack;
mod fasta;
mod fastq;
mod fastx;
mod parallel;
#[cfg(feature = "paraseq")]
mod paraseq_traits;
mod traits;

use core::ops::Range;

pub(crate) use bitpack::*;
pub use fasta::*;
pub use fastq::*;
pub use fastx::*;
pub use parallel::ParallelFastxParser;
#[cfg(feature = "zstd")]
pub(crate) use parallel::index_zstd_frames;
pub use traits::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

impl Format {
    /// Identifies the format from the initial record character (`>` for FASTA, `@` for FASTQ).
    #[inline]
    pub fn from_byte(byte: u8) -> Option<Self> {
        match byte {
            b'>' => Some(Format::Fasta),
            b'@' => Some(Format::Fastq),
            _ => None,
        }
    }

    /// Identifies the format from `byte`, returning an [`io::Error`] if neither FASTA nor FASTQ.
    #[inline]
    pub fn try_from_byte(byte: u8) -> std::io::Result<Self> {
        Self::from_byte(byte)
            .ok_or_else(|| std::io::Error::other("Invalid record start, expected '>' or '@'"))
    }
}

#[derive(Debug)]
pub enum Event {
    Record(Range<usize>),
    DnaChunk(Range<usize>),
    IoError(std::io::Error),
}
