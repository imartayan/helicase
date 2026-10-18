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
pub use traits::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

#[derive(Debug)]
pub enum Event {
    Record(Range<usize>),
    DnaChunk(Range<usize>),
    IoError(std::io::Error),
}
