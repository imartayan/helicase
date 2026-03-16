//! Parser types and traits.

mod fasta;
mod fastq;
mod fastx;
#[cfg(feature = "paraseq")]
mod paraseq_traits;
mod traits;

use core::ops::Range;

pub use fasta::*;
pub use fastq::*;
pub use fastx::*;
pub use traits::*;

#[derive(Debug)]
pub enum Format {
    Fasta,
    Fastq,
}

pub enum Event {
    Record(Range<usize>),
    DnaChunk(Range<usize>),
}
