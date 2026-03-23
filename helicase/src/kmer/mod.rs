mod available_kmer;
mod bitstring;
pub mod kmer_layouts;
pub mod extract;
pub use bitstring::{BitStorage, BitString, sliding_window};
pub use kmer_layouts::{Kmer, Mer};
