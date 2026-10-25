//! Bitpacked DNA formats.
//!
//! All lengths, capacities and counts are in bases, whatever the number of bits per base.

mod columnar;
mod mask;
mod packed;

pub use columnar::*;
pub use mask::*;
pub use packed::*;
