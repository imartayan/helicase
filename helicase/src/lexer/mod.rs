//! Lexer types, bitmasks, and scanning modes.

use crate::config::{advanced::*, *};

mod fasta;
mod fastq;

pub use fasta::*;
pub use fastq::*;

/// Context-aware scanning modes for 64-byte SIMD kernels.
///
/// Controls when the SIMD kernel can safely skip expensive DNA 2-bit packing
/// and ACTG validation lookups.
pub mod scan_mode {
    /// Scanning DNA sequence.
    ///
    /// Never skips DNA bitpacking: every block is expected to contain sequence bases.
    pub const DNA: u8 = 0;

    /// Scanning record headers (`>` in FASTA or `@` in FASTQ).
    ///
    /// Skips DNA bitpacking as long as `line_feeds == 0` (the entire 64B block is header text).
    /// If a newline is found, the header ends and sequence begins in that block, so DNA
    /// is bitpacked for the remainder.
    pub const HEADER: u8 = 1;

    /// Scanning FASTQ non-DNA sections (`+` separator and Phred quality scores).
    ///
    /// Skips DNA bitpacking when `line_feeds` has at most 1 newline (`(line_feeds & (line_feeds - 1)) == 0`).
    /// - 0 newlines: within quality scores.
    /// - 1 newline: quality line ends, transitioning to the next record's `@header` (still not DNA).
    ///
    /// Only bitpacks if 2 or more newlines occur in the same 64B block (quality ends + tiny header ends + sequence starts).
    pub const QUALITY: u8 = 2;
}

pub use scan_mode::{
    DNA as SCAN_MODE_DNA, HEADER as SCAN_MODE_HEADER, QUALITY as SCAN_MODE_QUALITY,
};

pub(crate) struct Bitmask {
    pub line_feeds: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
    pub mask_non_actg: u64,
    pub mask_n: u64,
}

#[derive(Copy, Clone, Default, PartialEq)]
pub struct Chunk {
    pub len: usize,
    pub newline: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
    pub mask_non_actg: u64,
    pub mask_n: u64,
}

impl Chunk {
    #[inline(always)]
    pub(crate) fn from_mask<const CONFIG: Config>(len: usize, mask: Bitmask) -> Self {
        let is_dna = if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            mask.is_dna
        } else if len < 64 {
            !mask.line_feeds & (!0 >> (64 - len))
        } else {
            !mask.line_feeds
        };
        let mask_non_actg = if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            mask.mask_non_actg & !mask.line_feeds
        } else {
            0
        };
        Self {
            len,
            newline: mask.line_feeds,
            is_dna,
            two_bits: mask.two_bits,
            high_bit: mask.high_bit,
            low_bit: mask.low_bit,
            mask_non_actg,
            mask_n: mask.mask_n,
        }
    }
}
