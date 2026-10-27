//! Lexer types, bitmasks, and scanning modes.

use crate::config::{advanced::*, *};
use crate::input::*;
use crate::parser::{DnaBitpackers, Result};
use crate::simd::extract_bitmask;
use core::marker::PhantomData;

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

pub struct Bitmask {
    pub line_feeds: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
    pub mask_non_actg: u64,
    pub mask_n: u64,
}

impl Bitmask {
    pub const ALL_DNA: Self = Self {
        line_feeds: 0,
        is_dna: !0,
        two_bits: 0,
        high_bit: 0,
        low_bit: 0,
        mask_non_actg: 0,
        mask_n: 0,
    };
}

#[derive(Copy, Clone, Default, PartialEq, Eq)]
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
    pub const fn from_mask<const CONFIG: Config>(len: usize, mask: Bitmask) -> Self {
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

    #[inline(always)]
    pub const fn set<const CONFIG: Config>(&mut self, other: Self) {
        self.len = other.len;
        self.newline = other.newline;
        self.is_dna = other.is_dna;
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.two_bits = other.two_bits;
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.high_bit = other.high_bit;
            self.low_bit = other.low_bit;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg = other.mask_non_actg;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n = other.mask_n;
        }
    }
}

pub struct Lexer<'a, const CONFIG: Config, I: InputData<'a>> {
    pub(crate) input: I,
    last_found_line_feed: usize,
    _phantom: PhantomData<&'a [u8]>,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I> for Lexer<'a, CONFIG, I> {
    fn from_input(input: I) -> Result<Self> {
        Ok(Self {
            input,
            last_found_line_feed: 0,
            _phantom: PhantomData,
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Lexer<'a, CONFIG, I> {
    #[inline(always)]
    pub(crate) fn scan_current_block(&self, len: usize) -> Chunk {
        Chunk::from_mask::<CONFIG>(
            len,
            extract_bitmask::<CONFIG, SCAN_MODE_DNA>(self.input.current_block()),
        )
    }

    // Block containing `pos`, in the current buffer, with the line feeds before it cleared.
    #[inline(always)]
    pub(crate) fn seek_chunk(&mut self, pos: usize) -> Chunk {
        let block = pos / 64;
        self.input.seek_block(block);
        let end = self.input.buffer_offset() + self.input.buffer().len();
        let mut chunk = self.scan_current_block(64.min(end - 64 * block));
        chunk.newline &= !0 << (pos % 64);
        chunk
    }

    #[inline(always)]
    pub(crate) fn next_block<const MODE: u8>(&mut self) -> Option<Result<Chunk>> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                Chunk::from_mask::<CONFIG>(chunk.len(), extract_bitmask::<CONFIG, MODE>(chunk))
            })
        })
    }

    #[inline(always)]
    pub(crate) fn skip_long_dna(
        &mut self,
        bitpackers: &mut DnaBitpackers,
        dna_len: &mut usize,
        block_counter: &mut usize,
        anchor: usize,
        keep_anchor: bool,
    ) {
        let mut curr_pos = (*block_counter + 1) * 64;
        loop {
            let buf_off = self.input.buffer_offset();
            let data = self.input.buffer();
            let data_len = data.len();
            if curr_pos < buf_off + data_len {
                let local_curr = curr_pos - buf_off;
                let nl_opt = if self.last_found_line_feed >= curr_pos
                    && self.last_found_line_feed < buf_off + data_len
                {
                    Some(self.last_found_line_feed - curr_pos)
                } else {
                    let nl_opt = memchr::memchr(b'\n', &data[local_curr..]);
                    if let Some(nl_off) = nl_opt {
                        self.last_found_line_feed = curr_pos + nl_off;
                    }
                    nl_opt
                };
                let rem_len = match nl_opt {
                    Some(nl_off) => nl_off,
                    None if I::RANDOM_ACCESS => data_len - local_curr,
                    None => (data_len - local_curr).saturating_sub(64),
                };
                let aligned_len = rem_len & !63;
                if aligned_len > 0 {
                    let packed_bytes =
                        bitpackers.pack_contiguous::<CONFIG>(&data[local_curr..], aligned_len);
                    if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                        *dna_len += packed_bytes;
                    }
                    curr_pos += packed_bytes;
                    *block_counter = curr_pos / 64 - 1;
                    self.input.seek_block(*block_counter);
                    if packed_bytes < aligned_len {
                        return;
                    }
                }
                if nl_opt.is_some() || I::RANDOM_ACCESS {
                    return;
                }
            }
            if I::RANDOM_ACCESS {
                return;
            }
            let prev_end = buf_off + data_len;
            self.input
                .make_room(if keep_anchor { anchor } else { curr_pos });
            if self.input.buffer_offset() + self.input.buffer().len() <= prev_end {
                return;
            }
        }
    }

    #[inline(always)]
    pub(crate) fn skip_quality_line(
        &mut self,
        target_pos: usize,
        block_counter: &mut usize,
        pos_in_block: &mut usize,
        block: &mut Chunk,
        anchor: usize,
        keep_anchor: bool,
    ) {
        if !I::RANDOM_ACCESS {
            let mut curr_pos = (*block_counter + 1) * 64;
            while target_pos >= self.input.buffer_offset() + self.input.buffer().len() {
                let buf_off = self.input.buffer_offset();
                let data_len = self.input.buffer().len();
                let skip_bytes = data_len.saturating_sub(curr_pos - buf_off + 64) & !63;
                if skip_bytes > 0 {
                    curr_pos += skip_bytes;
                    *block_counter = curr_pos / 64 - 1;
                    self.input.seek_block(*block_counter);
                }
                let prev_end = buf_off + data_len;
                self.input
                    .make_room(if keep_anchor { anchor } else { curr_pos });
                if self.input.buffer_offset() + self.input.buffer().len() <= prev_end {
                    return;
                }
            }
        }
        let buf_off = self.input.buffer_offset();
        let data = self.input.buffer();
        let data_len = data.len();
        if target_pos < buf_off + data_len && data[target_pos - buf_off] == b'\n' {
            *block_counter = target_pos / 64;
            *pos_in_block = target_pos % 64;
            self.input.seek_block(*block_counter);
            let block_len = 64.min(buf_off + data_len - *block_counter * 64);
            block.set::<CONFIG>(Chunk::from_mask::<CONFIG>(
                block_len,
                extract_bitmask::<CONFIG, SCAN_MODE_QUALITY>(self.input.current_block()),
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const CONFIG_COLUMNAR: Config = ParserOptions::default().dna_columnar().config();

    #[test]
    fn test_lexer() {
        let fasta = b">head\nACGT\n";
        let chunk = Lexer::<CONFIG_COLUMNAR, _>::from_slice(fasta)
            .unwrap()
            .next_block::<SCAN_MODE_DNA>()
            .unwrap()
            .unwrap();
        assert_eq!(chunk.len, 11);
        assert_eq!(chunk.newline, (1 << 5) | (1 << 10));
    }
}
