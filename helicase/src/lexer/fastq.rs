use super::*;
use crate::input::*;
use crate::simd::{extract_bitmask, extract_fastq_bitmask};

use core::marker::PhantomData;
use std::io;

pub(crate) type FastqBitmask = Bitmask;
pub type FastqChunk = Chunk;

pub struct FastqLexer<'a, const CONFIG: Config, I: InputData<'a>> {
    pub(crate) input: I,
    _phantom: PhantomData<&'a [u8]>,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastqLexer<'a, CONFIG, I>
{
    fn from_input(input: I) -> io::Result<Self> {
        Ok(Self {
            input,
            _phantom: PhantomData,
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqLexer<'a, CONFIG, I> {
    /// Re-run the SIMD scan on the current block and return the resulting
    /// [`FastqChunk`] with the given `len`.  Used after `skip_to_newline` to
    /// refresh stale bitmask fields (`is_dna`, `newline`, etc.) that were left
    /// over from the quality block.
    #[inline(always)]
    pub(crate) fn scan_current_block(&self, len: usize) -> FastqChunk {
        FastqChunk::from_mask::<CONFIG>(
            len,
            extract_bitmask::<CONFIG, SCAN_MODE_QUALITY>(self.input.current_block()),
        )
    }

    #[inline(always)]
    pub(crate) fn next_in_header(&mut self) -> Option<io::Result<FastqChunk>> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                FastqChunk::from_mask::<CONFIG>(
                    chunk.len(),
                    extract_bitmask::<CONFIG, SCAN_MODE_HEADER>(chunk),
                )
            })
        })
    }

    #[inline(always)]
    pub(crate) fn next_in_non_dna(&mut self) -> Option<io::Result<FastqChunk>> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                FastqChunk::from_mask::<CONFIG>(
                    chunk.len(),
                    extract_bitmask::<CONFIG, SCAN_MODE_QUALITY>(chunk),
                )
            })
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastqLexer<'a, CONFIG, I> {
    type Item = io::Result<FastqChunk>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                FastqChunk::from_mask::<CONFIG>(chunk.len(), extract_fastq_bitmask::<CONFIG>(chunk))
            })
        })
    }
}
