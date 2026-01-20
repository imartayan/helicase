use super::*;
use crate::config::*;
use crate::input::*;
use crate::simd::extract_fastq_bitmask;

use core::marker::PhantomData;
use std::io;

pub(crate) struct FastqBitmask {
    pub line_feeds: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
}

#[derive(Default)]
pub struct FastqChunk {
    pub len: usize,
    pub newline: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
}

impl Chunk for FastqChunk {}

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

impl<'a, const CONFIG: Config, I: InputData<'a>> Lexer for FastqLexer<'a, CONFIG, I> {
    type Input = I;

    #[inline(always)]
    fn input(&self) -> &I {
        &self.input
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastqLexer<'a, CONFIG, I> {
    type Item = FastqChunk;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.input.next().map(|chunk| {
            let mask = extract_fastq_bitmask::<CONFIG>(chunk);
            FastqChunk {
                len: chunk.len(),
                newline: mask.line_feeds,
                is_dna: mask.is_dna & !mask.line_feeds,
                two_bits: mask.two_bits,
                high_bit: mask.high_bit,
                low_bit: mask.low_bit,
            }
        })
    }
}
