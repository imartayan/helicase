use super::*;
use crate::input::*;
use crate::simd::{extract_bitmask, extract_fasta_bitmask};

use core::marker::PhantomData;
use std::io;

pub(crate) type FastaBitmask = Bitmask;
pub type FastaChunk = Chunk;

pub struct FastaLexer<'a, const CONFIG: Config, I: InputData<'a>> {
    pub(crate) input: I,
    _phantom: PhantomData<&'a [u8]>,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastaLexer<'a, CONFIG, I>
{
    fn from_input(input: I) -> io::Result<Self> {
        Ok(Self {
            input,
            _phantom: PhantomData,
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaLexer<'a, CONFIG, I> {
    #[inline(always)]
    pub(crate) fn scan_current_block(&self, len: usize) -> FastaChunk {
        FastaChunk::from_mask::<CONFIG>(
            len,
            extract_fasta_bitmask::<CONFIG>(self.input.current_block()),
        )
    }

    #[inline(always)]
    pub(crate) fn next_in_header(&mut self) -> Option<io::Result<FastaChunk>> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                FastaChunk::from_mask::<CONFIG>(
                    chunk.len(),
                    extract_bitmask::<CONFIG, SCAN_MODE_HEADER>(chunk),
                )
            })
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastaLexer<'a, CONFIG, I> {
    type Item = io::Result<FastaChunk>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.input.next().map(|res| {
            res.map(|chunk| {
                FastaChunk::from_mask::<CONFIG>(chunk.len(), extract_fasta_bitmask::<CONFIG>(chunk))
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const CONFIG_COLUMNAR: Config = ParserOptions::default().dna_columnar().config();

    #[test]
    fn test_lexer() {
        let fasta = b">head\nACGT\n";
        let chunk = FastaLexer::<CONFIG_COLUMNAR, _>::from_slice(fasta)
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(chunk.len, 11);
        assert_eq!(chunk.newline, (1 << 5) | (1 << 10));
    }
}
