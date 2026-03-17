use super::*;
use crate::carrying_add::Carry;
use crate::config::{advanced::*, *};
use crate::input::*;
use crate::simd::extract_fasta_bitmask;

use core::fmt;
use core::marker::PhantomData;
use std::io;

pub(crate) struct FastaBitmask {
    pub open_bracket: u64,
    pub line_feeds: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
    pub mask_non_actg: u64,
    pub mask_n: u64,
}

#[derive(Default, PartialEq)]
pub struct FastaChunk {
    pub len: usize,
    pub header: u64,
    pub split: u64,
    pub is_dna: u64,
    pub two_bits: u128,
    pub high_bit: u64,
    pub low_bit: u64,
    pub mask_non_actg: u64,
    pub mask_n: u64,
}

impl Chunk for FastaChunk {}

pub struct FastaLexer<'a, const CONFIG: Config, I: InputData<'a>> {
    pub(crate) input: I,
    carry: Carry,
    _phantom: PhantomData<&'a [u8]>,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastaLexer<'a, CONFIG, I>
{
    fn from_input(input: I) -> io::Result<Self> {
        Ok(Self {
            input,
            carry: Carry::new(false),
            _phantom: PhantomData,
        })
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Lexer for FastaLexer<'a, CONFIG, I> {
    type Input = I;

    #[inline(always)]
    fn input(&self) -> &I {
        &self.input
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastaLexer<'a, CONFIG, I> {
    type Item = FastaChunk;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.input.next().map(|chunk| {
            let mask = extract_fasta_bitmask::<CONFIG>(chunk);
            let non_lf = !mask.line_feeds;
            let c = self.carry.add(mask.open_bracket, non_lf);
            let header = (c ^ non_lf) | mask.open_bracket;
            let is_dna = mask.is_dna & !header & non_lf;
            let split = if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                !header & !is_dna & non_lf
            } else {
                0
            };
            FastaChunk {
                len: chunk.len(),
                header,
                split,
                is_dna,
                two_bits: mask.two_bits,
                high_bit: mask.high_bit,
                low_bit: mask.low_bit,
                mask_non_actg: mask.mask_non_actg & !mask.line_feeds,
                mask_n: mask.mask_n,
            }
        })
    }
}

impl fmt::Display for FastaChunk {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let escape = !self.header & !self.is_dna & !self.split;
        for i in 0..self.len {
            let hi: bool = (self.high_bit >> i) & 1 != 0;
            let lo: bool = (self.low_bit >> i) & 1 != 0;
            let escape: bool = (escape >> i) & 1 != 0;
            let split: bool = (self.split >> i) & 1 != 0;
            let header: bool = (self.header >> i) & 1 != 0;
            match (hi, lo, escape, split, header) {
                (_, _, _, _, true) => write!(f, ">"),
                (_, _, _, true, _) => write!(f, "|"),
                (_, _, true, _, _) => write!(f, " "),
                (false, false, _, _, _) => write!(f, "A"),
                (true, false, _, _, _) => write!(f, "T"),
                (false, true, _, _, _) => write!(f, "C"),
                (true, true, _, _, _) => write!(f, "G"),
            }?
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const CONFIG_COLUMNAR: Config = ParserOptions::default().dna_columnar().config();

    #[inline]
    fn parse_and_format(fasta: &str) -> String {
        let f = FastaLexer::<CONFIG_COLUMNAR, _>::from_slice(fasta.as_bytes())
            .unwrap()
            .next()
            .unwrap();
        format!("{f}")
    }

    #[test]
    fn test_parse() {
        let fasta = ">head\nTTTCTtaAAAA\nAGAAAA\nACAA\n>h>h\nCTCTTANNAAA\nCAAAnAGCTTT";
        let expected = ">>>>>>TTTCTTAAAAA AGAAAA ACAA >>>>>CTCTTA||AAA CAAA|AGCTTT";
        assert_eq!(parse_and_format(fasta), expected);
    }
}
