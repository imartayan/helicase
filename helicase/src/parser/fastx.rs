use super::*;
use crate::config::*;
use crate::dna_format::*;
use crate::input::*;
use std::io;

/// A wrapper for [`FastaParser`] / [`FastqParser`] detecting the format at runtime.
pub struct FastxParser<'a, const CONFIG: Config>(Box<dyn ParserIter + 'a>);

impl<'a, const CONFIG: Config, I: InputData<'a> + 'a> FromInputData<'a, I>
    for FastxParser<'a, CONFIG>
{
    fn from_input(input: I) -> io::Result<Self> {
        match input.first_byte() {
            b'>' => Ok(Self(Box::new(FastaParser::<CONFIG, I>::from_input(input)?))),
            b'@' => Ok(Self(Box::new(FastqParser::<CONFIG, I>::from_input(input)?))),
            _ => Err(io::Error::other(
                "Invalid record start, expected '>' or '@'",
            )),
        }
    }
}

impl<'a, const CONFIG: Config> Parser for FastxParser<'a, CONFIG> {
    #[inline(always)]
    fn format(&self) -> Format {
        self.0.format()
    }

    #[inline(always)]
    fn clear_record(&mut self) {
        self.0.clear_record();
    }

    #[inline(always)]
    fn clear_chunk(&mut self) {
        self.0.clear_chunk();
    }

    #[inline(always)]
    fn get_header(&self) -> &[u8] {
        self.0.get_header()
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        self.0.get_header_owned()
    }

    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        self.0.get_quality()
    }

    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        self.0.get_quality_owned()
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        self.0.get_dna_string()
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        self.0.get_dna_string_owned()
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        self.0.get_dna_columnar()
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        self.0.get_dna_columnar_owned()
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        self.0.get_dna_packed()
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        self.0.get_dna_packed_owned()
    }

    #[inline(always)]
    fn get_mask_non_actg(&self) -> &BitMask {
        self.0.get_mask_non_actg()
    }

    #[inline(always)]
    fn get_mask_non_actg_owned(&mut self) -> BitMask {
        self.0.get_mask_non_actg_owned()
    }

    #[inline(always)]
    fn get_mask_n(&self) -> &BitMask {
        self.0.get_mask_n()
    }

    #[inline(always)]
    fn get_mask_n_owned(&mut self) -> BitMask {
        self.0.get_mask_n_owned()
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        self.0.get_dna_len()
    }
}

impl<'a, const CONFIG: Config> Iterator for FastxParser<'a, CONFIG> {
    type Item = Event;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}
