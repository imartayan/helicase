use super::*;
use crate::config::*;
use crate::input::*;

use alloc::borrow::Cow;
use paraseq::Record;

impl<'a, const CONFIG: Config, I: InputData<'a>> Record for FastaParser<'a, CONFIG, I> {
    #[inline(always)]
    fn id(&self) -> &[u8] {
        self.get_header()
    }
    #[inline(always)]
    fn seq(&self) -> Cow<'_, [u8]> {
        Cow::Borrowed(self.get_dna_string())
    }
    #[inline(always)]
    fn seq_raw(&self) -> &[u8] {
        self.get_dna_string()
    }
    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        self.get_quality()
    }
    #[inline(always)]
    fn index(&self) -> u64 {
        self.record_index()
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Record for FastqParser<'a, CONFIG, I> {
    #[inline(always)]
    fn id(&self) -> &[u8] {
        self.get_header()
    }
    #[inline(always)]
    fn seq(&self) -> Cow<'_, [u8]> {
        Cow::Borrowed(self.get_dna_string())
    }
    #[inline(always)]
    fn seq_raw(&self) -> &[u8] {
        self.get_dna_string()
    }
    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        self.get_quality()
    }
    #[inline(always)]
    fn index(&self) -> u64 {
        self.record_index()
    }
}

impl<'a, const CONFIG: Config> Record for FastxParser<'a, CONFIG> {
    #[inline(always)]
    fn id(&self) -> &[u8] {
        self.get_header()
    }
    #[inline(always)]
    fn seq(&self) -> Cow<'_, [u8]> {
        Cow::Borrowed(self.get_dna_string())
    }
    #[inline(always)]
    fn seq_raw(&self) -> &[u8] {
        self.get_dna_string()
    }
    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        self.get_quality()
    }
    #[inline(always)]
    fn index(&self) -> u64 {
        self.record_index()
    }
}
