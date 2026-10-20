use super::*;
use crate::config::*;
use crate::input::*;

use alloc::borrow::Cow;
use paraseq::Record;

macro_rules! impl_record {
    (($($generics:tt)*), $type:ty) => {
        impl<$($generics)*> Record for $type {
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

        impl<$($generics)*> Record for &$type {
            #[inline(always)]
            fn id(&self) -> &[u8] {
                (*self).get_header()
            }
            #[inline(always)]
            fn seq(&self) -> Cow<'_, [u8]> {
                Cow::Borrowed((*self).get_dna_string())
            }
            #[inline(always)]
            fn seq_raw(&self) -> &[u8] {
                (*self).get_dna_string()
            }
            #[inline(always)]
            fn qual(&self) -> Option<&[u8]> {
                (*self).get_quality()
            }
            #[inline(always)]
            fn index(&self) -> u64 {
                (*self).record_index()
            }
        }
    };
}

impl_record!(('a, const CONFIG: Config, I: InputData<'a>), FastaParser<'a, CONFIG, I>);
impl_record!(('a, const CONFIG: Config, I: InputData<'a>), FastqParser<'a, CONFIG, I>);
impl_record!(('a, const CONFIG: Config), FastxParser<'a, CONFIG>);
