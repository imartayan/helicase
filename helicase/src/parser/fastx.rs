use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use alloc::boxed::Box;

/// A wrapper for [`FastaParser`] / [`FastqParser`] detecting the format at runtime.
/// Configured via [`ParserOptions`].
///
/// # Example
///
/// ```rust,no_run
/// use helicase::input::*;
/// use helicase::*;
///
/// // set the options of the parser (at compile-time)
/// const CONFIG: Config = ParserOptions::default().config();
///
/// # #[cfg(feature = "std")]
/// fn main() -> std::io::Result<()> {
///     let path = "...";
///
///     // create a parser with the desired options
///     let mut parser = FastxParser::<CONFIG>::from_file(&path)?;
///
///     // iterate over records
///     while let Some(event) = parser.next() {
///         event?;
///         let header = parser.get_header();
///         let seq = parser.get_dna_string();
///         // ...
///     }
///     Ok(())
/// }
/// # #[cfg(not(feature = "std"))]
/// # fn main() {}
/// ```
pub struct FastxParser<'a, const CONFIG: Config>(Box<dyn HelicaseParserIter + Send + 'a>);

impl<'a, const CONFIG: Config, I: InputData<'a> + Send + 'a> FromInputData<'a, I>
    for FastxParser<'a, CONFIG>
{
    fn from_input(input: I) -> Result<Self> {
        let format = if input.buffer().is_empty() {
            Format::Fasta
        } else {
            Format::try_from_byte(input.first_byte())?
        };
        match format {
            Format::Fasta => Ok(Self(Box::new(FastaParser::<CONFIG, I>::from_input(input)?))),
            Format::Fastq => Ok(Self(Box::new(FastqParser::<CONFIG, I>::from_input(input)?))),
        }
    }
}

// Checked at compile time, unlike the inner parsers whose methods are all instantiated by the vtable.
macro_rules! require {
    ($flag:expr, $msg:literal) => {
        const {
            assert!(
                flag_is_set(CONFIG, $flag),
                concat!("Parser config error: ", $msg)
            )
        }
    };
}

impl<const CONFIG: Config> crate::sealed::Sealed for FastxParser<'_, CONFIG> {}

impl<'a, const CONFIG: Config> HelicaseParser for FastxParser<'a, CONFIG> {
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
        require!(COMPUTE_HEADER, "headers are ignored");
        self.0.get_header()
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        require!(COMPUTE_HEADER, "headers are ignored");
        self.0.get_header_owned()
    }

    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        require!(COMPUTE_QUALITY, "quality is ignored");
        self.0.get_quality()
    }

    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        require!(COMPUTE_QUALITY, "quality is ignored");
        self.0.get_quality_owned()
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        require!(COMPUTE_DNA_STRING, "dna_string is not enabled");
        self.0.get_dna_string()
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        require!(COMPUTE_DNA_STRING, "dna_string is not enabled");
        self.0.get_dna_string_owned()
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        require!(COMPUTE_DNA_COLUMNAR, "dna_columnar is not enabled");
        self.0.get_dna_columnar()
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        require!(COMPUTE_DNA_COLUMNAR, "dna_columnar is not enabled");
        self.0.get_dna_columnar_owned()
    }

    #[inline(always)]
    fn swap_dna_columnar(&mut self, other: &mut ColumnarDNA) {
        require!(COMPUTE_DNA_COLUMNAR, "dna_columnar is not enabled");
        self.0.swap_dna_columnar(other)
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        require!(COMPUTE_DNA_PACKED, "dna_packed is not enabled");
        self.0.get_dna_packed()
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        require!(COMPUTE_DNA_PACKED, "dna_packed is not enabled");
        self.0.get_dna_packed_owned()
    }

    #[inline(always)]
    fn swap_dna_packed(&mut self, other: &mut PackedDNA) {
        require!(COMPUTE_DNA_PACKED, "dna_packed is not enabled");
        self.0.swap_dna_packed(other)
    }

    #[inline(always)]
    fn get_mask_non_actg(&self) -> &BitMask {
        require!(COMPUTE_MASK_NON_ACTG, "mask_non_actg is not enabled");
        self.0.get_mask_non_actg()
    }

    #[inline(always)]
    fn get_mask_non_actg_owned(&mut self) -> BitMask {
        require!(COMPUTE_MASK_NON_ACTG, "mask_non_actg is not enabled");
        self.0.get_mask_non_actg_owned()
    }

    #[inline(always)]
    fn swap_mask_non_actg(&mut self, other: &mut BitMask) {
        require!(COMPUTE_MASK_NON_ACTG, "mask_non_actg is not enabled");
        self.0.swap_mask_non_actg(other)
    }

    #[inline(always)]
    fn get_mask_n(&self) -> &BitMask {
        require!(COMPUTE_MASK_N, "mask_n is not enabled");
        self.0.get_mask_n()
    }

    #[inline(always)]
    fn get_mask_n_owned(&mut self) -> BitMask {
        require!(COMPUTE_MASK_N, "mask_n is not enabled");
        self.0.get_mask_n_owned()
    }

    #[inline(always)]
    fn swap_mask_n(&mut self, other: &mut BitMask) {
        require!(COMPUTE_MASK_N, "mask_n is not enabled");
        self.0.swap_mask_n(other)
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        const DNA: Config =
            COMPUTE_DNA_LEN | COMPUTE_DNA_STRING | COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED;
        require!(DNA, "dna is ignored");
        self.0.get_dna_len()
    }

    #[cfg(feature = "paraseq")]
    #[inline(always)]
    fn record_index(&self) -> u64 {
        self.0.record_index()
    }
}

impl<'a, const CONFIG: Config> Iterator for FastxParser<'a, CONFIG> {
    type Item = Result<Event>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}
