use super::*;
use crate::dna_format::*;

#[cfg(feature = "packed-seq")]
use packed_seq::{PackedSeq, PackedSeqVec};

/// Getters for the current record, implemented by all parsers.
pub trait HelicaseParser: crate::sealed::Sealed {
    /// Get the [`Format`] associated to this parser (FASTA or FASTQ).
    fn format(&self) -> Format;

    /// Get a reference to the current header.
    fn get_header(&self) -> &[u8];

    /// Get an owned version of the current header.
    /// This will trigger a new allocation and a copy.
    fn get_header_owned(&mut self) -> Vec<u8>;

    /// Get a reference to the current sequence as a slice of bytes.
    fn get_dna_string(&self) -> &[u8];

    /// Get an owned version of the current sequence as a `Vec<u8>`.
    /// This will trigger a new allocation and possibly a copy.
    fn get_dna_string_owned(&mut self) -> Vec<u8>;

    /// Get a reference to the current sequence as [`ColumnarDNA`].
    fn get_dna_columnar(&self) -> &ColumnarDNA;

    /// Get an owned version of the current sequence as [`ColumnarDNA`].
    /// This allocates a buffer of the same capacity for the parser, see [`swap_dna_columnar`](Self::swap_dna_columnar) to reuse one.
    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        let mut owned = ColumnarDNA::with_capacity(self.get_dna_columnar().capacity());
        self.swap_dna_columnar(&mut owned);
        owned
    }

    /// Swap the current sequence as [`ColumnarDNA`] with `other`, without allocating.
    /// The parser keeps parsing into the buffer of `other`, which is cleared at the next record
    /// (or appended to with [`merge_records`](crate::ParserOptions::merge_records)).
    fn swap_dna_columnar(&mut self, other: &mut ColumnarDNA);

    /// Get a reference to the current sequence as [`PackedDNA`].
    fn get_dna_packed(&self) -> &PackedDNA;

    /// Get an owned version of the current sequence as [`PackedDNA`].
    /// This allocates a buffer of the same capacity for the parser, see [`swap_dna_packed`](Self::swap_dna_packed) to reuse one.
    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        let mut owned = PackedDNA::with_capacity(self.get_dna_packed().capacity());
        self.swap_dna_packed(&mut owned);
        owned
    }

    /// Swap the current sequence as [`PackedDNA`] with `other`, without allocating.
    /// The parser keeps parsing into the buffer of `other`, which is cleared at the next record
    /// (or appended to with [`merge_records`](crate::ParserOptions::merge_records)).
    fn swap_dna_packed(&mut self, other: &mut PackedDNA);

    /// Get a [`PackedSeq`] based on the underlying [`PackedDNA`].
    ///
    /// Note this is equivalent to:
    /// ```ignore
    /// parser.get_dna_packed().as_packed_seq()
    /// ```
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    fn get_packed_seq(&self) -> PackedSeq<'_> {
        self.get_dna_packed().as_packed_seq()
    }

    /// Get a [`PackedSeqVec`] based on the underlying [`PackedDNA`].
    /// This will trigger a new allocation.
    ///
    /// Note this is equivalent to:
    /// ```ignore
    /// parser.get_dna_packed_owned().into_packed_seq_vec()
    /// ```
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    fn get_packed_seq_vec(&mut self) -> PackedSeqVec {
        self.get_dna_packed_owned().into_packed_seq_vec()
    }

    /// Swap the current sequence as [`PackedSeqVec`] with `other`, like [`swap_dna_packed`](Self::swap_dna_packed).
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    fn swap_packed_seq_vec(&mut self, other: &mut PackedSeqVec) {
        let mut packed = PackedDNA::from_packed_seq_vec(core::mem::take(other));
        self.swap_dna_packed(&mut packed);
        *other = packed.into_packed_seq_vec();
    }

    /// Get a reference to the bitmask marking non-ACTG bases.
    fn get_mask_non_actg(&self) -> &BitMask;

    /// Get an owned version of the bitmask marking non-ACTG bases.
    /// This allocates a buffer of the same capacity for the parser, see [`swap_mask_non_actg`](Self::swap_mask_non_actg) to reuse one.
    #[inline(always)]
    fn get_mask_non_actg_owned(&mut self) -> BitMask {
        let mut owned = BitMask::with_capacity(self.get_mask_non_actg().capacity());
        self.swap_mask_non_actg(&mut owned);
        owned
    }

    /// Swap the bitmask marking non-ACTG bases with `other`, like [`swap_dna_packed`](Self::swap_dna_packed).
    fn swap_mask_non_actg(&mut self, other: &mut BitMask);

    /// Get a reference to the bitmask marking N bases.
    fn get_mask_n(&self) -> &BitMask;

    /// Get an owned version of the bitmask marking N bases.
    /// This allocates a buffer of the same capacity for the parser, see [`swap_mask_n`](Self::swap_mask_n) to reuse one.
    #[inline(always)]
    fn get_mask_n_owned(&mut self) -> BitMask {
        let mut owned = BitMask::with_capacity(self.get_mask_n().capacity());
        self.swap_mask_n(&mut owned);
        owned
    }

    /// Swap the bitmask marking N bases with `other`, like [`swap_dna_packed`](Self::swap_dna_packed).
    fn swap_mask_n(&mut self, other: &mut BitMask);

    /// Get the length of the current sequence, or the total length so far with
    /// [`merge_records`](crate::ParserOptions::merge_records).
    fn get_dna_len(&self) -> usize;

    /// Get a reference to the current quality line.
    /// Returns `None` for FASTA.
    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        None
    }

    /// Get an owned version of the current quality line.
    /// This will trigger a new allocation and a copy.
    /// Returns `None` for FASTA.
    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        None
    }

    /// Manually clear the information of the current chunk.
    /// This is only useful when [`MERGE_DNA_CHUNKS`](crate::config::advanced::MERGE_DNA_CHUNKS)
    /// or [`merge_records`](crate::ParserOptions::merge_records) is enabled.
    fn clear_chunk(&mut self);

    /// Manually clear the information of the current record.
    /// This is only useful with [`merge_records`](crate::ParserOptions::merge_records).
    fn clear_record(&mut self);

    /// Get the 0-based index of the current record.
    #[cfg(feature = "paraseq")]
    fn record_index(&self) -> u64;
}

pub(crate) trait HelicaseParserIter:
    HelicaseParser + Iterator<Item = Result<Event>>
{
}

impl<T: HelicaseParser + Iterator<Item = Result<Event>>> HelicaseParserIter for T {}
