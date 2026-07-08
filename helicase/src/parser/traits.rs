use super::*;
use crate::dna_format::*;

#[cfg(feature = "packed-seq")]
use packed_seq::{PackedSeq, PackedSeqVec};

pub trait HelicaseParser {
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
    /// This will trigger a new allocation.
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA;

    /// Get a reference to the current sequence as [`PackedDNA`].
    fn get_dna_packed(&self) -> &PackedDNA;

    /// Get an owned version of the current sequence as [`PackedDNA`].
    /// This will trigger a new allocation.
    fn get_dna_packed_owned(&mut self) -> PackedDNA;

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

    /// Get a reference to the bitmask marking non-ACTG bases.
    fn get_mask_non_actg(&self) -> &BitMask;

    /// Get an owned version of the bitmask marking non-ACTG bases.
    /// This will trigger a new allocation.
    fn get_mask_non_actg_owned(&mut self) -> BitMask;

    /// Get a reference to the bitmask marking non-ACTG bases.
    fn get_mask_n(&self) -> &BitMask;

    /// Get an owned version of the bitmask marking non-ACTG bases.
    /// This will trigger a new allocation.
    fn get_mask_n_owned(&mut self) -> BitMask;

    /// Get the length of the current sequence.
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
    /// This is only useful when [`MERGE_DNA_CHUNKS`](crate::config::advanced::MERGE_DNA_CHUNKS) is enabled.
    fn clear_chunk(&mut self);

    /// Manually clear the information of the current record.
    /// This is only useful when [`MERGE_RECORDS`](crate::config::advanced::MERGE_RECORDS) is enabled.
    fn clear_record(&mut self);
}

pub trait HelicaseParserIter: HelicaseParser + Iterator<Item = Event> {}

impl<T: HelicaseParser + Iterator<Item = Event>> HelicaseParserIter for T {}
