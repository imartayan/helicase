//! Parser types and traits.

mod fasta;
mod fastq;
mod fastx;
#[cfg(feature = "paraseq")]
mod paraseq_traits;
mod traits;

use core::ops::Range;

pub use fasta::*;
pub use fastq::*;
pub use fastx::*;
pub use traits::*;

#[derive(Debug)]
pub enum Format {
    Fasta,
    Fastq,
}

#[derive(Debug)]
pub enum Event {
    Record(Range<usize>),
    DnaChunk(Range<usize>),
    IoError(std::io::Error),
}

/// Unrolled 4x bulk 2-bit bitpacking into a pre-allocated buffer pointer.
#[inline(always)]
pub(crate) unsafe fn bulk_pack_dna_blocks<const CONFIG: crate::config::Config>(
    data: &[u8],
    q: usize,
    out_ptr: *mut u64,
) {
    let q4 = q / 4;
    for k in 0..q4 {
        let base_offset = 256 * k;
        let m0 = crate::simd::extract_dna_bitmask::<CONFIG>(&data[base_offset..base_offset + 64]);
        let m1 =
            crate::simd::extract_dna_bitmask::<CONFIG>(&data[base_offset + 64..base_offset + 128]);
        let m2 =
            crate::simd::extract_dna_bitmask::<CONFIG>(&data[base_offset + 128..base_offset + 192]);
        let m3 =
            crate::simd::extract_dna_bitmask::<CONFIG>(&data[base_offset + 192..base_offset + 256]);
        unsafe {
            let p = out_ptr.add(8 * k);
            p.write(m0.two_bits as u64);
            p.add(1).write((m0.two_bits >> 64) as u64);
            p.add(2).write(m1.two_bits as u64);
            p.add(3).write((m1.two_bits >> 64) as u64);
            p.add(4).write(m2.two_bits as u64);
            p.add(5).write((m2.two_bits >> 64) as u64);
            p.add(6).write(m3.two_bits as u64);
            p.add(7).write((m3.two_bits >> 64) as u64);
        }
    }
    for k in (q4 * 4)..q {
        let block_start = 64 * k;
        let m = crate::simd::extract_dna_bitmask::<CONFIG>(&data[block_start..block_start + 64]);
        unsafe {
            let p = out_ptr.add(2 * k);
            p.write(m.two_bits as u64);
            p.add(1).write((m.two_bits >> 64) as u64);
        }
    }
}
