#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::x86_64::*;
use core::mem::transmute;

const GREATER_THAN: __m512i = unsafe { transmute([b'>'; 64]) };
const LINE_FEED: __m512i = unsafe { transmute([b'\n'; 64]) };
const ASCII_N: __m512i = unsafe { transmute([b'N'; 64]) };
const LUT_ACTG: __m512i =
    unsafe { transmute(*b"A_C_T_G_________A_C_T_G_________A_C_T_G_________A_C_T_G_________") };

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m512i;
        let v = _mm512_loadu_si512(ptr);

        let open_bracket = _mm512_cmpeq_epi8_mask(v, GREATER_THAN);
        let line_feeds = _mm512_cmpeq_epi8_mask(v, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) = bitpack_dna::<CONFIG>(v);

        FastaBitmask {
            open_bracket,
            line_feeds,
            is_dna,
            two_bits,
            high_bit,
            low_bit,
            mask_non_actg,
            mask_n,
        }
    }
}

#[inline(always)]
pub fn extract_fastq_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastqBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m512i;
        let v = _mm512_loadu_si512(ptr);

        let line_feeds = _mm512_cmpeq_epi8_mask(v, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) = bitpack_dna::<CONFIG>(v);

        FastqBitmask {
            line_feeds,
            is_dna,
            two_bits,
            high_bit,
            low_bit,
            mask_non_actg,
            mask_n,
        }
    }
}

#[inline(always)]
fn bitpack_dna<const CONFIG: Config>(v: __m512i) -> (u64, u128, u64, u64, u64, u64) {
    unsafe {
        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;
        let mut mask_non_actg = 0;
        let mut mask_n = 0;

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit = _mm512_movepi8_mask(_mm512_slli_epi16(v, 5));
            low_bit = _mm512_movepi8_mask(_mm512_slli_epi16(v, 6));
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            #[cfg(all(target_feature = "bmi2", not(feature = "no-pdep")))]
            {
                let hi = _mm512_movepi8_mask(_mm512_slli_epi16(v, 5));
                let lo = _mm512_movepi8_mask(_mm512_slli_epi16(v, 6));
                let mm_1 = _pdep_u64(hi & 0xFFFFFFFF, 0xAAAAAAAAAAAAAAAA)
                    | _pdep_u64(lo & 0xFFFFFFFF, 0x5555555555555555);
                let mm_2 = _pdep_u64(hi >> 32, 0xAAAAAAAAAAAAAAAA)
                    | _pdep_u64(lo >> 32, 0x5555555555555555);
                two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
            #[cfg(any(not(target_feature = "bmi2"), feature = "no-pdep"))]
            {
                // Adapted from https://github.com/Daniel-Liu-c0deb0t/cute-nucleotides/commit/007164bce68f671188fa5c607982fbd306112cb3
                let v_lo = _mm512_extracti64x4_epi64(v, 0);
                let v_hi = _mm512_extracti64x4_epi64(v, 1);
                let (iv_lo, iv_hi) = (
                    _mm256_permute4x64_epi64(v_lo, 0xD8),
                    _mm256_permute4x64_epi64(v_hi, 0xD8),
                );
                let (hi_1, lo_1, hi_2, lo_2) = (
                    _mm256_slli_epi16(iv_lo, 5),
                    _mm256_slli_epi16(iv_lo, 6),
                    _mm256_slli_epi16(iv_hi, 5),
                    _mm256_slli_epi16(iv_hi, 6),
                );
                let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2) = (
                    _mm256_movemask_epi8(_mm256_unpackhi_epi8(lo_1, hi_1)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_unpacklo_epi8(lo_1, hi_1)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_unpackhi_epi8(lo_2, hi_2)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_unpacklo_epi8(lo_2, hi_2)) as u32 as u64,
                );
                let mm_1 = (mm_hi_1 << 32) | mm_lo_1;
                let mm_2 = (mm_hi_2 << 32) | mm_lo_2;
                two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG) {
            let mask_two_bits = _mm512_set1_epi8(0b110i8);
            let mask_upper = _mm512_set1_epi8(0b11011111u8 as i8);
            let uv = _mm512_and_si512(v, mask_upper);

            let mask_actg = _mm512_cmpeq_epi8_mask(
                _mm512_shuffle_epi8(LUT_ACTG, _mm512_and_si512(v, mask_two_bits)),
                uv,
            );
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask_n = u8_mask(v, ASCII_N);
            }
        }

        (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n)
    }
}

#[inline(always)]
pub fn u8_mask(v: __m512i, vc: __m512i) -> u64 {
    unsafe { _mm512_cmpeq_epi8_mask(v, vc) }
}
