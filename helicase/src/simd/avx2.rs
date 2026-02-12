#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::x86_64::*;
use core::mem::transmute;

const GREATER_THAN: __m256i = unsafe { transmute([b'>'; 32]) };
const LINE_FEED: __m256i = unsafe { transmute([b'\n'; 32]) };
const ASCII_N: __m256i = unsafe { transmute([b'N'; 32]) };
const LUT_ACTG: __m256i = unsafe { transmute(*b"A_C_T_G_________A_C_T_G_________") };

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m256i;
        let v1 = _mm256_loadu_si256(ptr);
        let v2 = _mm256_loadu_si256(ptr.add(1));

        let open_bracket = u8_mask(v1, v2, GREATER_THAN);
        let line_feeds = u8_mask(v1, v2, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
            bitpack_dna::<CONFIG>(v1, v2);

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
        let ptr = buf.as_ptr() as *const __m256i;
        let v1 = _mm256_loadu_si256(ptr);
        let v2 = _mm256_loadu_si256(ptr.add(1));

        let line_feeds = u8_mask(v1, v2, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
            bitpack_dna::<CONFIG>(v1, v2);

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
fn bitpack_dna<const CONFIG: Config>(v1: __m256i, v2: __m256i) -> (u64, u128, u64, u64, u64, u64) {
    unsafe {
        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;
        let mut mask_non_actg = 0;
        let mut mask_n = 0;

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2) = (
                _mm256_movemask_epi8(_mm256_slli_epi16(v1, 5)) as u32 as u64,
                _mm256_movemask_epi8(_mm256_slli_epi16(v1, 6)) as u32 as u64,
                _mm256_movemask_epi8(_mm256_slli_epi16(v2, 5)) as u32 as u64,
                _mm256_movemask_epi8(_mm256_slli_epi16(v2, 6)) as u32 as u64,
            );
            high_bit = mm_hi_1 | (mm_hi_2 << 32);
            low_bit = mm_lo_1 | (mm_lo_2 << 32);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            #[cfg(all(target_feature = "bmi2", not(feature = "no-pdep")))]
            {
                let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2) = (
                    _mm256_movemask_epi8(_mm256_slli_epi16(v1, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v1, 6)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v2, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v2, 6)) as u32 as u64,
                );
                let mm_1 =
                    _pdep_u64(mm_hi_1, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_1, 0x5555555555555555);
                let mm_2 =
                    _pdep_u64(mm_hi_2, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_2, 0x5555555555555555);
                two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
            #[cfg(any(not(target_feature = "bmi2"), feature = "no-pdep"))]
            {
                // Adapted from https://github.com/Daniel-Liu-c0deb0t/cute-nucleotides/commit/007164bce68f671188fa5c607982fbd306112cb3
                let (iv1, iv2) = (
                    _mm256_permute4x64_epi64(v1, 0xD8),
                    _mm256_permute4x64_epi64(v2, 0xD8),
                );
                let (hi_1, lo_1, hi_2, lo_2) = (
                    _mm256_slli_epi16(iv1, 5),
                    _mm256_slli_epi16(iv1, 6),
                    _mm256_slli_epi16(iv2, 5),
                    _mm256_slli_epi16(iv2, 6),
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
            let mask_two_bits = _mm256_set1_epi8(0b110i8);
            let mask_upper = _mm256_set1_epi8(0b11011111u8 as i8);
            let uv1 = _mm256_and_si256(v1, mask_upper);
            let uv2 = _mm256_and_si256(v2, mask_upper);

            let mask_actg = movemask_64(
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v1, mask_two_bits)),
                    uv1,
                ),
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v2, mask_two_bits)),
                    uv2,
                ),
            );
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask_n = u8_mask(v1, v2, ASCII_N);
            }
        }

        (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n)
    }
}

#[inline(always)]
fn movemask_64(v1: __m256i, v2: __m256i) -> u64 {
    unsafe {
        (_mm256_movemask_epi8(v1) as u32 as u64) | ((_mm256_movemask_epi8(v2) as u32 as u64) << 32)
    }
}

#[inline(always)]
pub fn u8_mask(v1: __m256i, v2: __m256i, vc: __m256i) -> u64 {
    unsafe {
        let cmp_c1 = _mm256_cmpeq_epi8(v1, vc);
        let cmp_c2 = _mm256_cmpeq_epi8(v2, vc);
        let a = _mm256_movemask_epi8(cmp_c1) as u32 as u64;
        let b = _mm256_movemask_epi8(cmp_c2) as u32 as u64;
        a | (b << 32)
    }
}
