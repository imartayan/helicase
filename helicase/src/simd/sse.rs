#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::x86_64::*;
use core::mem::transmute;

const GREATER_THAN: __m128i = unsafe { transmute([b'>'; 16]) };
const LINE_FEED: __m128i = unsafe { transmute([b'\n'; 16]) };
const ASCII_N: __m128i = unsafe { transmute([b'N'; 16]) };
const LUT_ACTG: __m128i = unsafe { transmute(*b"A_C_T_G_________") };

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m128i;
        let v1 = _mm_loadu_si128(ptr);
        let v2 = _mm_loadu_si128(ptr.add(1));
        let v3 = _mm_loadu_si128(ptr.add(2));
        let v4 = _mm_loadu_si128(ptr.add(3));

        let open_bracket = u8_mask(v1, v2, v3, v4, GREATER_THAN);
        let line_feeds = u8_mask(v1, v2, v3, v4, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
            bitpack_dna::<CONFIG>(v1, v2, v3, v4);

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
        let ptr = buf.as_ptr() as *const __m128i;
        let v1 = _mm_loadu_si128(ptr);
        let v2 = _mm_loadu_si128(ptr.add(1));
        let v3 = _mm_loadu_si128(ptr.add(2));
        let v4 = _mm_loadu_si128(ptr.add(3));

        let line_feeds = u8_mask(v1, v2, v3, v4, LINE_FEED);

        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
            bitpack_dna::<CONFIG>(v1, v2, v3, v4);

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
fn bitpack_dna<const CONFIG: Config>(
    v1: __m128i,
    v2: __m128i,
    v3: __m128i,
    v4: __m128i,
) -> (u64, u128, u64, u64, u64, u64) {
    unsafe {
        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;
        let mut mask_non_actg = 0;
        let mut mask_n = 0;

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2, mm_hi_3, mm_lo_3, mm_hi_4, mm_lo_4) = (
                _mm_movemask_epi8(_mm_slli_epi16(v1, 5)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v1, 6)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v2, 5)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v2, 6)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v3, 5)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v3, 6)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v4, 5)) as u16 as u64,
                _mm_movemask_epi8(_mm_slli_epi16(v4, 6)) as u16 as u64,
            );
            high_bit = mm_hi_1 | (mm_hi_2 << 16) | (mm_hi_3 << 32) | (mm_hi_4 << 48);
            low_bit = mm_lo_1 | (mm_lo_2 << 16) | (mm_lo_3 << 32) | (mm_lo_4 << 48);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            #[cfg(all(target_feature = "bmi2", not(feature = "no-pdep")))]
            {
                let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2, mm_hi_3, mm_lo_3, mm_hi_4, mm_lo_4) = (
                    _mm_movemask_epi8(_mm_slli_epi16(v1, 5)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v1, 6)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v2, 5)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v2, 6)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v3, 5)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v3, 6)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v4, 5)) as u16 as u64,
                    _mm_movemask_epi8(_mm_slli_epi16(v4, 6)) as u16 as u64,
                );
                let mm_hi_12 = mm_hi_1 | (mm_hi_2 << 16);
                let mm_lo_12 = mm_lo_1 | (mm_lo_2 << 16);
                let mm_hi_34 = mm_hi_3 | (mm_hi_4 << 16);
                let mm_lo_34 = mm_lo_3 | (mm_lo_4 << 16);
                let mm_1 = _pdep_u64(mm_hi_12, 0xAAAAAAAAAAAAAAAA)
                    | _pdep_u64(mm_lo_12, 0x5555555555555555);
                let mm_2 = _pdep_u64(mm_hi_34, 0xAAAAAAAAAAAAAAAA)
                    | _pdep_u64(mm_lo_34, 0x5555555555555555);
                two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
            #[cfg(any(not(target_feature = "bmi2"), feature = "no-pdep"))]
            {
                // Adapted from https://github.com/Daniel-Liu-c0deb0t/cute-nucleotides/commit/007164bce68f671188fa5c607982fbd306112cb3
                let (hi_1, lo_1, hi_2, lo_2, hi_3, lo_3, hi_4, lo_4) = (
                    _mm_slli_epi16(v1, 5),
                    _mm_slli_epi16(v1, 6),
                    _mm_slli_epi16(v2, 5),
                    _mm_slli_epi16(v2, 6),
                    _mm_slli_epi16(v3, 5),
                    _mm_slli_epi16(v3, 6),
                    _mm_slli_epi16(v4, 5),
                    _mm_slli_epi16(v4, 6),
                );
                let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2, mm_hi_3, mm_lo_3, mm_hi_4, mm_lo_4) = (
                    _mm_movemask_epi8(_mm_unpackhi_epi8(lo_1, hi_1)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpacklo_epi8(lo_1, hi_1)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpackhi_epi8(lo_2, hi_2)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpacklo_epi8(lo_2, hi_2)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpackhi_epi8(lo_3, hi_3)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpacklo_epi8(lo_3, hi_3)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpackhi_epi8(lo_4, hi_4)) as u16 as u64,
                    _mm_movemask_epi8(_mm_unpacklo_epi8(lo_4, hi_4)) as u16 as u64,
                );
                let mm_1 = mm_lo_1 | (mm_hi_1 << 16) | (mm_lo_2 << 32) | (mm_hi_2 << 48);
                let mm_2 = mm_lo_3 | (mm_hi_3 << 16) | (mm_lo_4 << 32) | (mm_hi_4 << 48);
                two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG) {
            let mask_two_bits = _mm_set1_epi8(0b110i8);
            let mask_upper = _mm_set1_epi8(0b11011111u8 as i8);
            let uv1 = _mm_and_si128(v1, mask_upper);
            let uv2 = _mm_and_si128(v2, mask_upper);
            let uv3 = _mm_and_si128(v3, mask_upper);
            let uv4 = _mm_and_si128(v4, mask_upper);

            let mask_actg = movemask_64(
                _mm_cmpeq_epi8(
                    _mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v1, mask_two_bits)),
                    uv1,
                ),
                _mm_cmpeq_epi8(
                    _mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v2, mask_two_bits)),
                    uv2,
                ),
                _mm_cmpeq_epi8(
                    _mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v3, mask_two_bits)),
                    uv3,
                ),
                _mm_cmpeq_epi8(
                    _mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v4, mask_two_bits)),
                    uv4,
                ),
            );
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask_n = u8_mask(uv1, uv2, uv3, uv4, ASCII_N);
            }
        }

        (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n)
    }
}

#[inline(always)]
fn movemask_64(v1: __m128i, v2: __m128i, v3: __m128i, v4: __m128i) -> u64 {
    unsafe {
        (_mm_movemask_epi8(v1) as u16 as u64)
            | ((_mm_movemask_epi8(v2) as u16 as u64) << 16)
            | ((_mm_movemask_epi8(v3) as u16 as u64) << 32)
            | ((_mm_movemask_epi8(v4) as u16 as u64) << 48)
    }
}

#[inline(always)]
pub fn u8_mask(v1: __m128i, v2: __m128i, v3: __m128i, v4: __m128i, vc: __m128i) -> u64 {
    unsafe {
        let a = _mm_movemask_epi8(_mm_cmpeq_epi8(v1, vc)) as u16 as u64;
        let b = _mm_movemask_epi8(_mm_cmpeq_epi8(v2, vc)) as u16 as u64;
        let c = _mm_movemask_epi8(_mm_cmpeq_epi8(v3, vc)) as u16 as u64;
        let d = _mm_movemask_epi8(_mm_cmpeq_epi8(v4, vc)) as u16 as u64;
        a | (b << 16) | (c << 32) | (d << 48)
    }
}
