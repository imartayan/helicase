#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::x86_64::*;
use core::mem::transmute;

const LINE_FEED: __m512i = unsafe { transmute([b'\n'; 64]) };
const ASCII_N: __m512i = unsafe { transmute([b'N'; 64]) };
const LUT_ACTG: __m512i =
    unsafe { transmute(*b"A_C_T_G_________A_C_T_G_________A_C_T_G_________A_C_T_G_________") };

const NEEDS_DNA_MASKS: Config = COMPUTE_DNA_COLUMNAR
    | COMPUTE_DNA_PACKED
    | SPLIT_NON_ACTG
    | COMPUTE_MASK_NON_ACTG
    | COMPUTE_MASK_N;

#[inline(always)]
pub fn extract_bitmask<const CONFIG: Config, const MODE: u8>(buf: &[u8]) -> Bitmask {
    unsafe {
        let v = _mm512_loadu_si512(buf.as_ptr() as *const __m512i);

        if MODE == SCAN_MODE_DNA || flag_is_not_set(CONFIG, NEEDS_DNA_MASKS) {
            let mut mask = bitpack_dna::<CONFIG>(v);
            if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) || mask.is_dna != !0 {
                mask.line_feeds = u8_mask(v, LINE_FEED);
            }
            mask
        } else {
            let line_feeds = u8_mask(v, LINE_FEED);
            let skip_dna = if MODE == SCAN_MODE_QUALITY {
                (line_feeds & line_feeds.wrapping_sub(1)) == 0
            } else {
                line_feeds == 0
            };
            let mut mask = if skip_dna {
                Bitmask::ALL_DNA
            } else {
                bitpack_dna::<CONFIG>(v)
            };
            mask.line_feeds = line_feeds;
            mask
        }
    }
}

#[inline(always)]
pub fn extract_dna_bitmask<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    unsafe {
        let v = _mm512_loadu_si512(buf.as_ptr() as *const __m512i);
        bitpack_dna::<CONFIG>(v)
    }
}

#[inline(always)]
fn bitpack_dna<const CONFIG: Config>(v: __m512i) -> Bitmask {
    unsafe {
        let mut mask = Bitmask::ALL_DNA;

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            mask.high_bit = _mm512_movepi8_mask(_mm512_slli_epi16(v, 5));
            mask.low_bit = _mm512_movepi8_mask(_mm512_slli_epi16(v, 6));
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
                mask.two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
            #[cfg(any(not(target_feature = "bmi2"), feature = "no-pdep"))]
            {
                let iv = _mm512_permutex_epi64::<0xD8>(v);
                let hi = _mm512_slli_epi16(iv, 5);
                let lo = _mm512_slli_epi16(iv, 6);
                let mm_lo = _mm512_movepi8_mask(_mm512_unpacklo_epi8(lo, hi));
                let mm_hi = _mm512_movepi8_mask(_mm512_unpackhi_epi8(lo, hi));
                let mm_1 = ((mm_hi as u32 as u64) << 32) | (mm_lo as u32 as u64);
                let mm_2 = (mm_hi & 0xFFFFFFFF_00000000) | (mm_lo >> 32);
                mask.two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
            }
        }

        if flag_is_set(
            CONFIG,
            SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N,
        ) {
            let mask_two_bits = _mm512_set1_epi8(0b110i8);
            let mask_upper = _mm512_set1_epi8(0b11011111u8 as i8);
            let uv = _mm512_and_si512(v, mask_upper);

            let mask_actg = _mm512_cmpeq_epi8_mask(
                _mm512_shuffle_epi8(LUT_ACTG, _mm512_and_si512(v, mask_two_bits)),
                uv,
            );
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                mask.is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask.mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask.mask_n = u8_mask(uv, ASCII_N);
            }
        }

        mask
    }
}

#[inline(always)]
pub fn u8_mask(v: __m512i, vc: __m512i) -> u64 {
    unsafe { _mm512_cmpeq_epi8_mask(v, vc) }
}
