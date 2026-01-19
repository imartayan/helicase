#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::x86_64::*;
use core::mem::transmute;

const GREATER_THAN: __m256i = unsafe { transmute([b'>'; 32]) };
const LINE_FEED: __m256i = unsafe { transmute([b'\n'; 32]) };
const LUT_ACTG: __m256i = unsafe { transmute(*b"A_C_T_G_________A_C_T_G_________") };

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m256i;
        let v_buf1 = _mm256_loadu_si256(ptr);
        let v_buf2 = _mm256_loadu_si256(ptr.add(1));

        let open_bracket = u8_mask(v_buf1, v_buf2, GREATER_THAN);
        let line_feeds = u8_mask(v_buf1, v_buf2, LINE_FEED);

        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;

        let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2) =
            if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
                (
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf1, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf1, 6)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf2, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf2, 6)) as u32 as u64,
                )
            } else {
                (0, 0, 0, 0)
            };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit = mm_hi_1 | (mm_hi_2 << 32);
            low_bit = mm_lo_1 | (mm_lo_2 << 32);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            let mm_1 =
                _pdep_u64(mm_hi_1, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_1, 0x5555555555555555);
            let mm_2 =
                _pdep_u64(mm_hi_2, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_2, 0x5555555555555555);
            two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            let mask_two_bits = _mm256_set1_epi8(0b110i8);
            let mask_upper = _mm256_set1_epi8(0b11011111u8 as i8);

            is_dna = movemask_64(
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v_buf1, mask_two_bits)),
                    _mm256_and_si256(v_buf1, mask_upper),
                ),
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v_buf2, mask_two_bits)),
                    _mm256_and_si256(v_buf2, mask_upper),
                ),
            );
        }

        FastaBitmask {
            open_bracket,
            line_feeds,
            is_dna,
            two_bits,
            high_bit,
            low_bit,
        }
    }
}

#[inline(always)]
pub fn extract_fastq_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastqBitmask {
    unsafe {
        let ptr = buf.as_ptr() as *const __m256i;
        let v_buf1 = _mm256_loadu_si256(ptr);
        let v_buf2 = _mm256_loadu_si256(ptr.add(1));

        let line_feeds = u8_mask(v_buf1, v_buf2, LINE_FEED);

        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;

        let (mm_hi_1, mm_lo_1, mm_hi_2, mm_lo_2) =
            if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
                (
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf1, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf1, 6)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf2, 5)) as u32 as u64,
                    _mm256_movemask_epi8(_mm256_slli_epi16(v_buf2, 6)) as u32 as u64,
                )
            } else {
                (0, 0, 0, 0)
            };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit = mm_hi_1 | (mm_hi_2 << 32);
            low_bit = mm_lo_1 | (mm_lo_2 << 32);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            let mm_1 =
                _pdep_u64(mm_hi_1, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_1, 0x5555555555555555);
            let mm_2 =
                _pdep_u64(mm_hi_2, 0xAAAAAAAAAAAAAAAA) | _pdep_u64(mm_lo_2, 0x5555555555555555);
            two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            let mask_two_bits = _mm256_set1_epi8(0b110i8);
            let mask_upper = _mm256_set1_epi8(0b11011111u8 as i8);

            is_dna = movemask_64(
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v_buf1, mask_two_bits)),
                    _mm256_and_si256(v_buf1, mask_upper),
                ),
                _mm256_cmpeq_epi8(
                    _mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v_buf2, mask_two_bits)),
                    _mm256_and_si256(v_buf2, mask_upper),
                ),
            );
        }

        FastqBitmask {
            line_feeds,
            is_dna,
            two_bits,
            high_bit,
            low_bit,
        }
    }
}

#[inline(always)]
fn movemask_64(v1: __m256i, v2: __m256i) -> u64 {
    unsafe {
        (_mm256_movemask_epi8(v1) as u32 as u64) | ((_mm256_movemask_epi8(v2) as u32 as u64) << 32)
    }
}

#[inline(always)]
pub fn u8_mask(v_buf: __m256i, v_buf2: __m256i, v_c: __m256i) -> u64 {
    unsafe {
        let cmp_c = _mm256_cmpeq_epi8(v_buf, v_c);
        let cmp_c2 = _mm256_cmpeq_epi8(v_buf2, v_c);
        let a = _mm256_movemask_epi8(cmp_c) as u32 as u64;
        let b = _mm256_movemask_epi8(cmp_c2) as u32 as u64;
        a | (b << 32)
    }
}
