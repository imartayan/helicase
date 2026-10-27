#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::wasm32::*;
use core::mem::transmute;

const LINE_FEED: v128 = unsafe { transmute([b'\n'; 16]) };
const ASCII_N: v128 = unsafe { transmute([b'N'; 16]) };
const LUT_ACTG: v128 = unsafe { transmute(*b"A_C_T_G_________") };
const TWO_BITS: v128 = unsafe { transmute([0b110u8; 16]) };
const UPPERCASE: v128 = unsafe { transmute([0b11011111u8; 16]) };

type Vectors = (v128, v128, v128, v128);

#[inline(always)]
pub fn line_feeds((v1, v2, v3, v4): Vectors) -> u64 {
    u8_mask(v1, v2, v3, v4, LINE_FEED)
}

#[inline(always)]
pub fn scan_line_feeds(buf: &[u8]) -> u64 {
    line_feeds(load(buf))
}

#[inline(always)]
pub fn bitpack_dna<const CONFIG: Config>((v1, v2, v3, v4): Vectors) -> Bitmask {
    let mut mask = Bitmask::ALL_DNA;

    if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
        // Bit 2 (high) and bit 1 (low) of each byte moved to the sign bit.
        let (hi_1, lo_1, hi_2, lo_2, hi_3, lo_3, hi_4, lo_4) = (
            i16x8_shl(v1, 5),
            i16x8_shl(v1, 6),
            i16x8_shl(v2, 5),
            i16x8_shl(v2, 6),
            i16x8_shl(v3, 5),
            i16x8_shl(v3, 6),
            i16x8_shl(v4, 5),
            i16x8_shl(v4, 6),
        );

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            mask.high_bit = movemask_64(hi_1, hi_2, hi_3, hi_4);
            mask.low_bit = movemask_64(lo_1, lo_2, lo_3, lo_4);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            let mm_1 = movemask_64(
                interleave_lo(lo_1, hi_1),
                interleave_hi(lo_1, hi_1),
                interleave_lo(lo_2, hi_2),
                interleave_hi(lo_2, hi_2),
            );
            let mm_2 = movemask_64(
                interleave_lo(lo_3, hi_3),
                interleave_hi(lo_3, hi_3),
                interleave_lo(lo_4, hi_4),
                interleave_hi(lo_4, hi_4),
            );
            mask.two_bits = (mm_1 as u128) | ((mm_2 as u128) << 64);
        }
    }

    if flag_is_set(
        CONFIG,
        SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N,
    ) {
        let uv1 = v128_and(v1, UPPERCASE);
        let uv2 = v128_and(v2, UPPERCASE);
        let uv3 = v128_and(v3, UPPERCASE);
        let uv4 = v128_and(v4, UPPERCASE);

        let mask_actg = movemask_64(
            u8x16_eq(i8x16_swizzle(LUT_ACTG, v128_and(v1, TWO_BITS)), uv1),
            u8x16_eq(i8x16_swizzle(LUT_ACTG, v128_and(v2, TWO_BITS)), uv2),
            u8x16_eq(i8x16_swizzle(LUT_ACTG, v128_and(v3, TWO_BITS)), uv3),
            u8x16_eq(i8x16_swizzle(LUT_ACTG, v128_and(v4, TWO_BITS)), uv4),
        );
        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            mask.is_dna = mask_actg;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            mask.mask_non_actg = !mask_actg;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            mask.mask_n = u8_mask(uv1, uv2, uv3, uv4, ASCII_N);
        }
    }

    mask
}

#[inline(always)]
pub fn load(buf: &[u8]) -> Vectors {
    // Blocks are 64 bytes long, shorter blocks are padded by the input.
    unsafe {
        let ptr = buf.as_ptr() as *const v128;
        (
            v128_load(ptr),
            v128_load(ptr.add(1)),
            v128_load(ptr.add(2)),
            v128_load(ptr.add(3)),
        )
    }
}

#[inline(always)]
fn interleave_lo(a: v128, b: v128) -> v128 {
    i8x16_shuffle::<0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23>(a, b)
}

#[inline(always)]
fn interleave_hi(a: v128, b: v128) -> v128 {
    i8x16_shuffle::<8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31>(a, b)
}

#[inline(always)]
fn movemask_64(v1: v128, v2: v128, v3: v128, v4: v128) -> u64 {
    (u8x16_bitmask(v1) as u64)
        | ((u8x16_bitmask(v2) as u64) << 16)
        | ((u8x16_bitmask(v3) as u64) << 32)
        | ((u8x16_bitmask(v4) as u64) << 48)
}

#[inline(always)]
pub fn u8_mask(v1: v128, v2: v128, v3: v128, v4: v128, vc: v128) -> u64 {
    movemask_64(
        u8x16_eq(v1, vc),
        u8x16_eq(v2, vc),
        u8x16_eq(v3, vc),
        u8x16_eq(v4, vc),
    )
}
