#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use fearless_simd::prelude::*;
#[allow(unused_imports)]
use fearless_simd::{Level, Simd, u8x64};

const NEEDS_DNA_MASKS: Config = COMPUTE_DNA_COLUMNAR
    | COMPUTE_DNA_PACKED
    | SPLIT_NON_ACTG
    | COMPUTE_MASK_NON_ACTG
    | COMPUTE_MASK_N;

#[inline(always)]
#[allow(dead_code, clippy::missing_const_for_fn)]
fn interleave_32_to_64(lo: u32, hi: u32) -> u64 {
    #[cfg(all(
        target_arch = "x86_64",
        target_feature = "bmi2",
        not(feature = "no-pdep")
    ))]
    unsafe {
        core::arch::x86_64::_pdep_u64(hi as u64, 0xAAAAAAAAAAAAAAAA)
            | core::arch::x86_64::_pdep_u64(lo as u64, 0x5555555555555555)
    }
    #[cfg(not(all(
        target_arch = "x86_64",
        target_feature = "bmi2",
        not(feature = "no-pdep")
    )))]
    {
        let lo_spread = {
            let mut x = lo as u64;
            x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
            x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
            x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
            x = (x | (x << 2)) & 0x3333333333333333;
            x = (x | (x << 1)) & 0x5555555555555555;
            x
        };
        let hi_spread = {
            let mut x = hi as u64;
            x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
            x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
            x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
            x = (x | (x << 2)) & 0x3333333333333333;
            x = (x | (x << 1)) & 0x5555555555555555;
            x
        };
        lo_spread | (hi_spread << 1)
    }
}

#[inline(always)]
#[allow(dead_code)]
fn interleave_64_to_128(low_bit: u64, high_bit: u64) -> u128 {
    let w0 = interleave_32_to_64(low_bit as u32, high_bit as u32) as u128;
    let w1 = interleave_32_to_64((low_bit >> 32) as u32, (high_bit >> 32) as u32) as u128;
    w0 | (w1 << 64)
}

#[inline(always)]
pub fn extract_bitmask<const CONFIG: Config, const MODE: u8>(buf: &[u8]) -> Bitmask {
    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
    {
        extract_bitmask_impl::<_, CONFIG, MODE>(
            unsafe { fearless_simd::Neon::assume_supported() },
            buf,
        )
    }
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        extract_bitmask_impl::<_, CONFIG, MODE>(
            unsafe { fearless_simd::Avx2::assume_supported() },
            buf,
        )
    }
    #[cfg(not(any(
        all(target_arch = "aarch64", target_feature = "neon"),
        all(target_arch = "x86_64", target_feature = "avx2")
    )))]
    {
        fearless_simd::dispatch!(Level::new(), simd => {
            extract_bitmask_impl::<_, CONFIG, MODE>(simd, buf)
        })
    }
}

#[inline(always)]
pub fn extract_dna_bitmask<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
    {
        extract_dna_bitmask_impl::<_, CONFIG>(
            unsafe { fearless_simd::Neon::assume_supported() },
            buf,
        )
    }
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        extract_dna_bitmask_impl::<_, CONFIG>(
            unsafe { fearless_simd::Avx2::assume_supported() },
            buf,
        )
    }
    #[cfg(not(any(
        all(target_arch = "aarch64", target_feature = "neon"),
        all(target_arch = "x86_64", target_feature = "avx2")
    )))]
    {
        fearless_simd::dispatch!(Level::new(), simd => {
            extract_dna_bitmask_impl::<_, CONFIG>(simd, buf)
        })
    }
}

#[inline(always)]
fn extract_bitmask_impl<S: Simd, const CONFIG: Config, const MODE: u8>(
    simd: S,
    buf: &[u8],
) -> Bitmask {
    let mut tmp = [0u8; 64];
    let slice_64 = if buf.len() >= 64 {
        unsafe { core::slice::from_raw_parts(buf.as_ptr(), 64) }
    } else {
        tmp[..buf.len()].copy_from_slice(buf);
        &tmp[..]
    };
    let v = u8x64::from_slice(simd, slice_64);

    if MODE == SCAN_MODE_DNA || flag_is_not_set(CONFIG, NEEDS_DNA_MASKS) {
        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
            bitpack_dna_impl::<S, CONFIG>(simd, v);

        let line_feeds = if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            if is_dna == !0 {
                0
            } else {
                v.simd_eq(u8x64::splat(simd, b'\n')).to_bitmask()
            }
        } else {
            v.simd_eq(u8x64::splat(simd, b'\n')).to_bitmask()
        };

        Bitmask {
            line_feeds,
            is_dna,
            two_bits,
            high_bit,
            low_bit,
            mask_non_actg,
            mask_n,
        }
    } else {
        let line_feeds = v.simd_eq(u8x64::splat(simd, b'\n')).to_bitmask();
        let skip_dna = if MODE == SCAN_MODE_QUALITY {
            (line_feeds & line_feeds.wrapping_sub(1)) == 0
        } else {
            line_feeds == 0
        };
        let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) = if skip_dna {
            (!0, 0, 0, 0, 0, 0)
        } else {
            bitpack_dna_impl::<S, CONFIG>(simd, v)
        };

        Bitmask {
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
fn extract_dna_bitmask_impl<S: Simd, const CONFIG: Config>(simd: S, buf: &[u8]) -> Bitmask {
    let mut tmp = [0u8; 64];
    let slice_64 = if buf.len() >= 64 {
        unsafe { core::slice::from_raw_parts(buf.as_ptr(), 64) }
    } else {
        tmp[..buf.len()].copy_from_slice(buf);
        &tmp[..]
    };
    let v = u8x64::from_slice(simd, slice_64);
    let (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n) =
        bitpack_dna_impl::<S, CONFIG>(simd, v);

    Bitmask {
        line_feeds: 0,
        is_dna,
        two_bits,
        high_bit,
        low_bit,
        mask_non_actg,
        mask_n,
    }
}

#[inline(always)]
fn bitpack_dna_impl<S: Simd, const CONFIG: Config>(
    simd: S,
    v: u8x64<S>,
) -> (u64, u128, u64, u64, u64, u64) {
    let mut is_dna = !0;
    let mut two_bits = 0;
    let mut high_bit = 0;
    let mut low_bit = 0;
    let mut mask_non_actg = 0;
    let mut mask_n = 0;

    if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
        let hi_mask = (v & u8x64::splat(simd, 0b100)).simd_ne(u8x64::splat(simd, 0));
        let lo_mask = (v & u8x64::splat(simd, 0b010)).simd_ne(u8x64::splat(simd, 0));
        high_bit = hi_mask.to_bitmask();
        low_bit = lo_mask.to_bitmask();

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            {
                use core::arch::aarch64::*;
                unsafe {
                    let raw: [u8; 64] = v.into();
                    let native = vld4q_u8(raw.as_ptr());
                    let s0 = vshlq_n_u8::<5>(native.0);
                    let s1 = vshlq_n_u8::<5>(native.1);
                    let s2 = vshlq_n_u8::<5>(native.2);
                    let s3 = vshlq_n_u8::<5>(native.3);
                    let bitpacked = vsriq_n_u8(vsriq_n_u8(s3, s2, 2), vsriq_n_u8(s1, s0, 2), 4);
                    two_bits = core::mem::transmute(bitpacked);
                }
            }
            #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
            {
                two_bits = interleave_64_to_128(low_bit, high_bit);
            }
        }
    }

    if flag_is_set(
        CONFIG,
        SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N,
    ) {
        let uppercase = v & u8x64::splat(simd, 0b11011111);
        let indices = v & u8x64::splat(simd, 0b110);
        let lut = u8x64::from_slice(
            simd,
            b"A_C_T_G_________A_C_T_G_________A_C_T_G_________A_C_T_G_________",
        );
        let lookup = lut.swizzle_dyn_within_blocks(indices);
        let eq = uppercase.simd_eq(lookup);

        let mask_actg = eq.to_bitmask();
        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            is_dna = mask_actg;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            mask_non_actg = !mask_actg;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            mask_n = uppercase.simd_eq(u8x64::splat(simd, b'N')).to_bitmask();
        }
    }

    (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n)
}
