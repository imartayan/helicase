#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::aarch64::*;
use core::mem::transmute;

const LINE_FEED: uint8x16_t = unsafe { transmute([b'\n'; 16]) };
const ASCII_N: uint8x16_t = unsafe { transmute([b'N'; 16]) };
const UPPERCASE: uint8x16_t = unsafe { transmute([0b11011111u8; 16]) };

// 0xFF at `c ^ 0x40` for c in ACGT, indices >= 32 return 0.
const LUT_ACTG_XOR: uint8x16x2_t = {
    let mut lut = [0u8; 32];
    lut[(b'A' ^ 0x40) as usize] = 0xFF;
    lut[(b'C' ^ 0x40) as usize] = 0xFF;
    lut[(b'G' ^ 0x40) as usize] = 0xFF;
    lut[(b'T' ^ 0x40) as usize] = 0xFF;
    unsafe { transmute(lut) }
};

const NEEDS_DNA_MASKS: Config = COMPUTE_DNA_COLUMNAR
    | COMPUTE_DNA_PACKED
    | SPLIT_NON_ACTG
    | COMPUTE_MASK_NON_ACTG
    | COMPUTE_MASK_N;

#[inline(always)]
pub fn extract_bitmask<const CONFIG: Config, const MODE: u8>(buf: &[u8]) -> Bitmask {
    unsafe {
        let ptr = buf.as_ptr();

        if flag_is_not_set(CONFIG, NEEDS_DNA_MASKS) {
            return Bitmask {
                line_feeds: extract_bitmask_contiguous(ptr),
                ..Bitmask::ALL_DNA
            };
        }

        if MODE == SCAN_MODE_DNA {
            let v = vld4q_u8(ptr);
            let mut mask = bitpack_dna::<CONFIG>(v);
            if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) || mask.is_dna != !0 {
                mask.line_feeds = movemask_64(map_8x16x4(v, |v| vceqq_u8(v, LINE_FEED)));
            }
            mask
        } else {
            let line_feeds = extract_bitmask_contiguous(ptr);
            let skip_dna = if MODE == SCAN_MODE_QUALITY {
                (line_feeds & line_feeds.wrapping_sub(1)) == 0
            } else {
                line_feeds == 0
            };
            let mut mask = if skip_dna {
                Bitmask::ALL_DNA
            } else {
                bitpack_dna::<CONFIG>(vld4q_u8(ptr))
            };
            mask.line_feeds = line_feeds;
            mask
        }
    }
}

// Lower latency, line feed only
#[inline(always)]
unsafe fn extract_bitmask_contiguous(ptr: *const u8) -> u64 {
    unsafe {
        let q0 = vld1q_u8(ptr);
        let q1 = vld1q_u8(ptr.add(16));
        let q2 = vld1q_u8(ptr.add(32));
        let q3 = vld1q_u8(ptr.add(48));

        let eq0 = vceqq_u8(q0, LINE_FEED);
        let eq1 = vceqq_u8(q1, LINE_FEED);
        let eq2 = vceqq_u8(q2, LINE_FEED);
        let eq3 = vceqq_u8(q3, LINE_FEED);

        movemask_64_contiguous(eq0, eq1, eq2, eq3)
    }
}

#[inline(always)]
pub fn extract_dna_bitmask<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    unsafe { bitpack_dna::<CONFIG>(vld4q_u8(buf.as_ptr())) }
}

#[inline(always)]
fn bitpack_dna<const CONFIG: Config>(v: uint8x16x4_t) -> Bitmask {
    unsafe {
        let mut mask = Bitmask::ALL_DNA;

        let shift_5 = if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
            map_8x16x4(v, |v| vshlq_n_u8::<5>(v))
        } else {
            v
        };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            let shift_6 = map_8x16x4(v, |v| vshlq_n_u8::<6>(v));
            mask.high_bit = movemask_64(shift_5);
            mask.low_bit = movemask_64(shift_6);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            let bitpacked = vsriq_n_u8(
                vsriq_n_u8(shift_5.3, shift_5.2, 2),
                vsriq_n_u8(shift_5.1, shift_5.0, 2),
                4,
            );
            mask.two_bits = transmute(bitpacked);
        }

        if flag_is_set(
            CONFIG,
            SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N,
        ) {
            let eq = map_8x16x4(v, |v| actg_mask(v));
            if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                && flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N)
            {
                let all_and = vandq_u8(vandq_u8(eq.0, eq.1), vandq_u8(eq.2, eq.3));
                if vminvq_u8(all_and) != 0xFF {
                    mask.is_dna = movemask_64(eq);
                }
            } else {
                let mask_actg = movemask_64(eq);
                if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                    mask.is_dna = mask_actg;
                }
                if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                    mask.mask_non_actg = !mask_actg;
                }
                if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                    mask.mask_n =
                        movemask_64(map_8x16x4(v, |v| vceqq_u8(vandq_u8(v, UPPERCASE), ASCII_N)));
                }
            }
        }

        mask
    }
}

#[inline(always)]
pub fn all_actg(buf: &[u8]) -> bool {
    unsafe {
        let ptr = buf[..256].as_ptr();
        let is_actg = |i: usize| actg_mask(vld1q_u8(ptr.add(16 * i)));
        let mut acc = [is_actg(0), is_actg(1), is_actg(2), is_actg(3)];
        for i in (4..16).step_by(4) {
            for (j, a) in acc.iter_mut().enumerate() {
                *a = vandq_u8(*a, is_actg(i + j));
            }
        }
        let all = vandq_u8(vandq_u8(acc[0], acc[1]), vandq_u8(acc[2], acc[3]));
        vminvq_u8(all) == 0xFF
    }
}

#[inline(always)]
unsafe fn actg_mask(v: uint8x16_t) -> uint8x16_t {
    unsafe {
        vqtbl2q_u8(
            LUT_ACTG_XOR,
            veorq_u8(vandq_u8(v, UPPERCASE), vdupq_n_u8(0x40)),
        )
    }
}

#[inline(always)]
fn map_8x16x4<F>(v: uint8x16x4_t, mut f: F) -> uint8x16x4_t
where
    F: FnMut(uint8x16_t) -> uint8x16_t,
{
    uint8x16x4_t(f(v.0), f(v.1), f(v.2), f(v.3))
}

// computing movemask is significantly more expensive than on x86

/// Movemask for 4 de-interleaved 16-byte vector lanes (from `vld4q_u8`).
#[inline(always)]
fn movemask_64(v: uint8x16x4_t) -> u64 {
    // https://stackoverflow.com/questions/74722950/convert-vector-compare-mask-into-bit-mask-in-aarch64-simd-or-arm-neon/74748402#74748402
    unsafe {
        let acc = vsriq_n_u8(vsriq_n_u8(v.3, v.2, 1), vsriq_n_u8(v.1, v.0, 1), 2);
        vget_lane_u64(
            vreinterpret_u64_u8(vshrn_n_u16(
                vreinterpretq_u16_u8(vsriq_n_u8(acc, acc, 4)),
                4,
            )),
            0,
        )
    }
}

/// Movemask for 4 contiguous 16-byte vector lanes (from `vld1q_u8` / `ldp`).
#[inline(always)]
unsafe fn movemask_64_contiguous(
    eq0: uint8x16_t,
    eq1: uint8x16_t,
    eq2: uint8x16_t,
    eq3: uint8x16_t,
) -> u64 {
    unsafe {
        let lo = vsriq_n_u8::<1>(vuzp2q_u8(eq0, eq1), vuzp1q_u8(eq0, eq1));
        let hi = vsriq_n_u8::<1>(vuzp2q_u8(eq2, eq3), vuzp1q_u8(eq2, eq3));
        let packed = vsriq_n_u8::<2>(vuzp2q_u8(lo, hi), vuzp1q_u8(lo, hi));
        let lo_half = vget_low_u8(packed);
        let hi_half = vget_high_u8(packed);
        let final_packed = vsri_n_u8::<4>(vuzp2_u8(lo_half, hi_half), vuzp1_u8(lo_half, hi_half));
        vget_lane_u64::<0>(vreinterpret_u64_u8(final_packed))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GREATER_THAN: uint8x16_t = unsafe { transmute([b'>'; 16]) };

    #[test]
    fn movemask_test() {
        unsafe {
            let buf = *b">CGT>CGT>CG>CG>CGTACGTACGT>CGT>CGT>CGT>CGT>CGT>CGTACGTACGT>CGT>A";
            let v = vld4q_u8(buf.as_ptr());
            let mask = map_8x16x4(v, |v| vceqq_u8(v, GREATER_THAN));
            let res = movemask_64(mask);
            assert_eq!(
                0b0100010000000000010001000100010001000100000000000100100100010001,
                res
            );
        }
    }

    #[test]
    fn movemask_contiguous_test() {
        unsafe {
            let buf = *b">CGT>CGT>CG>CG>CGTACGTACGT>CGT>CGT>CGT>CGT>CGT>CGTACGTACGT>CGT>A";
            let eq0 = vceqq_u8(vld1q_u8(buf.as_ptr()), GREATER_THAN);
            let eq1 = vceqq_u8(vld1q_u8(buf.as_ptr().add(16)), GREATER_THAN);
            let eq2 = vceqq_u8(vld1q_u8(buf.as_ptr().add(32)), GREATER_THAN);
            let eq3 = vceqq_u8(vld1q_u8(buf.as_ptr().add(48)), GREATER_THAN);
            let res = movemask_64_contiguous(eq0, eq1, eq2, eq3);
            assert_eq!(
                0b0100010000000000010001000100010001000100000000000100100100010001,
                res
            );
        }
    }
}
