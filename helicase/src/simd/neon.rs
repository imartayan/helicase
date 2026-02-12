#![allow(clippy::missing_transmute_annotations)]

use crate::config::{advanced::*, *};
use crate::lexer::*;
use core::arch::aarch64::*;
use core::mem::transmute;

const GREATER_THAN: uint8x16_t = unsafe { transmute([b'>'; 16]) };
const LINE_FEED: uint8x16_t = unsafe { transmute([b'\n'; 16]) };
const ASCII_N: uint8x16_t = unsafe { transmute([b'N'; 16]) };
const UPPERCASE: uint8x16_t = unsafe { transmute([0b11011111u8; 16]) };
const TWO_BITS: uint8x16_t = unsafe { transmute([0b110u8; 16]) };
const LUT_ACTG: uint8x16_t = unsafe { transmute(*b"A_C_T_G_________") };

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    unsafe {
        let ptr = buf.as_ptr();
        let v = vld4q_u8(ptr);

        let open_bracket = movemask_64(map_8x16x4(v, |v| vceqq_u8(v, GREATER_THAN)));
        let line_feeds = movemask_64(map_8x16x4(v, |v| vceqq_u8(v, LINE_FEED)));

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
        let ptr = buf.as_ptr();
        let v = vld4q_u8(ptr);

        let line_feeds = movemask_64(map_8x16x4(v, |v| vceqq_u8(v, LINE_FEED)));

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
fn bitpack_dna<const CONFIG: Config>(v: uint8x16x4_t) -> (u64, u128, u64, u64, u64, u64) {
    unsafe {
        let mut is_dna = !0;
        let mut two_bits = 0;
        let mut high_bit = 0;
        let mut low_bit = 0;
        let mut mask_non_actg = 0;
        let mut mask_n = 0;

        let shift_5 = if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED) {
            map_8x16x4(v, |v| vshlq_n_u8::<5>(v))
        } else {
            v
        };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            let shift_6 = map_8x16x4(v, |v| vshlq_n_u8::<6>(v));
            high_bit = movemask_64(shift_5);
            low_bit = movemask_64(shift_6);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            let bitpacked = vsriq_n_u8(
                vsriq_n_u8(shift_5.3, shift_5.2, 2),
                vsriq_n_u8(shift_5.1, shift_5.0, 2),
                4,
            );
            two_bits = transmute(bitpacked);
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG) {
            let lookup = map_8x16x4(v, |v| vqtbl1q_u8(LUT_ACTG, vandq_u8(v, TWO_BITS)));
            let uppercase = map_8x16x4(v, |v| vandq_u8(v, UPPERCASE));
            let mask_actg =
                movemask_64(map_two_8x16x4(lookup, uppercase, |v1, v2| vceqq_u8(v1, v2)));
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask_n = movemask_64(map_8x16x4(uppercase, |v| vceqq_u8(v, ASCII_N)));
            }
        }

        (is_dna, two_bits, high_bit, low_bit, mask_non_actg, mask_n)
    }
}

#[inline(always)]
fn map_8x16x4<F>(v: uint8x16x4_t, mut f: F) -> uint8x16x4_t
where
    F: FnMut(uint8x16_t) -> uint8x16_t,
{
    uint8x16x4_t(f(v.0), f(v.1), f(v.2), f(v.3))
}

#[inline(always)]
fn map_two_8x16x4<F>(v1: uint8x16x4_t, v2: uint8x16x4_t, mut f: F) -> uint8x16x4_t
where
    F: FnMut(uint8x16_t, uint8x16_t) -> uint8x16_t,
{
    uint8x16x4_t(f(v1.0, v2.0), f(v1.1, v2.1), f(v1.2, v2.2), f(v1.3, v2.3))
}

// computing movemask is significantly more expensive than on x86

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ld4_test() {
        unsafe {
            let buf: [u8; 256] = core::array::from_fn(|i| i as u8);
            let vbuf = vld4q_u8(buf.as_ptr());
            dbg!(vbuf.0);
        }
    }

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
}
