// Compares the SIMD backend with the scalar fallback on the same blocks.
use crate::config::{advanced::*, *};
use crate::lexer::*;
use crate::scalar;

const FULL: Config =
    COMPUTE_DNA_PACKED | COMPUTE_DNA_COLUMNAR | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N;
const FULL_SPLIT: Config = FULL | SPLIT_NON_ACTG;
const PACKED_SPLIT: Config = COMPUTE_DNA_PACKED | SPLIT_NON_ACTG;
const COLUMNAR_SPLIT: Config = COMPUTE_DNA_COLUMNAR | SPLIT_NON_ACTG;

fn blocks() -> Vec<[u8; 64]> {
    let alphabet = b"ACGTACGTACGTacgtNNRYKM\n\n\r>@+\0";
    let mut state = 42u64;
    let mut next = |n: usize| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize % n
    };
    (0..4000)
        .map(|k| {
            let mut block = [0; 64];
            for byte in &mut block {
                // A quarter of the blocks only contain ACGT.
                let n = if k % 4 == 0 { 4 } else { alphabet.len() };
                *byte = alphabet[next(n)];
            }
            block
        })
        .collect()
}

fn check<const CONFIG: Config>() {
    let set = |flag| flag_is_set(CONFIG, flag);
    for block in blocks() {
        let id = String::from_utf8_lossy(&block);
        let simd = crate::simd::extract_bitmask::<CONFIG, SCAN_MODE_DNA>(&block);
        let expected = scalar::extract_bitmask::<CONFIG, SCAN_MODE_DNA>(&block);
        assert_eq!(simd.line_feeds, expected.line_feeds, "line feeds {id:?}");
        let line_feeds = expected.line_feeds;
        for (simd, expected) in [
            (simd, expected),
            (
                crate::simd::extract_dna_bitmask::<CONFIG>(&block),
                scalar::extract_dna_bitmask::<CONFIG>(&block),
            ),
        ] {
            if set(SPLIT_NON_ACTG) {
                assert_eq!(simd.is_dna, expected.is_dna, "is_dna {id:?}");
            }
            if set(COMPUTE_DNA_PACKED) {
                assert_eq!(simd.two_bits, expected.two_bits, "two_bits {id:?}");
            }
            if set(COMPUTE_DNA_COLUMNAR) {
                assert_eq!(simd.high_bit, expected.high_bit, "high_bit {id:?}");
                assert_eq!(simd.low_bit, expected.low_bit, "low_bit {id:?}");
            }
            if set(COMPUTE_MASK_NON_ACTG) {
                assert_eq!(
                    simd.mask_non_actg, expected.mask_non_actg,
                    "mask_non_actg {id:?}"
                );
            }
            if set(COMPUTE_MASK_N) {
                assert_eq!(simd.mask_n, expected.mask_n, "mask_n {id:?}");
            }
        }
        for simd in [
            crate::simd::extract_bitmask::<CONFIG, SCAN_MODE_HEADER>(&block),
            crate::simd::extract_bitmask::<CONFIG, SCAN_MODE_QUALITY>(&block),
        ] {
            assert_eq!(simd.line_feeds, line_feeds, "line feeds {id:?}");
        }
    }
}

#[test]
fn test_simd_matches_scalar() {
    check::<0>();
    check::<FULL>();
    check::<FULL_SPLIT>();
    check::<PACKED_SPLIT>();
    check::<COLUMNAR_SPLIT>();
}
