//! SIMD kernels computing the bitmasks of 64-byte blocks.
//!
//! Each backend loads a block into vectors, and extracts its line feeds and DNA bitmasks.

use crate::config::{advanced::*, *};
use crate::lexer::*;

const NEEDS_DNA_MASKS: Config = COMPUTE_DNA_COLUMNAR
    | COMPUTE_DNA_PACKED
    | SPLIT_NON_ACTG
    | COMPUTE_MASK_NON_ACTG
    | COMPUTE_MASK_N;

#[cfg(all(target_feature = "avx512f", target_feature = "avx512bw"))]
#[path = "avx512.rs"]
mod backend;
#[cfg(all(
    not(all(target_feature = "avx512f", target_feature = "avx512bw")),
    target_feature = "avx2"
))]
#[path = "avx2.rs"]
mod backend;
#[cfg(all(
    not(all(target_feature = "avx512f", target_feature = "avx512bw")),
    not(target_feature = "avx2"),
    target_feature = "ssse3"
))]
#[deprecated(
    note = "Helicase currently uses SSE3 instead of AVX2/AVX-512 instructions. Compile using `-C target-cpu=native` to get better performances."
)]
#[path = "sse.rs"]
mod backend;
#[cfg(target_feature = "neon")]
#[path = "neon.rs"]
mod backend;
#[cfg(target_feature = "simd128")]
#[path = "simd128.rs"]
mod backend;
#[cfg(not(any(
    all(target_feature = "avx512f", target_feature = "avx512bw"),
    target_feature = "avx2",
    target_feature = "ssse3",
    target_feature = "neon",
    target_feature = "simd128"
)))]
#[deprecated(
    note = "Helicase currently uses (slow) non-vectorized instructions. Compile using `-C target-cpu=native` to get better performances."
)]
#[path = "fallback.rs"]
mod backend;

// The scalar fallback is a reference for the SIMD backends in tests.
#[cfg(all(
    test,
    any(
        target_feature = "neon",
        target_feature = "ssse3",
        target_feature = "simd128"
    )
))]
#[path = "fallback.rs"]
#[allow(dead_code)]
mod scalar;
#[cfg(all(
    test,
    any(
        target_feature = "neon",
        target_feature = "ssse3",
        target_feature = "simd128"
    )
))]
mod tests;

#[cfg(target_feature = "neon")]
pub use backend::all_actg;

#[inline(always)]
pub fn extract_bitmask<const CONFIG: Config, const MODE: u8>(buf: &[u8]) -> Bitmask {
    if flag_is_not_set(CONFIG, NEEDS_DNA_MASKS) {
        return Bitmask {
            line_feeds: backend::scan_line_feeds(buf),
            ..Bitmask::ALL_DNA
        };
    }
    if MODE == SCAN_MODE_DNA {
        let v = backend::load(buf);
        let mut mask = backend::bitpack_dna::<CONFIG>(v);
        if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) || mask.is_dna != !0 {
            mask.line_feeds = backend::line_feeds(v);
        }
        return mask;
    }
    let line_feeds = backend::scan_line_feeds(buf);
    let skip_dna = if MODE == SCAN_MODE_QUALITY {
        (line_feeds & line_feeds.wrapping_sub(1)) == 0
    } else {
        line_feeds == 0
    };
    let mut mask = if skip_dna {
        Bitmask::ALL_DNA
    } else {
        backend::bitpack_dna::<CONFIG>(backend::load(buf))
    };
    mask.line_feeds = line_feeds;
    mask
}

#[inline(always)]
pub fn extract_dna_bitmask<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    backend::bitpack_dna::<CONFIG>(backend::load(buf))
}
