mod carrying_add;
pub mod config;
pub mod dna_format;
pub mod input;
pub mod lexer;
pub mod parser;

pub use config::{Config, ParserOptions};
pub use parser::{Event, FastaParser, FastqParser, FastxParser, Parser};

#[cfg(target_feature = "avx2")]
pub(crate) mod simd {
    mod avx2;
    pub use avx2::*;
}
#[cfg(all(not(target_feature = "avx2"), target_feature = "ssse3"))]
#[deprecated(
    note = "Helicase currently uses SSE3 instead of AVX2 instructions. Compile using `-C target-cpu=native` to get better performances."
)]
pub(crate) mod simd {
    mod sse;
    pub use sse::*;
}
#[cfg(target_feature = "neon")]
pub(crate) mod simd {
    mod neon;
    pub use neon::*;
}
#[cfg(not(any(
    target_feature = "avx2",
    target_feature = "sse3",
    target_feature = "neon"
)))]
#[deprecated(
    note = "Helicase currently uses (slow) non-vectorized instructions. Compile using `-C target-cpu=native` to get better performances."
)]
pub(crate) mod simd {
    mod fallback;
    pub use fallback::*;
}
