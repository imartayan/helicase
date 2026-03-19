//! A vectorized library for FASTA/FASTQ parsing and bitpacking.
//!
//! # Requirements
//!
//! This library requires AVX2, SSE3, or NEON instruction sets. Enable `target-cpu=native` when
//! building:
//!
//! ```sh
//! RUSTFLAGS="-C target-cpu=native" cargo run --release
//! ```
//!
//! If your CPU has poor support for the
//! [PDEP instruction](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract)
//! (e.g. AMD CPUs prior to 2020), use the `no-pdep` feature:
//!
//! ```sh
//! RUSTFLAGS="-C target-cpu=native" cargo run --release -F no-pdep
//! ```
//!
//! # Minimal example
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! // set the options of the parser (at compile-time)
//! const CONFIG: Config = ParserOptions::default().config();
//!
//! fn main() {
//!     let path = "...";
//!
//!     // create a parser with the desired options
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
//!
//!     // iterate over records
//!     while let Some(_event) = parser.next() {
//!         // get a reference to the header
//!         let header = parser.get_header();
//!
//!         // get a reference to the sequence (without newlines)
//!         let seq = parser.get_dna_string();
//!
//!         // ...
//!     }
//! }
//! ```
//!
//! # Adjusting the configuration
//!
//! The parser is configured at compile-time via [`ParserOptions`]. For example, to skip headers
//! and ignore non-ACTG bases:
//!
//! ```rust
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default()
//!     .ignore_headers()
//!     .skip_non_actg()
//!     .config();
//! ```
//!
//! # Bitpacked DNA formats
//!
//! The parser can output a bitpacked representation of the sequence in two formats:
//!
//! - [`dna_format::PackedDNA`] — maps each base to two bits and packs them (compatible with
//!   [packed-seq](https://github.com/rust-seq/packed-seq) via the `packed-seq` feature).
//! - [`dna_format::ColumnarDNA`] — separates the high bit and the low bit of each base into two
//!   bitmasks.
//!
//! Since each base is encoded in two bits, non-ACTG bases must be handled explicitly. Three
//! options are available via [`ParserOptions`]:
//!
//! - [`split_non_actg`](ParserOptions::split_non_actg) — splits the sequence at non-ACTG bases,
//!   yielding multiple [`parser::Event::DnaChunk`] events.
//! - [`skip_non_actg`](ParserOptions::skip_non_actg) — skips non-ACTG bases and merges the
//!   remaining chunks, yielding one [`parser::Event::Record`] event per record.
//! - [`keep_non_actg`](ParserOptions::keep_non_actg) — keeps non-ACTG bases with a lossy
//!   two-bit encoding, yielding one [`parser::Event::Record`] event per record.
//!
//! # Iterating over chunks of packed DNA
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default()
//!     .dna_packed()
//!     // don't stop the iterator at the end of a record
//!     .return_record(false)
//!     .config();
//!
//! fn main() {
//!     let path = "...";
//!
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
//!
//!     // iterate over each chunk of ACTG bases
//!     while let Some(_event) = parser.next() {
//!         // headers are still accessible between chunks
//!         let header = parser.get_header();
//!
//!         // get a reference to the packed sequence
//!         let seq = parser.get_dna_packed();
//!
//!         // or directly get a PackedSeq (requires the packed-seq feature)
//!         // let packed_seq = parser.get_packed_seq();
//!     }
//! }
//! ```
//!
//! # Crate features
//!
//! | Feature    | Default | Description |
//! |------------|---------|-------------|
//! | `packed-seq` | no    | conversion to [packed-seq](https://github.com/rust-seq/packed-seq) types |
//! | `no-pdep`  | no      | disable PDEP instruction (recommended for AMD CPUs prior to 2020) |
//! | `gz`       | yes     | gzip decompression via [deko](https://github.com/igankevich/deko) |
//! | `zstd`     | yes     | zstd decompression |
//! | `bz2`      | no      | bzip2 decompression |
//! | `xz`       | no      | xz decompression |

pub(crate) mod carrying_add;
pub mod config;
pub mod dna_format;
pub mod input;
pub(crate) mod lexer;
pub mod parser;

pub use config::{Config, ParserOptions};
pub use parser::{FastaParser, FastqParser, FastxParser, HelicaseParser};

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
