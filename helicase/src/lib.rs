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
//! The main entry point is to define a configuration via [`ParserOptions`]
//! and build a [`FastxParser`] with this configuration.
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
//! The parser is configured at compile-time via [`ParserOptions`].
//! For example, to ignore headers and split non-ACTG bases:
//!
//! ```rust
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default()
//!     .ignore_headers()
//!     .split_non_actg()
//!     .config();
//! ```
//!
//! # Bitpacked DNA formats
//!
//! The parser can output a bitpacked representation of the sequence in two formats:
//! - [`PackedDNA`](dna_format::PackedDNA) maps each base to two bits and packs them
//!   (compatible with [packed-seq](https://github.com/rust-seq/packed-seq) via the `packed-seq` feature).
//! - [`ColumnarDNA`](dna_format::ColumnarDNA) separates the high bit and the low bit of each base into two bitmasks.
//!
//! Since each base is encoded in two bits, non-ACTG bases must be handled explicitly. Three
//! options are available via [`ParserOptions`]:
//! - [`split_non_actg`](ParserOptions::split_non_actg) splits the sequence at non-ACTG bases,
//!   yielding one [`DnaChunk`](parser::Event::DnaChunk) event per contiguous ACTG run (default for bitpacked formats),
//!   and disables `Record` events by default to avoid processing each sequence twice.
//! - [`skip_non_actg`](ParserOptions::skip_non_actg) skips non-ACTG bases and merges the remaining chunks,
//!   yielding one [`Record`](parser::Event::Record) event per record.
//! - [`keep_non_actg`](ParserOptions::keep_non_actg) keeps non-ACTG bases and encodes them lossily,
//!   yielding one [`Record`](parser::Event::Record) event per record (default for string format).
//!
//! # Events
//!
//! The parser is an iterator that yields [`Event`](parser::Event) values.
//! An event signals a record boundary or a contiguous DNA chunk,
//! but the data is always read from the parser itself via [`get_header`](HelicaseParser::get_header), [`get_dna_string`](HelicaseParser::get_dna_string), etc.
//!
//! There are two kinds of event:
//! - [`Event::Record`](parser::Event::Record) emitted once per record, after all of its DNA chunks.
//!   Enabled by [`return_record`](ParserOptions::return_record) (on by default, off by default
//!   whenever non-ACTG bases are split, to avoid processing each sequence twice).
//! - [`Event::DnaChunk`](parser::Event::DnaChunk) emitted for each contiguous ACTG run.
//!   Enabled by [`return_dna_chunk`](ParserOptions::return_dna_chunk) (on by default with
//!   [`dna_packed`](ParserOptions::dna_packed) and [`dna_columnar`](ParserOptions::dna_columnar)).
//!
//! When both are active you need to match on the event to distinguish them:
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::parser::Event;
//! use helicase::*;
//!
//! // dna_packed enables DnaChunk events; explicitly turn Record events back on.
//! const CONFIG: Config = ParserOptions::default().dna_packed().return_record(true).config();
//!
//! fn main() {
//!     let path = "...";
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
//!
//!     while let Some(event) = parser.next() {
//!         match event {
//!             Event::DnaChunk(_) => {
//!                 // one contiguous ACTG run is ready
//!                 let seq = parser.get_dna_packed();
//!             }
//!             Event::Record(_) => {
//!                 // all chunks of this record have been processed
//!             }
//!         }
//!     }
//! }
//! ```
//!
//! When only one type of event is active, the event value can be safely ignored:
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! // Default config: only Record events, one per record.
//! const CONFIG: Config = ParserOptions::default().config();
//!
//! fn main() {
//!     let path = "...";
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");
//!
//!     while let Some(_event) = parser.next() {
//!         let header = parser.get_header();
//!         let seq = parser.get_dna_string();
//!     }
//! }
//! ```
//!
//! It is even possible to disable all events to process the entire file in one go, for instance if you simply want to count bases.
//!
//! # Iterating over chunks of packed DNA
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default()
//!     // by default, dna_packed splits non-ACTG bases, stops after each chunk, and doesn't stop at the end of a record
//!     .dna_packed()
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
//! | Feature      | Default | Description |
//! |--------------|---------|-------------|
//! | `packed-seq` | no      | conversion to [packed-seq](https://github.com/rust-seq/packed-seq) types |
//! | `no-pdep`    | no      | disable PDEP instruction (recommended for AMD CPUs prior to 2020) |
//! | `gz`         | yes     | gzip decompression |
//! | `zstd`       | yes     | zstd decompression |
//! | `bz2`        | no      | bzip2 decompression |
//! | `xz`         | no      | xz decompression |

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
