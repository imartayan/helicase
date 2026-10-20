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
//! and build a [`FastxParser`]/[`ParallelFastxParser`] with this configuration.
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
//! # Input
//!
//! Files are opened with [`from_file`](input::FromFile::from_file) (which supports [transparent decompression](#crate-features)) or [`from_file_mmap`](input::FromMmap::from_file_mmap) (which only supports uncompressed data but can be faster).
//! In-memory slices ([`from_slice`](input::FromSlice::from_slice)) and generic streams ([`from_reader`](input::FromReader::from_reader)) are also supported.
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
//! An event signals a record boundary, a contiguous DNA chunk, or an I/O error:
//! - [`Event::Record(range)`](parser::Event::Record): emitted once per record boundary, carrying the byte offset `Range<usize>` of the entire record. Enabled by [`return_record`](ParserOptions::return_record) (on by default for string and lossy formats).
//! - [`Event::DnaChunk(range)`](parser::Event::DnaChunk): emitted for each contiguous ACTG run when non-ACTG characters are split, carrying the byte offset `Range<usize>` of the chunk. Enabled by [`return_dna_chunk`](ParserOptions::return_dna_chunk) (on by default with [`split_non_actg`](ParserOptions::split_non_actg), [`dna_packed`](ParserOptions::dna_packed), and [`dna_columnar`](ParserOptions::dna_columnar)).
//! - [`Event::IoError(err)`](parser::Event::IoError): emitted if an I/O or stream decompression error occurs during iteration.
//!
//! Sequence, header, and quality data are read directly from the parser methods ([`get_header`](HelicaseParser::get_header), [`get_dna_string`](HelicaseParser::get_dna_string), [`get_dna_packed`](HelicaseParser::get_dna_packed), [`get_quality`](HelicaseParser::get_quality), etc.).
//!
//! When both records and chunks are active (e.g. `dna_packed().return_record(true)`), match on the event to distinguish them:
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::parser::Event;
//! use helicase::*;
//!
//! // dna_packed enables DnaChunk events; explicitly turn Record events back on.
//! const CONFIG: Config = ParserOptions::default().dna_packed().return_record(true).config();
//!
//! fn main() -> std::io::Result<()> {
//!     let mut parser = FastxParser::<CONFIG>::from_file("reads.fastq.gz")?;
//!
//!     while let Some(event) = parser.next() {
//!         match event {
//!             Event::DnaChunk(range) => {
//!                 // One contiguous ACTG run is ready
//!                 let seq = parser.get_dna_packed();
//!             }
//!             Event::Record(range) => {
//!                 // Record boundary reached (all chunks of this record processed)
//!                 let header = parser.get_header();
//!             }
//!             Event::IoError(err) => {
//!                 return Err(err);
//!             }
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! When only one type of event is active (such as the default configuration yielding `Event::Record`), the event value can be safely ignored:
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! // Default config: only Record events, one per record.
//! const CONFIG: Config = ParserOptions::default().config();
//!
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path)?;
//!
//!     while let Some(_event) = parser.next() {
//!         let header = parser.get_header();
//!         let seq = parser.get_dna_string();
//!     }
//!     Ok(())
//! }
//! ```
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
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path)?;
//!
//!     // iterate over each chunk of ACTG bases
//!     while let Some(_event) = parser.next() {
//!         // headers are still accessible between chunks
//!         let header = parser.get_header();
//!         let seq = parser.get_dna_packed();
//!
//!         // or directly get a PackedSeq (requires the packed-seq feature)
//!         // let packed_seq = parser.get_packed_seq();
//!     }
//!     Ok(())
//! }
//! ```
//!
//! # Parallel parsing
//!
//! For multi-threaded processing, [`ParallelFastxParser`] can parse files in parallel across worker threads.
//! Like [`FastxParser`], it can be created using [`from_file`](ParallelFastxParser::from_file) or [`from_file_mmap`](ParallelFastxParser::from_file_mmap) (see [Input](#input)).
//!
//! You drive the parsing loop directly on each thread with [`for_each`](ParallelFastxParser::for_each) or [`reduce`](ParallelFastxParser::reduce):
//!
//! ```rust,no_run
//! use helicase::config::*;
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default().ignore_headers().config();
//!
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!
//!     // Process records in parallel across 8 threads using for_each
//!     ParallelFastxParser::<CONFIG>::from_file(&path)?
//!         .threads(8)
//!         .for_each(|mut parser| {
//!             while let Some(_event) = parser.next() {
//!                 let header = parser.get_header();
//!                 let seq = parser.get_dna_string();
//!                 // ...
//!             }
//!         })?;
//!
//!     // Sum total sequence length across 8 threads using reduce
//!     let total_len: usize = ParallelFastxParser::<CONFIG>::from_file(&path)?
//!         .threads(8)
//!         .reduce(
//!             |mut parser| {
//!                 let mut local_len = 0;
//!                 while let Some(_event) = parser.next() {
//!                     local_len += parser.get_dna_string().len();
//!                 }
//!                 local_len
//!             },
//!             |acc, val| acc + val,
//!         )?;
//!
//!     Ok(())
//! }
//! ```
//!
//! # Crate features
//!
//! ### Packed-seq
//!
//! The [PackedDNA format](#bitpacked-dna-formats) is compatible with [packed-seq](https://github.com/rust-seq/packed-seq) and can be converted when the `packed-seq` feature is enabled (disabled by default).
//!
//! This can be useful for [hashing *k*-mers](https://github.com/rust-seq/seq-hash) or [computing minimizers & syncmers](https://github.com/rust-seq/simd-minimizers).
//!
//! ### Decompression
//!
//! This library supports transparent file decompression using [deko](https://github.com/igankevich/deko), you can choose the supported formats using the following features:
//! - `gz` for gzip (default)
//! - `zstd` for zstd (default)
//! - `bz2` for bzip2
//! - `xz` for xz
//!
//! ### No PDEP
//!
//! By default, this library uses [PDEP](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) to compute the [PackedDNA format](#bitpacked-dna-formats).
//! However, this instruction can be slow on some older CPUs (especially AMD CPUs prior to 2020).
//! If you want an efficient implementation for these CPUs, we recommend using the `no-pdep` feature.

extern crate alloc;

pub mod config;
pub mod dna_format;
pub mod input;
pub(crate) mod lexer;
pub mod parser;

pub use config::{Config, ParserOptions};
pub use parser::{FastaParser, FastqParser, FastxParser, HelicaseParser, ParallelFastxParser};

#[cfg(feature = "needletail")]
pub mod needletail_helpers;
#[cfg(feature = "needletail")]
pub use needletail_helpers::*;

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
        target_feature = "ssse3",
        target_feature = "neon"
)))]
#[deprecated(
    note = "Helicase currently uses (slow) non-vectorized instructions. Compile using `-C target-cpu=native` to get better performances."
)]
pub(crate) mod simd {
    mod fallback;
    pub use fallback::*;
}
