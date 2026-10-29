//! A carefully optimized, SIMD-accelerated FASTA/FASTQ parser for DNA sequences.
//!
//! It is designed for speed and flexibility: zero-overhead compile-time configuration, robust handling of non-ACTG bases, and on-the-fly 2-bit bitpacking.
//!
//! It is typically 1.1-4× faster than `needletail` on a single thread, and 3-6× faster than `paraseq` across 4 threads.
//!
//! The underlying algorithm is described in [Helicase: Vectorized parsing and bitpacking of genomic sequences](https://doi.org/10.64898/2026.03.19.712912), please cite it if you use this library.
//!
//! # Requirements
//!
//! This library requires AVX-512, AVX2, SSE3 or NEON instruction sets, make sure to enable `target-cpu=native` when using it:
//!
//! ```sh
//! RUSTFLAGS="-C target-cpu=native" cargo run --release
//! ```
//!
//! Note: if your CPU has a bad support for the [PDEP instruction](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) (e.g. AMD CPUs prior to 2020), it is recommended to use the `no-pdep` [feature](#crate-features).
//!
//! # Usage
//!
//! If you are migrating from `needletail` or `paraseq`, see the [migration guide](https://github.com/imartayan/helicase/blob/main/MIGRATION.md).
//!
//! ## Minimal example
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
//! # #[cfg(feature = "std")]
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!
//!     // create a parser with the desired options
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path)?;
//!
//!     // iterate over records
//!     while let Some(event) = parser.next() {
//!         event?;
//!         // get a reference to the header
//!         let header = parser.get_header();
//!
//!         // get a reference to the sequence (without newlines)
//!         let seq = parser.get_dna_string();
//!
//!         // ...
//!     }
//!     Ok(())
//! }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//! ## Input
//!
//! Files are opened with [`from_file`](input::FromFile::from_file) (which supports [transparent decompression](#crate-features)) or [`from_file_mmap`](input::FromMmap::from_file_mmap) (which only supports uncompressed data but can be faster).
//! In-memory slices ([`from_slice`](input::FromSlice::from_slice)) and generic streams ([`from_reader`](input::FromReader::from_reader)) are also supported.
//!
//! ## Adjusting the configuration
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
//! ## Bitpacked DNA formats
//!
//! The parser can output a bitpacked representation of the sequence in two formats:
//! - [`PackedDNA`](dna_format::PackedDNA) maps each base to two bits and packs them
//!   (compatible with [packed-seq](https://github.com/rust-seq/packed-seq) via the `packed-seq` [feature](#crate-features)).
//! - [`ColumnarDNA`](dna_format::ColumnarDNA) separates the high bit and the low bit of each base into two bitmasks.
//!
//! Since 2-bit encoding only represents four bases, non-ACTG characters are handled explicitly.
//! Four options are available via [`ParserOptions`]:
//! - [`split_non_actg`](ParserOptions::split_non_actg) splits the sequence at non-ACTG bases,
//!   yielding one [`DnaChunk`](parser::Event::DnaChunk) event per contiguous ACTG run (default for bitpacked formats),
//!   and disables `Record` events by default to avoid processing each sequence twice.
//! - [`skip_non_actg`](ParserOptions::skip_non_actg) skips non-ACTG bases and merges the remaining chunks,
//!   yielding one [`Record`](parser::Event::Record) event per record.
//! - [`keep_non_actg`](ParserOptions::keep_non_actg) keeps the non-ACTG bases and encodes them lossily,
//!   yielding one [`Record`](parser::Event::Record) event per record (default for string format).
//! - [`reject_non_actg`](ParserOptions::reject_non_actg) returns an [error](parser::ErrorKind::InvalidBase) at the first non-ACTG base,
//!   yielding one [`Record`](parser::Event::Record) event per record.
//!
//! ## Events
//!
//! The parser is an iterator that yields [`Result<Event, parser::Error>`](parser::Result) values.
//! An event signals a record boundary or a contiguous DNA chunk:
//! - [`Event::Record(range)`](parser::Event::Record): emitted once per record boundary. The `range` gives the byte offsets `start..end` of the entire record. Enabled by [`return_record`](ParserOptions::return_record) (on by default for string and lossy formats).
//! - [`Event::DnaChunk(range)`](parser::Event::DnaChunk): emitted for each contiguous run of ACTG bases when non-ACTG characters are split. The `range` gives the byte offsets `start..end` of the chunk. Enabled by [`return_dna_chunk`](ParserOptions::return_dna_chunk) (on by default with [`split_non_actg`](ParserOptions::split_non_actg), [`dna_packed`](ParserOptions::dna_packed), and [`dna_columnar`](ParserOptions::dna_columnar)).
//!
//! Sequence, header, and quality data are read directly from the parser methods ([`get_header`](HelicaseParser::get_header), [`get_dna_string`](HelicaseParser::get_dna_string), [`get_dna_packed`](HelicaseParser::get_dna_packed), [`get_quality`](HelicaseParser::get_quality), etc.).
//!
//! When only one type of event is active (such as the default configuration yielding `Event::Record`),
//! matching is unnecessary and errors can simply be propagated with `event?;`:
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default().config();
//!
//! # #[cfg(feature = "std")]
//! fn main() -> std::io::Result<()> {
//!     let mut parser = FastxParser::<CONFIG>::from_file("reads.fastq.gz")?;
//!
//!     while let Some(event) = parser.next() {
//!         event?;
//!         let header = parser.get_header();
//!         let seq = parser.get_dna_string();
//!     }
//!     Ok(())
//! }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//!
//! When both records and chunks are active (e.g. `dna_packed().return_record(true)`), match on the event using `?` to distinguish them:
//!
//! ```rust,no_run
//! use helicase::input::*;
//! use helicase::parser::Event;
//! use helicase::*;
//!
//! // dna_packed enables DnaChunk events; explicitly turn Record events back on.
//! const CONFIG: Config = ParserOptions::default().dna_packed().return_record(true).config();
//!
//! # #[cfg(feature = "std")]
//! fn main() -> std::io::Result<()> {
//!     let mut parser = FastxParser::<CONFIG>::from_file("reads.fastq.gz")?;
//!
//!     while let Some(event) = parser.next() {
//!         match event? {
//!             Event::DnaChunk(range) => {
//!                 // one contiguous ACTG run is ready
//!                 let seq = parser.get_dna_packed();
//!             }
//!             Event::Record(range) => {
//!                 // all chunks of this record have been processed
//!                 let header = parser.get_header();
//!             }
//!         }
//!     }
//!     Ok(())
//! }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//!
//! ## Iterating over chunks of packed DNA
//!
//! To parse sequences into contiguous 2-bit chunks:
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
//! # #[cfg(feature = "std")]
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!     let mut parser = FastxParser::<CONFIG>::from_file(&path)?;
//!
//!     // iterate over each chunk of ACTG bases
//!     while let Some(event) = parser.next() {
//!         event?;
//!         // headers are still accessible between chunks
//!         let header = parser.get_header();
//!
//!         // get a reference to the packed sequence
//!         let seq = parser.get_dna_packed();
//!
//!         // or directly get a PackedSeq (requires the packed-seq feature)
//!         // let packed_seq = parser.get_packed_seq();
//!     }
//!     Ok(())
//! }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//!
//! ## Parallel parsing
//!
//! For multi-threaded processing, [`ParallelFastxParser`] can parse files in parallel across worker threads.
//! Like [`FastxParser`], it can be created using [`from_file`](ParallelFastxParser::from_file) or [`from_file_mmap`](ParallelFastxParser::from_file_mmap) (see [Input](#input)).
//!
//! For compressed data, multi-block BGZF and multi-frame zstd files are automatically decompressed in parallel.
//!
//! Process records across threads with [`for_each`](ParallelFastxParser::for_each):
//!
//! ```rust,no_run
//! use helicase::config::*;
//! use helicase::*;
//!
//! const CONFIG: Config = ParserOptions::default().ignore_headers().config();
//!
//! # #[cfg(feature = "std")]
//! fn main() -> std::io::Result<()> {
//!     let path = "...";
//!
//!     ParallelFastxParser::<CONFIG>::from_file(&path)?
//!         .threads(8)
//!         .for_each(|mut parser| {
//!             while let Some(event) = parser.next() {
//!                 event?;
//!                 let header = parser.get_header();
//!                 let seq = parser.get_dna_string();
//!                 // ...
//!             }
//!             Ok(())
//!         })?;
//!
//!     Ok(())
//! }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//!
//! Or keep a state per thread with [`for_each_with`](ParallelFastxParser::for_each_with),
//! which returns the state of each thread:
//!
//! ```rust,no_run
//! # use helicase::config::*;
//! # use helicase::*;
//! # const CONFIG: Config = ParserOptions::default().ignore_headers().config();
//! # #[cfg(feature = "std")]
//! # fn main() -> std::io::Result<()> {
//! # let path = "...";
//! let total_len: usize = ParallelFastxParser::<CONFIG>::from_file(&path)?
//!     .threads(8)
//!     .for_each_with(
//!         || 0,
//!         |len, mut parser| {
//!             while let Some(event) = parser.next() {
//!                 event?;
//!                 *len += parser.get_dna_string().len();
//!             }
//!             Ok(())
//!         },
//!     )?
//!     .into_iter()
//!     .sum();
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "std"))]
//! # fn main() {}
//! ```
//!
//! # Crate features
//!
//! ## Standard library (`std`)
//!
//! Enabled by default and required for now. With `default-features = false`, enable it explicitly, e.g. `features = ["std", "zstd"]`.
//!
//! ## Packed-seq (`packed-seq`)
//!
//! Enables conversion from [`PackedDNA`](dna_format::PackedDNA) to `PackedSeq` from the [packed-seq](https://github.com/rust-seq/packed-seq) crate (disabled by default).
//! This is useful for [hashing *k*-mers](https://github.com/rust-seq/seq-hash) or [computing minimizers and syncmers](https://github.com/rust-seq/simd-minimizers).
//!
//! ## Decompression
//!
//! Helicase supports transparent decompression via [`from_file`](input::FromFile::from_file), automatically detecting the following formats:
//! - `gz` (default) for gzip and BGZF using [sirocco](https://github.com/srijs/sirocco)
//! - `zstd` (default) for zstd, including multi-frame files
//! - `bz2` for bzip2
//! - `xz` for xz
//!
//! ## No PDEP (`no-pdep`)
//!
//! By default, this library uses [PDEP](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) to compute the [PackedDNA format](#bitpacked-dna-formats).
//! However, this instruction can be slow on some older CPUs (especially AMD CPUs prior to 2020).
//! If you want an efficient implementation for these CPUs, we recommend using the `no-pdep` feature.
//!
//! ## Compatibility (`needletail`, `paraseq`)
//!
//! Helicase provides optional features to ease migration from `needletail` and `paraseq`.
//!
//! For `needletail`, it provides constructors like `parse_fastx_file`, whose `next()` yields a record with the usual methods (`id`, `seq`, `raw_seq`, `qual`, `num_bases`).
//!
//! For `paraseq`, it implements [`paraseq::Record`](https://docs.rs/paraseq/latest/paraseq/trait.Record.html) and provides [`process_parallel`](ParallelFastxParser::process_parallel) to run existing `ParallelProcessor` directly.
//!
//! See the [migration guide](https://github.com/imartayan/helicase/blob/main/MIGRATION.md) for details.

#![cfg_attr(docsrs, feature(doc_cfg))]
#![cfg_attr(not(any(feature = "std", test)), no_std)]
#![warn(missing_docs)]

extern crate alloc;

pub mod config;
pub mod dna_format;
pub mod input;
pub(crate) mod lexer;
#[cfg(feature = "std")]
pub mod parallel;
pub mod parser;

pub use config::{Config, ParserOptions};
#[cfg(feature = "std")]
pub use parallel::ParallelFastxParser;
pub use parser::{FastxParser, HelicaseParser};

#[cfg(feature = "needletail")]
pub use parser::needletail_helpers::{parse_fastx_file, parse_fastx_reader, parse_fastx_slice};

// Traits only implemented in this crate, so that they can gain methods in minor versions.
mod sealed {
    pub trait Sealed {}
}

pub(crate) mod simd;
