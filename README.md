# Helicase

[![crates.io](https://img.shields.io/crates/v/helicase)](https://crates.io/crates/helicase)
[![docs](https://img.shields.io/docsrs/helicase)](https://docs.rs/helicase)

Helicase is a carefully optimized, SIMD-accelerated FASTA/FASTQ parser for DNA sequences.

It is designed for speed and flexibility: zero-overhead compile-time configuration, robust handling of non-ACTG bases, and on-the-fly 2-bit bitpacking.

It is typically 1.2-2.4× faster than `needletail` on a single thread, and 2.5-6× faster than `paraseq` across 4 threads.

The underlying algorithm is described in [Helicase: Vectorized parsing and bitpacking of genomic sequences](https://doi.org/10.64898/2026.03.19.712912), please [cite it](#citation) if you use this library.

## Requirements

This library requires AVX-512, AVX2, SSE3 or NEON instruction sets, make sure to enable `target-cpu=native` when using it:
``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

Note: if your CPU has a bad support for the [PDEP instruction](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) (e.g. AMD CPUs prior to 2020), it is recommended to use the `no-pdep` [feature](#crate-features).

## Usage

If you are migrating from `needletail` or `paraseq`, see the [migration guide](MIGRATION.md).

### Minimal example

The main entry point is to define a configuration via `ParserOptions` and build a `FastxParser`/`ParallelFastxParser` with this configuration.

```rust
use helicase::input::*;
use helicase::*;

// set the options of the parser (at compile-time)
const CONFIG: Config = ParserOptions::default().config();

fn main() -> std::io::Result<()> {
    let path = "...";

    // create a parser with the desired options
    let mut parser = FastxParser::<CONFIG>::from_file(&path)?;

    // iterate over records
    while let Some(event) = parser.next() {
        event?;
        // get a reference to the header
        let header = parser.get_header();

        // get a reference to the sequence (without newlines)
        let seq = parser.get_dna_string();

        // ...
    }
    Ok(())
}
```

### Input

Files are opened with `from_file` (which supports [transparent decompression](#crate-features)) or `from_file_mmap` (which only supports uncompressed data but can be faster).
In-memory slices (`from_slice`) and generic streams (`from_reader`) are also supported.

### Adjusting the configuration

The parser is configured at compile-time via `ParserOptions`.
For example, to ignore headers and split non-ACTG bases:
```rust
const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .split_non_actg()
    .config();
```

### Bitpacked DNA formats

The parser can output a bitpacked representation of the sequence in two formats:
- `PackedDNA` maps each base to two bits and packs them (compatible with [packed-seq](https://github.com/rust-seq/packed-seq) via the `packed-seq` [feature](#crate-features)).
- `ColumnarDNA` separates the high bit and the low bit of each base into two bitmasks.

Since 2-bit encoding only represents four bases, non-ACTG characters are handled explicitly.
Three options are available via `ParserOptions`:
- `split_non_actg` splits the sequence at non-ACTG bases, yielding one `DnaChunk` event per contiguous ACTG run (default for bitpacked formats), and disables `Record` events by default to avoid processing each sequence twice.
- `skip_non_actg` skips non-ACTG bases and merges the remaining chunks, yielding one `Record` event per record.
- `keep_non_actg` keeps the non-ACTG bases and encodes them lossily, yielding one `Record` event per record (default for string format).

### Events

The parser is an iterator that yields `Result<Event, parser::Error>` values.
An event signals a record boundary or a contiguous DNA chunk:
- `Event::Record(range)`: emitted once per record boundary. The `range` gives the byte offsets `start..end` of the entire record. Enabled by `return_record` (on by default for string and lossy formats).
- `Event::DnaChunk(range)`: emitted for each contiguous run of ACTG bases when non-ACTG characters are split. The `range` gives the byte offsets `start..end` of the chunk. Enabled by `return_dna_chunk` (on by default with `split_non_actg`, `dna_packed`, and `dna_columnar`).

Sequence, header, and quality data are read directly from the parser methods (`get_header`, `get_dna_string`, `get_dna_packed`, `get_quality`, etc.).

When only one type of event is active (such as the default configuration yielding `Event::Record`), matching is unnecessary and errors can simply be propagated with `event?;`:

```rust
while let Some(event) = parser.next() {
    event?;
    let header = parser.get_header();
    let seq = parser.get_dna_string();
}
```

When both records and chunks are active (e.g. `dna_packed().return_record(true)`), match on the event using `?` to distinguish them:

```rust
use helicase::input::*;
use helicase::parser::Event;
use helicase::*;

// dna_packed enables DnaChunk events; explicitly turn Record events back on.
const CONFIG: Config = ParserOptions::default().dna_packed().return_record(true).config();

fn main() -> std::io::Result<()> {
    let mut parser = FastxParser::<CONFIG>::from_file("reads.fastq.gz")?;

    while let Some(event) = parser.next() {
        match event? {
            Event::DnaChunk(range) => {
                // one contiguous ACTG run is ready
                let seq = parser.get_dna_packed();
            }
            Event::Record(range) => {
                // all chunks of this record have been processed
                let header = parser.get_header();
            }
        }
    }
    Ok(())
}
```

### Iterating over chunks of packed DNA

To parse sequences into contiguous 2-bit chunks:

```rust
use helicase::input::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default()
    // by default, dna_packed splits non-ACTG bases, stops after each chunk, and doesn't stop at the end of a record
    .dna_packed()
    .config();

fn main() -> std::io::Result<()> {
    let path = "...";
    let mut parser = FastxParser::<CONFIG>::from_file(&path)?;

    // iterate over each chunk of ACTG bases
    while let Some(event) = parser.next() {
        event?;
        // headers are still accessible between chunks
        let header = parser.get_header();

        // get a reference to the packed sequence
        let seq = parser.get_dna_packed();

        // or directly get a PackedSeq (requires the packed-seq feature)
        // let packed_seq = parser.get_packed_seq();
    }
    Ok(())
}
```

### Parallel parsing

For multi-threaded processing, `ParallelFastxParser` can parse files in parallel across worker threads.
Like `FastxParser`, it can be created using `from_file` or `from_file_mmap` (see [Input](#input)).

For compressed data, multi-block BGZF and multi-frame zstd files are automatically decompressed in parallel.

Process records across threads with `for_each`:

```rust
use helicase::config::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default().ignore_headers().config();

fn main() -> std::io::Result<()> {
    let path = "...";

    ParallelFastxParser::<CONFIG>::from_file(&path)?
        .threads(8)
        .for_each(|mut parser| {
            while let Some(event) = parser.next() {
                event?;
                let header = parser.get_header();
                let seq = parser.get_dna_string();
                // ...
            }
            Ok(())
        })?;

    Ok(())
}
```

Or aggregate results across threads with `reduce`:

```rust
let total_len: usize = ParallelFastxParser::<CONFIG>::from_file(&path)?
    .threads(8)
    .reduce(
        |mut parser| {
            let mut local_len = 0;
            while let Some(event) = parser.next() {
                event?;
                local_len += parser.get_dna_string().len();
            }
            Ok(local_len)
        },
        |acc, val| acc + val,
    )?;
```

## Crate features

### Packed-seq (`packed-seq`)

Enables conversion from [`PackedDNA`](#bitpacked-dna-formats) to `PackedSeq` from the [packed-seq](https://github.com/rust-seq/packed-seq) crate (disabled by default).
This is useful for [hashing *k*-mers](https://github.com/rust-seq/seq-hash) or [computing minimizers and syncmers](https://github.com/rust-seq/simd-minimizers).

### Decompression

Helicase supports transparent decompression via `from_file`, automatically detecting the following formats:
- `gz` (default) for gzip and BGZF using [sirocco](https://github.com/srijs/sirocco)
- `zstd` (default) for zstd, including multi-frame files
- `bz2` for bzip2
- `xz` for xz

### No PDEP (`no-pdep`)

By default, this library uses [PDEP](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) to compute the [PackedDNA format](#bitpacked-dna-formats).
However, this instruction can be slow on some older CPUs (especially AMD CPUs prior to 2020).
If you want an efficient implementation for these CPUs, we recommend using the `no-pdep` feature.

### Compatibility (`needletail`, `paraseq`)

Helicase provides optional features to ease migration from `needletail` and `paraseq`.

For `needletail`, it provides drop-in method aliases (`id`, `seq`, `raw_seq`, `qual`) and constructors like `parse_fastx_file`.

For `paraseq`, it implements `paraseq::Record` and provides `process_parallel` to run existing `ParallelProcessor` directly.

See the [migration guide](MIGRATION.md) for details.

## Benchmarks

A benchmark suite comparing Helicase against [needletail](https://github.com/onecodex/needletail) and [paraseq](https://github.com/noamteyssier/paraseq) is included.

Run without arguments to execute benchmarks on built-in synthetic datasets:
```sh
RUSTFLAGS="-C target-cpu=native" cargo bench
```

Or benchmark on any custom (optionally compressed) FASTA/FASTQ file:
```sh
RUSTFLAGS="-C target-cpu=native" cargo bench -- <file>
```

You can test it with public datasets such as:
- the [T2T Human Genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz) (FASTA)
- some [Illumina short reads](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/NIST_Illumina_2x250bps/D1_S1_L001_R2_007.fastq.gz) (FASTQ)
- some [PacBio HiFi long reads](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/PacBio_HiFi/15kb/m54328_180928_230446.Q20.fastq) (FASTQ)

More options (parallel scaling, bitpacking comparisons, decompression, format) are documented in the [bench README](bench/README.md).

## Acknowledgements

This project was initially started by [Loup Lobet](https://lplt.net/) during his internship with [Charles Paperman](https://paperman.name/).

## Citation

> Helicase: Vectorized parsing and bitpacking of genomic sequences. Igor Martayan, Loup Lobet, Camille Marchet, and Charles Paperman. https://doi.org/10.64898/2026.03.19.712912

```bibtex
@article{helicase,
  title     = {Helicase: Vectorized parsing and bitpacking of genomic sequences},
  author    = {Martayan, Igor and Lobet, Loup and Marchet, Camille and Paperman, Charles},
  year      = {2026},
  month     = {03},
  publisher = {Cold Spring Harbor Laboratory},
  journal   = {bioRxiv},
  doi       = {10.64898/2026.03.19.712912}
}
```
