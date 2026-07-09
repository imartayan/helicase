# Helicase

[![crates.io](https://img.shields.io/crates/v/helicase)](https://crates.io/crates/helicase)
[![docs](https://img.shields.io/docsrs/helicase)](https://docs.rs/helicase)

Helicase is a carefully optimized FASTA/FASTQ parser that extensively uses vectorized instructions.

It is designed for three main goals: being highly configurable, handling non-ACTG bases and computing bitpacked representations of DNA.

The underlying algorithm is described in [Helicase: Vectorized parsing and bitpacking of genomic sequences](https://doi.org/10.64898/2026.03.19.712912), please [cite it](#citation) if you use this library.

## Requirements

This library requires AVX2, SSE3 or NEON instruction sets, make sure to enable `target-cpu=native` when using it:
``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

Note: if your CPU has a bad support for the [PDEP instruction](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) (e.g. AMD CPUs prior to 2020), it is recommended to use the `no-pdep` [feature](#crate-features):
``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release -F no-pdep
```

## Usage

### Minimal example

The main entry point is to define a configuration via `ParserOptions` and build a `FastxParser` with this configuration.

```rust
use helicase::input::*;
use helicase::*;

// set the options of the parser (at compile-time)
const CONFIG: Config = ParserOptions::default().config();

fn main() {
    let path = "...";

    // create a parser with the desired options
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    // iterate over records
    while let Some(_event) = parser.next() {
        // get a reference to the header
        let header = parser.get_header();

        // get a reference to the sequence (without newlines)
        let seq = parser.get_dna_string();

        // ...
    }
}
```

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

The parser can output a bitpacked representation of the sequence in two different formats:
- `PackedDNA` which maps each base to two bits and packs them (compatible with [packed-seq](https://github.com/rust-seq/packed-seq) using the corresponding [feature](#crate-features)).
- `ColumnarDNA` which separates the high bit and the low bit of each base, and store them in two bitmasks.

Since each base is encoded using two bits, we have to handle non-ACTG bases differently.
Three options are available via `ParserOptions`:
- `split_non_actg` splits the sequence at non-ACTG bases, yielding one `DnaChunk` event per contiguous ACTG run (default for bitpacked formats), and disables `Record` events by default to avoid processing each sequence twice.
- `skip_non_actg` skips non-ACTG bases and merges the remaining chunks, yielding one `Record` event per record.
- `keep_non_actg` keeps the non-ACTG bases and encodes them lossily, yielding one [`Record`](parser::Event::Record) event per record (default for string format).

### Events

The parser is an iterator that yields `Event` values.
An event signals a record boundary or a contiguous DNA chunk,
but the data is always read from the parser itself via `get_header`, `get_dna_string`, etc.

There are two kinds of event:
- `Event::Record`: emitted once per record, after all of its DNA chunks. Enabled by `return_record` (on by default, off by default whenever non-ACTG bases are split, to avoid processing each sequence twice).
- `Event::DnaChunk`: emitted for each contiguous ACTG run. Enabled by `return_dna_chunk` (on by default with `split_non_actg`).

When both are active you need to match on the event to distinguish them:
```rust
use helicase::input::*;
use helicase::parser::Event;
use helicase::*;

// dna_packed enables DnaChunk events; explicitly turn Record events back on.
const CONFIG: Config = ParserOptions::default().dna_packed().return_record(true).config();

fn main() {
    let path = "...";
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    while let Some(event) = parser.next() {
        match event {
            Event::DnaChunk(_) => {
                // one contiguous ACTG run is ready
                let seq = parser.get_dna_packed();
            }
            Event::Record(_) => {
                // all chunks of this record have been processed
            }
        }
    }
}
```

When only one type of event is active, the event value can be safely ignored:

```rust
use helicase::input::*;
use helicase::*;

// Default config: only Record events, one per record.
const CONFIG: Config = ParserOptions::default().config();

fn main() {
    let path = "...";
    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    while let Some(_event) = parser.next() {
        let header = parser.get_header();
        let seq = parser.get_dna_string();
    }
}
```

It is even possible to disable all events to process the entire file in one go, for instance if you simply want to count bases.

### Iterating over chunks of packed DNA

```rust
use helicase::input::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default()
    // by default, dna_packed splits non-ACTG bases, stops after each chunk, and doesn't stop at the end of a record
    .dna_packed()
    .config();

fn main() {
    let path = "...";

    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    // iterate over each chunk of ACTG bases
    while let Some(_event) = parser.next() {
        // headers are still accessible between chunks
        let header = parser.get_header();

        // get a reference to the packed sequence
        let seq = parser.get_dna_packed();

        // or directly get a PackedSeq (requires the packed-seq feature)
        // let packed_seq = parser.get_packed_seq();
    }
}
```

## Crate features

### Packed-seq

The [PackedDNA format](#bitpacked-dna-formats) is compatible with [packed-seq](https://github.com/rust-seq/packed-seq) and can be converted when the `packed-seq` feature is enabled (disabled by default).

This can be useful for [hashing *k*-mers](https://github.com/rust-seq/seq-hash) or [computing minimizers & syncmers](https://github.com/rust-seq/simd-minimizers).

### No PDEP

By default, this library uses [PDEP](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Parallel_bit_deposit_and_extract) to compute the [PackedDNA format](#bitpacked-dna-formats).
However, this instruction can be very slow on some CPUs (especially AMD CPUs prior to 2020).
If you want an efficient implementation for these CPUs, we recommend using the `no-pdep` feature.

### Decompression

This library supports transparent file decompression using [deko](https://github.com/igankevich/deko), you can choose the supported formats using the following features:
- `bz2` for bzip2
- `gz` for gzip (default)
- `xz` for xz
- `zstd` for zstd (default)

## Benchmarks

Benchmarks against [needletail](https://github.com/onecodex/needletail) and [paraseq](https://github.com/noamteyssier/paraseq) are available in the `bench` directory.
You can run them on any (possibly compressed) FASTA/FASTQ file using:
```sh
RUSTFLAGS="-C target-cpu=native" cargo r -r --bin bench -- <file>
```

For instance, you can run it on [this human genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz), [these short reads](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/NIST_Illumina_2x250bps/D1_S1_L001_R2_007.fastq.gz) or [these long reads](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/PacBio_HiFi/15kb/m54328_180928_230446.Q20.fastq).

Note that the FASTQ files can easily be converted to FASTA using:
```sh
RUSTFLAGS="-C target-cpu=native" cargo r -r --example fq_to_fa -- <file.fastq>
```

More information in the [bench README](bench/README.md).

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
