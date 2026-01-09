# Helicase

Helicase is a carefully optimized FASTA/FASTQ parser that extensively uses vectorized instructions.

It is designed for three main goals: being highly configurable, handling non-ACTG bases and computing bitpacked representations of DNA.

[Documentation](https://imartayan.github.io/helicase/)

## Requirements

This library requires AVX2 or NEON instruction sets, make sure to enable `target-cpu=native` when using it:
``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

## Usage

### Minimal example

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

The parser supports options that can be adjusted in the `ParserOptions`.
For instance, if you don't need to look at the headers and you want to skip non-ACTG bases, you can change to configuration to:
```rust
const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .skip_non_actg()
    .config();
```

### Bitpacked DNA formats

The parser can output a bitpacked representation of the sequence in two different formats:
- `PackedDNA` which maps each base to two bits and packs them (compatible with [packed-seq](https://github.com/rust-seq/packed-seq) using the corresponding [feature](#crate-features)).
- `ColumnarDNA` which separates the high bit and the low bit of each base, and store them in two bitmasks.

Since each base is encoded using two bits, we have to handle non-ACTG bases differently.
Three options are available for that:
- `split_non_actg` splits the sequence into contiguous chunks of ACTG bases, stopping the iterator at each chunk.
- `skip_non_actg` skips the non-ACTG bases and merge the remaining chunks together, stopping once at the end of the record.
- `keep_non_actg` keeps the non-ACTG bases and encodes them with a lossy representation.

### Iterating over chunks of packed DNA

```rust
use helicase::input::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default()
    .dna_packed()
    // don't stop the iterator at the end of a record
    .return_record(false)
    .config();

fn main() {
    let path = "...";

    let mut parser = FastxParser::<CONFIG>::from_file(&path).expect("Cannot open file");

    // iterate over each chunk of ACTG
    while let Some(_event) = parser.next() {
        // we still have access to the header
        let header = parser.get_header();

        // get a reference to the packed sequence
        let seq = parser.get_dna_packed();

        // or directly get a PackedSeq (requires the packed-seq feature)
        let packed_seq = parser.get_packed_seq();

        // ...
    }
}
```

## Crate features

### Decompression

This library supports transparent file decompression using [deko](https://github.com/igankevich/deko), you can choose the supported formats using the following features:
- `bz2` for bzip2 (disabled by default)
- `gz` for gzip (enabled by default)
- `xz` for xz (disabled by default)
- `zstd` for zstd (enabled by default)

### Packed-seq

The [PackedDNA format](#bitpacked-dna-formats) is compatible with [packed-seq](https://github.com/rust-seq/packed-seq) and can be converted when the `packed-seq` feature is enabled (disabled by default).

This can be useful for [hashing *k*-mers](https://github.com/rust-seq/seq-hash) or [computing minimizers & syncmers](https://github.com/rust-seq/simd-minimizers).

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

## Acknowledgements

This project was initially started by [Loup Lobet](https://lplt.net/) during his internship with [Charles Paperman](https://paperman.name/).
