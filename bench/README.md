# Helicase Benchmark Suite

Benchmarks comparing Helicase against Needletail and Paraseq across in-memory slices, streaming file I/O, memory-mapped files, and multi-threaded scaling.

## Running benchmarks

Run without arguments to execute the benchmark suite on built-in synthetic datasets (Illumina reads, PacBio long reads, wrapped FASTA, and large contigs):

```sh
RUSTFLAGS="-C target-cpu=native" cargo bench
# or
RUSTFLAGS="-C target-cpu=native" cargo r -r --bin bench
```

Or benchmark on a custom FASTA or FASTQ file:

```sh
RUSTFLAGS="-C target-cpu=native" cargo bench -- <FILE>
```

### Common options

```
Arguments:
  [FILE]                  Input FASTA or FASTQ file (optional; runs synthetic suite if omitted)

Options:
  -r, --repeat <REPEAT>   Number of repetitions [default: 5]
  -w, --warmup <WARMUP>   Number of warmup iterations [default: 1]
  -p, --parallel [N]      Run multi-threaded benchmarks (e.g. -p 4 or -p 2,4,8) vs Paraseq
  -b, --bitpack           Detailed 2-bit bitpacking comparison vs Needletail + PackedSeqVec
  -F, --file              Add streaming FileInput benchmark (with 256 KiB buffer)
  -M, --mmap              Add memory-mapped MmapInput benchmark
  -B, --no-baseline       Disable baseline parsers (Needletail & Paraseq)
  -f, --format <FORMAT>   Output format: text, markdown, csv, json [default: text]
  -h, --help              Print help
```

### Examples

```sh
# Benchmark 2-bit bitpacking throughput
cargo bench -- -b

# Benchmark parallel scaling on 2, 4, and 8 threads
cargo bench -- -p 2,4,8

# Benchmark streaming file and mmap I/O on a custom file
cargo bench -- -F -M /path/to/reads.fastq.gz
```

## Minimizers benchmark

You can measure the performance gain of using `helicase` in [simd-minimizers](https://github.com/rust-seq/simd-minimizers) with:
```sh
RUSTFLAGS="-C target-cpu=native" cargo r -r -p bench --features simd-minimizers --example simd_mini -- <FILE>
```

## Legacy benchmark mode

The legacy benchmark runner evaluates individual files with optional CSV output and Linux hardware performance counters:

```sh
cargo bench -- legacy [OPTIONS] <FILE>
# or
RUSTFLAGS="-C target-cpu=native" cargo r -r --bin bench -- legacy [OPTIONS] <FILE>
```

### Options

```
Arguments:
  <FILE>                 Input FASTA or FASTQ file

Options:
  -r, --repeat <REPEAT>  Number of repetitions [default: 5]
  -w, --warmup <WARMUP>  Number of warmup iterations [default: 1]
  -c, --csv              Output in CSV format (for pipelines/plotting)
  -H, --no-csv-header    Suppress CSV header
  -B, --no-baseline      Disable baseline parsers (Needletail & Paraseq)
  -S, --no-slice         Disable in-memory slice parsing
  -F, --file             Enable streaming FileInput benchmark
  -M, --mmap             Enable memory-mapped MmapInput benchmark
  -P, --no-perf          Disable Linux hardware perf metrics (instructions, branches)
  -v, --show-val         Show result values (length, record count, verification)
  -h, --help             Print help
```

### Perf counters permissions (Linux only)

On Linux, hardware performance counters can be used in legacy mode to collect low-level CPU metrics (instructions, CPU cycles, branch mispredictions).
If you don't have permission to read them, either disable them with `-P` or grant permission (requires root):
```sh
sysctl -w kernel.perf_event_paranoid=1
```
