# Benchmarks against other parsers

```
RUSTFLAGS="-C target-cpu=native" cargo r -r --bin bench -- [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file (FASTA/FASTQ)

Options:
  -r, --repeat <REPEAT>  Number of repetitions [default: 10]
  -c, --csv              Output in csv format
  -H, --no-csv-header    Do not include csv header
  -P, --no-perf          Disable perf metrics (Linux only)
  -B, --no-baseline      Disable baseline (needletail & paraseq)
  -S, --no-slice         Disable slice bench
  -f, --file             Enable (compressed) file bench
  -m, --mmap             Enable mmap bench
  -v, --show-val         Show result values (length, #records...)
  -h, --help             Print help
  -V, --version          Print version
```

## Perf counters permissions (Linux only)

On Linux, perf counters are used by default to collect more metrics.
If you don't have the permission to use them, you can either disable them (`-P`) or grant the corresponding permission (requires root):
```sh
sysctl -w kernel.perf_event_paranoid=1
```

## Output validation

If you'd like to verify that `helicase` produces the same output as [needletail](https://github.com/onecodex/needletail), a validation script is provided:
```sh
RUSTFLAGS="-C target-cpu=native" cargo r -r --example validate -- <INPUT>
```

## Minimizers benchmark

You can measure the performance gain of using `helicase` in [simd-minimizers](https://github.com/rust-seq/simd-minimizers) with:
```sh
RUSTFLAGS="-C target-cpu=native" cargo r -r -F simd-minimizers --example simd_mini -- <INPUT>
```
