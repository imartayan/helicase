# Migrating to Helicase

This guide explains how to migrate code from `needletail` or `paraseq` to `helicase`.

Both libraries have an optional drop-in feature that lets you reuse existing code with minimal changes, as well as a native API.

## Needletail

### Drop-in (`needletail` feature)

Add `helicase` with the `needletail` feature:

```sh
cargo add helicase -F needletail
```

This provides `parse_fastx_file`, `parse_fastx_reader`, `parse_fastx_slice`, and aliases for common record methods (`id()`, `seq()`, `raw_seq()`, `qual()`, `num_bases()`):

```rust
use helicase::parse_fastx_file;
use helicase::*;

fn main() -> std::io::Result<()> {
    let mut reader = parse_fastx_file("reads.fastq.gz")?;
    let mut total_bases = 0;

    while let Some(_event) = reader.next() {
        let header = reader.id();
        let seq = reader.seq();
        total_bases += reader.num_bases();

        if let Some(qual) = reader.qual() {
            // ...
        }
    }

    println!("Total bases: {total_bases}");
    Ok(())
}
```

Unlike needletail, `next()` advances the parser in-place rather than allocating a separate record struct, so methods are called on the reader itself.

### Native API

With native helicase, options are configured at compile-time via `ParserOptions`:

```rust
use helicase::config::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default().config();

fn main() -> std::io::Result<()> {
    let mut parser = FastxParser::<CONFIG>::from_file("reads.fastq.gz")?;
    let mut total_bases = 0;

    while let Some(_event) = parser.next() {
        let header = parser.get_header();
        let seq = parser.get_dna_string();
        total_bases += parser.get_dna_len();

        if let Some(qual) = parser.get_quality() {
            // ...
        }
    }

    println!("Total bases: {total_bases}");
    Ok(())
}
```

Key points:
- Sequence and header methods return borrowed `&[u8]` slices directly into the parser buffer.
- Unused fields can be disabled in `CONFIG` (e.g. `.ignore_headers()`, `.ignore_quality()`) so the compiler skips parsing them entirely.
- 2-bit bitpacking is built-in via `.dna_packed()`:
  ```rust
  const CONFIG: Config = ParserOptions::default().dna_packed().config();

  while let Some(_event) = parser.next() {
      let packed = parser.get_dna_packed();
      let words = packed.as_slice();
  }
  ```

---

## Paraseq

### Drop-in (`paraseq` feature)

Add `helicase` with the `paraseq` feature:

```sh
cargo add helicase -F paraseq
```

Existing `ParallelProcessor` implementations work without changes via `process_parallel`:

```rust
use helicase::config::*;
use helicase::*;
use paraseq::parallel::ParallelProcessor;
use paraseq::Record;

#[derive(Clone, Default)]
struct MyProcessor { /* ... */ }

impl<R: Record> ParallelProcessor<R> for MyProcessor {
    fn process_record(&mut self, record: R) -> Result<(), paraseq::Error> {
        let seq = record.seq_raw();
        // ...
        Ok(())
    }
}

fn main() -> std::io::Result<()> {
    const CONFIG: Config = ParserOptions::default().config();
    let mut processor = MyProcessor::default();

    ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
        .threads(8)
        .process_parallel(&mut processor)?;

    Ok(())
}
```

Helicase implements `paraseq::Record` for `&FastxParser`, so `record.seq_raw()` returns borrowed zero-copy slices.

### Native API (`for_each` and `reduce`)

Helicase native parallel parsing uses closures instead of processor traits.

#### `reduce` for aggregations

Instead of atomic counters or shared mutexes, use `reduce` with a thread-local accumulator:

```rust
use helicase::config::*;
use helicase::*;

const CONFIG: Config = ParserOptions::default().ignore_headers().config();

fn main() -> std::io::Result<()> {
    let (total_records, total_len) = ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
        .threads(8)
        .reduce(
            |mut parser| {
                let mut records = 0;
                let mut len = 0;
                while let Some(_event) = parser.next() {
                    records += 1;
                    len += parser.get_dna_len();
                }
                (records, len)
            },
            |acc, val| (acc.0 + val.0, acc.1 + val.1),
        )?;

    println!("Records: {total_records}, Total len: {total_len}");
    Ok(())
}
```

#### `for_each` for batching

To collect results into a shared buffer across threads:

```rust
use helicase::config::*;
use helicase::*;
use std::sync::Mutex;

const CONFIG: Config = ParserOptions::default().ignore_headers().config();

fn main() -> std::io::Result<()> {
    let dest = Mutex::new(Vec::new());

    ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
        .threads(8)
        .for_each(|mut parser| {
            let mut local = Vec::with_capacity(1024);
            while let Some(_event) = parser.next() {
                local.push(parser.get_dna_len());
            }
            dest.lock().unwrap().extend_from_slice(&local);
        })?;

    let all_lengths = dest.into_inner().unwrap();
    Ok(())
}
```
