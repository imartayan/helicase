# Changelog

<!-- next-header -->
## [Unreleased]

Helicase is now usable in WebAssembly and without `std`.

Breaking changes:
- internal input methods are now private

### What's New

- Add support for `no_std`
- Add support for simd128 (WASM)
- Show feature requirements on gated items in the documentation

## 1.0.0

On file inputs, the single-threaded parser is now 1.1-4× faster than needletail, while the new multi-threaded parser is 3-6× faster than paraseq using 4 threads.

Breaking change: parsers now yield `Result<Event, parser::Error>` (instead of `Event`), enabling standard `?` error propagation.

### What's New

- Built-in multi-threaded parser (`ParallelFastxParser`) with multiple readers.
- Streaming multi-reader BGZF (`.gz`) and multi-frame Zstandard (`.zst`) parallel decompression with constant memory (~2 MB).
- Accelerated 2-bit bitpacking (`PackedDNA` and `ColumnarDNA`) for contiguous reads and multiline sequences.
- AVX-512 vector kernel for 512-bit scanning on supported x86 CPUs.
- Zero-overhead CRLF (`\r\n`) line ending support for Windows.
- `reject_non_actg` option, returning an `ErrorKind::InvalidBase` error at the first non-ACTG base.
- `swap_dna_packed`, `swap_dna_columnar`, `swap_mask_non_actg`, `swap_mask_n` and `swap_packed_seq_vec` to get the current sequence without allocating and reuse the given buffer.
- Calling a `FastxParser` getter disabled by the config (e.g. `get_header` with `ignore_headers`) is now a compile-time error.
- `FastxParser` is `Send`.
- Needletail and paraseq migration helpers.

## 0.2.0

- Bump `packed-seq` to 5.0, associated types are incompatible with the previous version.
- Breaking change: `ParserOptions::dna_packed`/`dna_columnar`/`and_dna_packed`/`and_dna_columnar`/`split_non_actg` now disable `Record` events when splitting non-ACTG bases to avoid processing each sequence twice. This can be turned back on using `return_record`.

## 0.1.1

- Improve documentation

## 0.1.0

- Initial release
