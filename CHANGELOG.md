# Changelog

<!-- next-header -->

## git

## 0.2.0

- Bump `packed-seq` to 5.0, associated types are incompatible with the previous version.
- Breaking change: `ParserOptions::dna_packed`/`dna_columnar`/`and_dna_packed`/`and_dna_columnar`/`split_non_actg` now disable `Record` events when splitting non-ACTG bases to avoid processing each sequence twice. This can be turned back on using `return_record`.

## 0.1.1

Improve documentation

## 0.1.0

Initial release
