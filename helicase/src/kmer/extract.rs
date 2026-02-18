use crate::FastxParser;
use crate::config::*;
use crate::input::*;
use crate::kmer::*;
use crate::parser::*;

use std::io;

pub fn kmer_from_fastx_slice<const K: usize, B: BitStorage>(
    data: &[u8],
) -> io::Result<MerChunk<K, B>> {
    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;
    let mut chunks = MerChunk::<_, B>::new();
    while let Some(_) = parser.next() {
        // unsafe unwrap:: this unwrap actually panic if the DNA contains non-nuc symbols. This
        // can't happen per the semantic of the parser.
        let mut nchunks: MerChunk<_, B> = parser.get_dna_string().try_into().unwrap();
        chunks.append(&mut nchunks);
    }
    Ok(chunks)
}
