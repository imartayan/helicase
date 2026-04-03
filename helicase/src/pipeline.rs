use crate::Config;
use crate::FastxParser;
use crate::HelicaseParser;
use crate::ParserOptions;
use crate::config::advanced::*;
use crate::input::FromSlice;
use crate::kmer::BitStorage;
use crate::kmer_collection::MerChunk;
use std::io;
use tracing::info;

pub fn kmer_from_fastx_slice<const K: usize, B: BitStorage>(
    data: &[u8],
) -> io::Result<MerChunk<K, B>> {
    const DNA_LEN: Config = COMPUTE_DNA_LEN | SPLIT_NON_ACTG | RETURN_DNA_CHUNK;
    let mut parser = FastxParser::<DNA_LEN>::from_slice(data).unwrap();
    let mut kmer_nb = 0;
    while let Some(_) = parser.next() {
        kmer_nb += parser.get_dna_len() - K + 1;
    }
    info!("Computed number of kmer before parsing: {}", kmer_nb);

    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .return_dna_chunk(true)
        .return_record(false)
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;
    let mut chunks = MerChunk::<_, B>::with_capacity(kmer_nb);
    while let Some(_) = parser.next() {
        // unsafe unwrap:: this unwrap actually panic if the DNA contains non-nuc symbols. This
        // can't happen per the semantic of the parser.
        let dna_str = parser.get_dna_string();
        chunks.append_from_ascii(dna_str).unwrap()
    }
    Ok(chunks)
}

pub fn chunk_process_from_fastx_slice<const K: usize, B: BitStorage>(
    data: &[u8],
    mut closure: impl FnMut(&mut MerChunk<K, B>) -> io::Result<()>,
) -> io::Result<()> {
    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .return_dna_chunk(true)
        .return_record(false)
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;
    let mut chunks = MerChunk::<_, B>::new();
    while let Some(_) = parser.next() {
        // unsafe unwrap:: this unwrap actually panic if the DNA contains non-nuc symbols. This
        // can't happen per the semantic of the parser.
        let dna_str = parser.get_dna_string();
        chunks
            .append_from_ascii(dna_str)
            .map_err(std::io::Error::other)?;
        closure(&mut chunks)?;
        chunks.clear();
    }
    Ok(())
}
