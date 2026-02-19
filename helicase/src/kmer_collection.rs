use crate::*;
use std::io;

#[derive(Clone, Debug)]
pub struct MerSlice<'a, const K: usize, B: BitStorage>(
    pub &'a [BitString<B, K>],
    pub &'a [BitString<B, K>],
);

impl<'a, const K: usize, B: BitStorage> MerSlice<'a, K, B> {
    pub fn get(&self, i: usize) -> Mer<K, B> {
        Mer::<K, B>(self.0[i], self.1[i])
    }

    pub fn iter(&self) -> impl Iterator<Item = Mer<K, B>> + '_ {
        self.0
            .iter()
            .zip(self.1.iter())
            .map(|(u, v)| Mer::<K, B>(*u, *v))
    }
}

#[derive(Clone, Debug)]
pub struct MerChunk<const K: usize, B: BitStorage>(
    pub Vec<BitString<B, K>>,
    pub Vec<BitString<B, K>>,
);

impl<const K: usize, B: BitStorage> Default for MerChunk<K, B> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const K: usize, T: BitStorage> MerChunk<K, T> {
    pub fn new() -> Self {
        Self(vec![], vec![])
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self(
            Vec::<_>::with_capacity(capacity),
            Vec::<_>::with_capacity(capacity),
        )
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn as_slice(&self) -> MerSlice<'_, K, T> {
        MerSlice::<K, T>(&self.0, &self.1)
    }

    pub fn iter(&self) -> impl Iterator<Item = Mer<K, T>> + '_ {
        self.0
            .iter()
            .zip(self.1.iter())
            .map(|(u, v)| Mer::<K, T>(*u, *v))
    }

    pub fn get(&self, i: usize) -> Mer<K, T> {
        self.as_slice().get(i)
    }

    pub fn append(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
        self.1.append(&mut other.1);
        debug_assert_eq!(self.0.len(), self.1.len());
    }

    pub fn push(&mut self, el: Mer<K, T>) {
        self.0.push(el.0);
        self.1.push(el.1);
    }
}

impl<const K: usize, B: BitStorage> TryFrom<&[u8]> for MerChunk<K, B> {
    type Error = String;
    #[inline(always)]
    fn try_from(seq: &[u8]) -> Result<Self, Self::Error> {
        if seq.len() < K {
            return Ok(MerChunk::new());
        }
        let mut mers = MerChunk::<K, B>::with_capacity(seq.len() - K + 1);
        // check bound are already done, so the following is safe
        let first_chunk: &[u8; K] = &seq[0..K].try_into().unwrap();
        let mut mer: Mer<K, B> = first_chunk.try_into()?;
        mers.push(mer);
        for ch in &seq[K..] {
            mer = mer.append_right_ascii(*ch)?;
            mers.push(mer)
        }
        Ok(mers)
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_merchunk_from_slice() {
        let seq: &[u8] = b"ACGTAC";
        let chunk: MerChunk<4, u8> = seq.try_into().unwrap();
        assert_eq!(chunk.len(), 3); // "ACGT", "CGTA", "GTAC"
        let strings: Vec<String> = chunk.iter().map(|m| m.to_string()).collect();
        assert_eq!(strings.len(), 3);
    }

    #[test]
    fn test_merchunk_append_push() {
        let seq1: &[u8] = b"ACGT";
        let mut chunk1: MerChunk<4, u8> = seq1.try_into().unwrap();
        let len1 = chunk1.len();
        assert_eq!(len1, 1);
        let seq2: &[u8] = b"TGCA";
        let mut chunk2: MerChunk<4, u8> = seq2.try_into().unwrap();
        assert_eq!(chunk2.len(), 1);

        chunk1.append(&mut chunk2);
        assert_eq!(chunk1.len(), 2);

        let mer = Mer::<4, u8>::some_kmer();
        chunk1.push(mer);
        assert_eq!(chunk1.len(), 3);
    }
}
