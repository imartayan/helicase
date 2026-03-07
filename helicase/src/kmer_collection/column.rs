use crate::config::advanced::*;
use crate::dna_format::ColumnarDNA;
use crate::*;
use std::cmp::min;
use std::io;
use std::thread;
use tracing::info;

#[derive(Clone, Debug)]
pub struct MerSlice<'a, const K: usize, B: BitStorage>(
    pub &'a [BitString<B, K>],
    pub &'a [BitString<B, K>],
);

impl<'a, const K: usize, B: BitStorage> MerSlice<'a, K, B> {
    pub fn get(&self, i: usize) -> Mer<K, B> {
        Mer::<K, B>(self.0[i], self.1[i])
    }

    /// Expose underlying storage slices
    pub fn as_slices(&self) -> (&'a [B], &'a [B]) {
        // UNSAFE: this is safe because a BitString is simply a wrapper around B.
        unsafe {
            let left: &[B] = std::slice::from_raw_parts(self.0.as_ptr() as *const B, self.0.len());
            let right: &[B] = std::slice::from_raw_parts(self.1.as_ptr() as *const B, self.1.len());
            (left, right)
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = Mer<K, B>> + '_ {
        self.0
            .iter()
            .zip(self.1.iter())
            .map(|(u, v)| Mer::<K, B>(*u, *v))
    }

    pub fn chunks(&self, size: usize) -> impl Iterator<Item = MerSlice<'_, K, B>> {
        let high = self.0.chunks(size);
        let low = self.1.chunks(size);
        high.zip(low).map(|(h, l)| MerSlice::<K, B>(h, l))
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

    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
        self.1.reserve(additional);
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

    #[inline(always)]
    pub fn append(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
        self.1.append(&mut other.1);
        debug_assert_eq!(self.0.len(), self.1.len());
    }

    pub fn clear(&mut self) {
        self.0.clear();
        self.1.clear();
    }

    pub fn truncate(&mut self, size: usize) {
        self.0.truncate(size);
        self.1.truncate(size);
    }

    pub fn push(&mut self, el: Mer<K, T>) {
        self.0.push(el.0);
        self.1.push(el.1);
    }

    // this could probably be improved !
    pub fn compute_hashes(&self) -> impl Iterator<Item = u64> {
        self.iter().map(|k| k.hash())
    }

    #[inline(always)]
    pub fn append_from_ascii(&mut self, seq: &[u8]) -> Result<(), String> {
        if seq.len() < K {
            return Ok(());
        }
        // check bound are already done, so the following is safe
        let first_chunk: &[u8; K] = &seq[0..K].try_into().unwrap();
        self.reserve(seq.len() - K + 1);
        let mut mer: Mer<K, T> = first_chunk.try_into()?;
        self.push(mer);
        for ch in &seq[K..] {
            mer = mer.append_right_ascii(*ch)?;
            self.push(mer)
        }
        Ok(())
    }

    #[inline(always)]
    pub fn append_from_columnar(&mut self, seq: &ColumnarDNA) -> Result<(), String> {
        if seq.len() < K {
            return Ok(());
        }
        // THIS IS FOR LOUP
        todo!();
    }

    pub fn par_sort(self, significant_bits: u8, dedup: bool, mut thread_count: usize) -> Self
    where
        T: Send + 'static,
    {
        let bucket_nb = 1 << significant_bits;
        thread_count = min(thread_count, bucket_nb);
        let shift = 64 - significant_bits;
        let mut buckets = vec![MerChunk::with_capacity(self.len() / bucket_nb); bucket_nb];
        for kmer in self.iter() {
            let buck_offset: usize = (kmer.hash() as usize) >> shift;
            buckets[buck_offset].push(kmer);
        }
        thread::scope(|s| {
            for chunk in buckets.chunks_mut(thread_count) {
                s.spawn(move || {
                    for merchunk in chunk {
                        merchunk.sort(dedup);
                    }
                });
            }
        });
        let mut chunks = MerChunk::new();
        for mut buc in buckets {
            chunks.append(&mut buc);
        }
        chunks
    }

    pub fn sort(&mut self, dedup: bool) {
        if self.is_empty() {
            return;
        }
        use super::permutation::*;
        info!("building the hash");
        let mut hash: Vec<(u64, usize)> =
            (0..self.len()).map(|a| (self.get(a).hash(), a)).collect();
        info!("Starting the sort");
        hash.sort_unstable();
        info!("Build the permutation");
        let mut onto_map: Vec<usize>; // This map the ith hash values to its offset 
        // in the original column.
        // So this is morally a representation of a permutation.
        //
        // If kmer are duplicated and dedup set to true, it becomes an onto-map.
        //

        let out_size;
        if dedup {
            info!("Dedup algorithm");
            onto_map = vec![0; self.len()];
            let mut curr_offset;
            let mut representatives: Vec<(_, usize)> = Vec::new();

            let first_idx = hash[0].1;
            let first_val = self.get(first_idx);
            onto_map[first_idx] = 0;
            representatives.push((first_val, 0));
            curr_offset = 0;

            for i in 1..hash.len() {
                let curr_idx = hash[i].1;
                let val = self.get(curr_idx);

                // If hash changed, clear representatives
                if hash[i - 1].0 != hash[i].0 {
                    representatives.clear();
                }

                // Check for duplicate in current hash run
                let mut found = None;
                for (rep_val, rep_off) in representatives.iter() {
                    if *rep_val == val {
                        found = Some(*rep_off);
                        break;
                    }
                }

                if let Some(rep_off) = found {
                    // duplicate → map to existing group
                    onto_map[curr_idx] = rep_off;
                    continue;
                }

                // new group
                curr_offset += 1;
                onto_map[curr_idx] = curr_offset;
                representatives.push((val, curr_offset));
            }

            out_size = curr_offset + 1;
        } else {
            onto_map = hash.into_iter().map(|(_, i)| i).collect();
            out_size = self.len();
        };
        info!("Applying mapping");
        self.0 = onto_map_copy(&mut self.0, &onto_map, out_size);
        self.1 = onto_map_copy(&mut self.1, &onto_map, out_size);
    }
}

impl<const K: usize, B: BitStorage> Extend<Mer<K, B>> for MerChunk<K, B> {
    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = Mer<K, B>>,
    {
        let iterator = iter.into_iter();
        let (lower, upper) = iterator.size_hint();
        if let Some(u) = upper {
            self.reserve(u);
        } else {
            self.reserve(lower);
        }
        for kmer in iterator {
            self.push(kmer);
        }
    }
}

impl<const K: usize, B: BitStorage> TryFrom<&[u8]> for MerChunk<K, B> {
    type Error = String;
    #[inline(always)]
    fn try_from(seq: &[u8]) -> Result<Self, Self::Error> {
        let mut mers = MerChunk::new();
        mers.append_from_ascii(seq)?;
        Ok(mers)
    }
}

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
    fn test_merchunk_append_utf8() {
        let seq1: &[u8] = b"ACGTA";
        let mut chunk = MerChunk::<4, u8>::new();
        chunk.append_from_ascii(seq1).unwrap();
        assert_eq!(chunk.len(), 2);
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
