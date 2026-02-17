//use crate::dna::ColumnarDNA;
use crate::kmer::bitstring::*;
use std::cmp::min;
use std::fmt::Debug;

pub trait Kmer: Sized + Ord + Copy {
    fn append_left_ascii(&self, x: u8) -> Result<Self, String>;
    fn append_right_ascii(&self, x: u8) -> Result<Self, String>;
    fn complement(&self) -> Self;
    fn hash(&self) -> u64;
    fn normalize(&self) -> Self {
        min(*self, self.rc())
    }
    fn reverse(&self) -> Self;
    fn rc(&self) -> Self {
        self.reverse().complement()
    }
    fn to_string(&self) -> String;

    /// Return a kmer whatever is the content. Usefull for starting iterating
    fn some_kmer() -> Self;
}

#[derive(Eq, PartialEq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct Mer<const K: usize, B: BitStorage>(BitString<B, K>, BitString<B, K>);
impl<const K: usize, B: BitStorage> Mer<K, B> {
    const ZERO: BitString<B, K> = BitString::<B, K>::ZERO;
    const ONE: BitString<B, K> = BitString::<B, K>::ONE;

    fn push_left<const L: usize>(&self, other: &Mer<L, B>) -> Self {
        Self(
            self.0.push_left::<_>(&other.0),
            self.1.push_left::<_>(&other.1),
        )
    }

    fn push_right<const L: usize>(&self, other: &Mer<L, B>) -> Self {
        Self(
            self.0.push_right::<_>(&other.0),
            self.1.push_right::<_>(&other.1),
        )
    }

    fn into_nucleotids(&self) -> impl Iterator<Item = Mer<1, B>> {
        (0..K).map(|i| Mer::<1, B>(self.0.get_bitstring::<1>(i), self.1.get_bitstring::<1>(i)))
    }
}

impl<B: BitStorage> Mer<1, B> {
    const A: Self = Self(Self::ZERO, Self::ZERO);
    const C: Self = Self(Self::ZERO, Self::ONE);
    const G: Self = Self(Self::ONE, Self::ONE);
    const T: Self = Self(Self::ONE, Self::ZERO);

    fn to_ascii(&self) -> u8 {
        match (self.0.get(0), self.1.get(0)) {
            (false, false) => b'A',
            (false, true) => b'C',
            (true, true) => b'G',
            (true, false) => b'T',
        }
    }

    fn from_ascii(ch: u8) -> Result<Self, String> {
        match ch {
            b'A' | b'a' => Ok(Self::A),
            b'C' | b'c' => Ok(Self::C),
            b'G' | b'g' => Ok(Self::G),
            b'T' | b't' => Ok(Self::T),
            _ => Err(format!("Invalid nucleotide: {}", ch as char)),
        }
    }
}

impl<const K: usize, B: BitStorage> Kmer for Mer<K, B> {
    fn append_left_ascii(&self, x: u8) -> Result<Self, String> {
        let nuc = Mer::<1, B>::from_ascii(x)?;
        Ok(self.push_left::<_>(&nuc))
    }

    fn append_right_ascii(&self, x: u8) -> Result<Self, String> {
        let nuc = Mer::<1, B>::from_ascii(x)?;
        Ok(self.push_right::<_>(&nuc))
    }

    #[inline(always)]
    fn complement(&self) -> Self {
        Self(self.0.not(), self.1)
    }

    #[inline(always)]
    fn hash(&self) -> u64 {
        self.0.hash() ^ (self.1.hash() << 3) // BOUUH
    }

    #[inline(always)]
    fn reverse(&self) -> Self {
        Self(self.0.reverse(), self.1.reverse())
    }

    #[inline(always)]
    fn to_string(&self) -> String {
        // unsafe unchecked utf8 is safe here because we produce only ascii symbols
        unsafe {
            String::from_utf8_unchecked(self.into_nucleotids().map(|nuc| nuc.to_ascii()).collect())
        }
    }

    #[inline(always)]
    fn some_kmer() -> Self {
        Self(Self::ZERO, Self::ZERO)
    }
}

impl<const K: usize, B: BitStorage> TryFrom<&[u8; K]> for Mer<K, B> {
    type Error = String;
    #[inline(always)]
    fn try_from(seq: &[u8; K]) -> Result<Self, Self::Error> {
        let mut mer = Mer::<K, B>::some_kmer();
        for c in seq {
            mer = mer.append_right_ascii(*c)?;
        }
        Ok(mer)
    }
}

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

//pub fn kmer_from_slice_fastx<const K: usize, B: BitStorage>(data: &[u8]) -> MerChunk<K, B> {
//    let parser = FastxParser::<MINIMAL>::from_slice(data).unwrap();
//    // HERE
//}

#[cfg(test)]
mod tests {
    use super::*;

    type B5 = Mer<5, u8>;
    type B8 = Mer<8, u8>;
    type B16 = Mer<16, u16>;

    #[test]
    fn test_mer_from_ascii() {
        let a = Mer::<1, u8>::from_ascii(b'A').unwrap();
        assert_eq!(a.to_ascii(), b'A');
        let g = Mer::<1, u8>::from_ascii(b'G').unwrap();
        assert_eq!(g.to_ascii(), b'G');
        assert!(Mer::<1, u8>::from_ascii(b'X').is_err());
    }

    #[test]
    fn test_mer_append() {
        let mer = Mer::<4, u8>::some_kmer();
        let mer = mer.append_right_ascii(b'A').unwrap();
        let mer = mer.append_left_ascii(b'T').unwrap();
        let s = mer.to_string();
        assert!(s.contains('A') || s.contains('T'));
    }

    #[test]
    fn test_mer_complement() {
        let mer = Mer::<4, u8>::some_kmer();
        let c = mer.complement();
        assert_ne!(mer.0, c.0);
    }

    #[test]
    fn test_mer_reverse() {
        let mer = Mer::<4, u8>::some_kmer();
        let rev = mer.reverse();
        assert_eq!(rev.reverse(), mer);
    }

    #[test]
    fn test_mer_rc_and_normalize() {
        let mer = Mer::<4, u8>::some_kmer();
        let rc = mer.rc();
        let norm = mer.normalize();
        assert_eq!(norm, std::cmp::min(mer, rc));
    }

    #[test]
    fn test_mer_into_nucleotids() {
        let mer: Mer<4, u8> = Mer::<4, u8>::some_kmer();
        let count = mer.into_nucleotids().count();
        assert_eq!(count, 4);
    }

    #[test]
    fn test_mer_try_from_array() {
        let seq = [b'A', b'C', b'G', b'T'];
        let mer: Mer<4, u8> = (&seq).try_into().unwrap();
        let s = mer.to_string();
        assert_eq!(s.len(), 4);
        assert!(s.contains('A') || s.contains('C') || s.contains('G') || s.contains('T'));
    }

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
