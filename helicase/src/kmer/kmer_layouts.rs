use crate::dna_format::ColumnarDNA;
use crate::kmer::bitstring::*;
use std::cmp::min;
use std::fmt::Debug;
use std::marker::PhantomData;

/// Owned container for producing a kmer iterator from owned DNA data
pub struct DNAOwnedIter<I> {
    dna: ColumnarDNA,
    _phantom: PhantomData<I>, // tie the iterator type to the struct
}

impl<I> DNAOwnedIter<I> {
    /// Create a new owned wrapper
    pub fn new(dna: ColumnarDNA) -> Self {
        Self {
            dna,
            _phantom: PhantomData,
        }
    }

    /// Borrow the underlying DNA for zero-copy iteration
    pub fn as_dna(&self) -> &ColumnarDNA {
        &self.dna
    }
}

// Implement IntoIterator for borrowed DNAOwnedIter to produce zero-copy iterator
impl<'a, K> IntoIterator for &'a DNAOwnedIter<K>
where
    K: Kmer,
{
    type Item = K;
    type IntoIter = K::Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        K::from_dna_columns(&self.dna)
    }
}

pub trait Kmer: Sized + Ord + Eq + Copy {
    type Iter<'a>: Iterator<Item = Self>;
    /// The hash type used by this kmer
    type Hash: Copy + Ord + Eq;

    fn complement(&self) -> Self;
    /// Hash implementation that will be encoded on the 2*K first bits
    /// where K is the kmer size.
    fn hash(&self) -> Self::Hash;
    fn normalize(&self) -> Self {
        min(*self, self.rc())
    }
    fn reverse(&self) -> Self;
    fn rc(&self) -> Self {
        self.reverse().complement()
    }
    fn to_string(&self) -> String;

    /// Produce an iterator by converting the input into columnar format.
    /// The iterator own its data, hence the 'static lifetime.
    fn from_str_iter<S: AsRef<str>>(input: S) -> DNAOwnedIter<Self::Iter<'static>> {
        let input_ref = input.as_ref();
        let mut dna = ColumnarDNA::with_capacity(input_ref.len());
        dna.push_str(input_ref);
        DNAOwnedIter::new(dna)
    }

    fn from_dna_columns<'a>(input: &'a ColumnarDNA) -> Self::Iter<'a>;

    /// Return a kmer whatever is the content. Usefull for starting iterating
    fn any() -> Self;
}

fn reversible_hash<const K: usize, I: Into<u64>>(high: I, low: I) -> u64 {
    let mut to_hash: u64;

    let high_64: u64 = high.into();
    let low_64: u64 = low.into();
    to_hash = (high_64 << K) ^ low_64;
    let mask: u64 = (!0) >> (64 - 2 * K);
    to_hash ^= to_hash >> 30;
    to_hash *= 0xbf58476d1ce4e5b9 & mask;
    to_hash ^= to_hash >> 27;
    to_hash *= 0x94d049bb133111eb & mask;
    to_hash ^= to_hash >> 31;
    to_hash & mask
}

#[derive(Eq, PartialEq, PartialOrd, Ord, Copy, Clone, Debug, Hash)]
pub struct ReversibleHashMer<const K: usize>(u64);

impl<const K: usize> From<Mer<K, u32>> for ReversibleHashMer<K> {
    fn from(mer: Mer<K, u32>) -> Self {
        Self(mer.hash())
    }
}

pub struct DNAIterReversibleHash<'a, const K: usize> {
    dna: &'a ColumnarDNA,
    offset: usize,
}

impl<'a, const K: usize> DNAIterReversibleHash<'a, K> {
    fn new(dna: &'a ColumnarDNA) -> Self {
        Self { dna, offset: 0 }
    }
}

impl<'a, const K: usize> Iterator for DNAIterReversibleHash<'a, K> {
    type Item = ReversibleHashMer<K>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.offset + K <= self.dna.len() {
            let (high, low) = self.dna.extract_u32(self.offset, K);
            Some(ReversibleHashMer(reversible_hash::<K, u32>(high, low)))
        } else {
            None
        }
    }
}

impl<const K: usize> From<ReversibleHashMer<K>> for Mer<K, u32> {
    fn from(_item: ReversibleHashMer<K>) -> Self {
        // YOAAAAAAAAAANNNNNNNNNNNNN
        todo!()
    }
}

impl<const K: usize> Kmer for ReversibleHashMer<K>
where
    Self: Into<Mer<K, u32>>,
{
    type Iter<'a> = DNAIterReversibleHash<'a, K>;
    type Hash = u64;
    fn hash(&self) -> Self::Hash {
        self.0
    }

    fn complement(&self) -> Self {
        let mer: Mer<K, u32> = (*self).into();
        mer.complement().into()
    }

    fn reverse(&self) -> Self {
        let mer: Mer<K, u32> = (*self).into();
        mer.reverse().into()
    }

    fn to_string(&self) -> String {
        let mer: Mer<K, u32> = (*self).into();
        mer.to_string()
    }

    fn from_dna_columns<'a>(input: &'a ColumnarDNA) -> Self::Iter<'a> {
        Self::Iter::new(input)
    }

    fn any() -> Self {
        Self(0)
    }
}

#[derive(Eq, PartialEq, PartialOrd, Ord, Copy, Clone, Debug, Hash)]
pub struct Mer<const K: usize, B: BitStorage>(pub BitString<B, K>, pub BitString<B, K>);
impl<const K: usize, B: BitStorage> Mer<K, B> {
    const ZERO: BitString<B, K> = BitString::<B, K>::ZERO;
    const ONE: BitString<B, K> = BitString::<B, K>::ONE;
    #[inline(always)]
    fn push_left<const L: usize>(&self, other: &Mer<L, B>) -> Self {
        Self(
            self.0.push_left::<_>(&other.0),
            self.1.push_left::<_>(&other.1),
        )
    }

    #[inline(always)]
    fn push_right<const L: usize>(&self, other: &Mer<L, B>) -> Self {
        Self(
            self.0.push_right::<_>(&other.0),
            self.1.push_right::<_>(&other.1),
        )
    }

    fn nucleotids(&self) -> impl Iterator<Item = Mer<1, B>> {
        (0..K).map(|i| Mer::<1, B>(self.0.get_bitstring::<1>(i), self.1.get_bitstring::<1>(i)))
    }

    #[inline(always)]
    pub fn append_left_ascii(&self, x: u8) -> Result<Self, String> {
        let nuc = Mer::<1, B>::from_ascii(x)?;
        Ok(self.push_left::<_>(&nuc))
    }

    #[inline(always)]
    pub fn append_right_ascii(&self, x: u8) -> Result<Self, String> {
        let nuc = Mer::<1, B>::from_ascii(x)?;
        Ok(self.push_right::<_>(&nuc))
    }
}

impl<B: BitStorage> Mer<1, B> {
    const A: Self = Self(Self::ZERO, Self::ZERO);
    const C: Self = Self(Self::ZERO, Self::ONE);
    const G: Self = Self(Self::ONE, Self::ONE);
    const T: Self = Self(Self::ONE, Self::ZERO);

    fn to_ascii(self) -> u8 {
        match (self.0.get(0), self.1.get(0)) {
            (false, false) => b'A',
            (false, true) => b'C',
            (true, true) => b'G',
            (true, false) => b'T',
        }
    }

    #[inline(always)]
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

pub struct DNAIterMer<'a, const K: usize, B: BitStorage> {
    dna: &'a ColumnarDNA,
    offset: usize,
    _phantom: PhantomData<B>,
}

impl<'a, const K: usize, B: BitStorage> DNAIterMer<'a, K, B> {
    fn new(dna: &'a ColumnarDNA) -> Self {
        Self {
            dna,
            offset: 0,
            _phantom: PhantomData,
        }
    }
}

impl<'a, const K: usize, B: BitStorage> Iterator for DNAIterMer<'a, K, B> {
    type Item = Mer<K, B>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.offset + K <= self.dna.len() {
            let (high, low) = BitString::<B, K>::dna_extract(self.dna, self.offset);
            Some(Mer::<K, B>(high, low))
        } else {
            None
        }
    }
}

impl<const K: usize, B: BitStorage> Kmer for Mer<K, B> {
    type Hash = u64;
    type Iter<'a> = DNAIterMer<'a, K, B>;

    #[inline(always)]
    fn complement(&self) -> Self {
        Self(self.0.not(), self.1)
    }

    #[inline(always)]
    fn hash(&self) -> u64 {
        let high = self.0.value.fold_to_u64();
        let low = self.1.value.fold_to_u64();
        reversible_hash::<K, u64>(high, low)
    }

    #[inline(always)]
    fn reverse(&self) -> Self {
        Self(self.0.reverse(), self.1.reverse())
    }

    fn from_dna_columns<'a>(input: &'a ColumnarDNA) -> Self::Iter<'a> {
        Self::Iter::new(input)
    }

    #[inline(always)]
    fn to_string(&self) -> String {
        // unsafe unchecked utf8 is safe here because we produce only ascii symbols
        unsafe {
            String::from_utf8_unchecked(self.nucleotids().map(|nuc| nuc.to_ascii()).collect())
        }
    }

    #[inline(always)]
    fn any() -> Self {
        Self(Self::ZERO, Self::ZERO)
    }
}

impl<const K: usize, B: BitStorage> TryFrom<&[u8; K]> for Mer<K, B> {
    type Error = String;
    #[inline(always)]
    fn try_from(seq: &[u8; K]) -> Result<Self, Self::Error> {
        let mut mer = Mer::<K, B>::any();
        for c in seq {
            mer = mer.append_right_ascii(*c)?;
        }
        Ok(mer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let mer = Mer::<4, u8>::any();
        let mer = mer.append_right_ascii(b'A').unwrap();
        let mer = mer.append_left_ascii(b'T').unwrap();
        let s = mer.to_string();
        assert!(s.contains('A') || s.contains('T'));
    }

    #[test]
    fn test_mer_complement() {
        let mer = Mer::<4, u8>::any();
        let c = mer.complement();
        assert_ne!(mer.0, c.0);
    }

    #[test]
    fn test_mer_reverse() {
        let mer = Mer::<4, u8>::any();
        let rev = mer.reverse();
        assert_eq!(rev.reverse(), mer);
    }

    #[test]
    fn test_mer_rc_and_normalize() {
        let mer = Mer::<4, u8>::any();
        let rc = mer.rc();
        let norm = mer.normalize();
        assert_eq!(norm, std::cmp::min(mer, rc));
    }

    #[test]
    fn test_mer_nucleotids() {
        let mer: Mer<4, u8> = Mer::<4, u8>::any();
        let count = mer.nucleotids().count();
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
    fn test_hash_uses_only_2k_bits_small_k() {
        fn check<const K: usize, T: BitStorage>() {
            let seq = [b'A'; K];
            let mer: Mer<K, T> = (&seq).try_into().unwrap();
            let h = mer.hash();

            let mask = if 2 * K == 64 {
                u64::MAX
            } else {
                (1u64 << (2 * K)) - 1
            };

            assert_eq!(h & !mask, 0, "hash has bits outside 2K region for K={}", K);
        }

        check::<5, u8>();
        check::<5, u128>();
        check::<10, u16>();
        check::<16, u16>();
        check::<16, u32>();
        check::<31, u32>();
        check::<31, u64>();
        check::<31, u128>();
    }

    #[test]
    fn test_hash_uses_only_2k_bits_varied_sequences() {
        fn check<const K: usize>(seq: [u8; K]) {
            let mer: Mer<K, u8> = (&seq).try_into().unwrap();
            let h = mer.hash();
            let offset: usize = 64 - 2 * K;

            let mask = if 2 * K == 64 {
                u64::MAX
            } else {
                (1u64 << (2 * K)) - 1
            };

            assert_eq!(
                h & !mask,
                0,
                "hash:\n {h:064b}\n{:>offset$}^\n has bits outside 2K region for sequence {:?} on K={K}",
                offset,
                seq,
            );
        }

        check::<4>([b'A', b'C', b'G', b'T']);
        check::<4>([b'T', b'T', b'G', b'A']);
        check::<8>([b'A', b'A', b'C', b'G', b'T', b'C', b'G', b'T']);
    }
}
