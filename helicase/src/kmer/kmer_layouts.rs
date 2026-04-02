use crate::kmer::bitstring::*;
use std::cmp::min;
use std::fmt::Debug;

pub trait Kmer: Sized + Ord + Copy {
    fn append_left_ascii(&self, x: u8) -> Result<Self, String>;
    fn append_right_ascii(&self, x: u8) -> Result<Self, String>;
    fn complement(&self) -> Self;
    /// Hash implementation that will be encoded on the 2*K first bits
    /// where K is the kmer size.
    fn hash(&self) -> u64;
    fn normalize(&self) -> Self {
        min(*self, self.rc())
    }
    fn reverse(&self) -> Self;
    fn rc(&self) -> Self {
        self.reverse().complement()
    }
    fn to_string(self) -> String;

    /// Return a kmer whatever is the content. Usefull for starting iterating
    fn some_kmer() -> Self;
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

impl<const K: usize, B: BitStorage> Kmer for Mer<K, B> {
    #[inline(always)]
    fn append_left_ascii(&self, x: u8) -> Result<Self, String> {
        let nuc = Mer::<1, B>::from_ascii(x)?;
        Ok(self.push_left::<_>(&nuc))
    }

    #[inline(always)]
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
        let high = self.0.value.fold_to_u64();
        let low = self.1.value.fold_to_u64();
        let mut to_hash: u64;
        let mask: u64;
        if K < 32 {
            to_hash = (high << K) ^ low;
            mask = (!0) >> (64 - 2 * K);
        } else if K < 64 {
            to_hash = (high << (64 - K)) ^ low;
            mask = !0;
        } else {
            to_hash = high ^ (low.rotate_left(32));
            mask = !0;
        }

        to_hash ^= to_hash >> 30;
        to_hash *= 0xbf58476d1ce4e5b9 & mask;
        to_hash ^= to_hash >> 27;
        to_hash *= 0x94d049bb133111eb & mask;
        to_hash ^= to_hash >> 31;
        to_hash & mask
    }

    #[inline(always)]
    fn reverse(&self) -> Self {
        Self(self.0.reverse(), self.1.reverse())
    }

    #[inline(always)]
    fn to_string(self) -> String {
        // unsafe unchecked utf8 is safe here because we produce only ascii symbols
        unsafe {
            String::from_utf8_unchecked(self.nucleotids().map(|nuc| nuc.to_ascii()).collect())
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
    fn test_mer_nucleotids() {
        let mer: Mer<4, u8> = Mer::<4, u8>::some_kmer();
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
