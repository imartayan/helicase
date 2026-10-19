//! Deterministic FASTA/FASTQ generator for testing and benchmarking.
//!
//! Generates synthetic genomic sequences with a seeded PRNG,
//! supporting standard formats and intentionally corrupted edge cases.

use rand::rngs::Xoshiro256PlusPlus;
use rand::{RngExt, SeedableRng};
use std::io::{self, Write};

/// Sequencing file format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

/// Character alphabet for sequence generation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Alphabet {
    /// Standard uppercase DNA bases: A, C, G, T
    Dna,
    /// Mixed-case DNA bases: A, C, G, T, a, c, g, t
    DnaMixed,
    /// Uppercase DNA bases including unknown base: A, C, G, T, N
    DnaWithN,
    /// Full IUPAC degenerate nucleotide codes: A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N
    DnaIupac,
    /// Standard 20 canonical amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    Protein,
    /// Custom user-defined alphabet
    Custom(Vec<u8>),
}

impl Alphabet {
    pub const DNA: &'static [u8] = b"ACGT";
    pub const DNA_WITH_N: &'static [u8] = b"ACGTN";
    pub const IUPAC: &'static [u8] = b"ACGTRYSWKMBDHVN";
    pub const DNA_MIXED: &'static [u8] = b"ACGTacgt";
    pub const PROTEIN: &'static [u8] = b"ACDEFGHIKLMNPQRSTVWY";

    pub const fn as_bytes(&self) -> &[u8] {
        match self {
            Self::Dna => Self::DNA,
            Self::DnaMixed => Self::DNA_MIXED,
            Self::DnaWithN => Self::DNA_WITH_N,
            Self::DnaIupac => Self::IUPAC,
            Self::Protein => Self::PROTEIN,
            Self::Custom(bytes) => bytes.as_slice(),
        }
    }

    /// Non-ACTG degenerate or unknown characters in this alphabet (`alphabet[4..]`).
    /// Returns empty slice for `Alphabet::Dna`, `Alphabet::Protein`, and `Alphabet::DnaMixed`.
    pub fn non_actg_bytes(&self) -> &[u8] {
        match self {
            Self::DnaWithN | Self::DnaIupac => &self.as_bytes()[4..],
            _ => b"",
        }
    }
}

/// Types of intentional corruption for error-handling tests, carrying occurrence probability `p` in `[0.0, 1.0]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Corrupt {
    /// First character of the record is invalid (e.g. `!` instead of `>` or `@`)
    InvalidHeaderPrefix(f64),
    /// FASTQ quality string length differs from sequence length
    MismatchedQualityLength(f64),
    /// FASTQ missing the `+` line
    MissingPlusLine(f64),
    /// Truncate the file abruptly inside this record
    TruncateMidRecord(f64),
}

impl Corrupt {
    #[inline(always)]
    pub const fn prob(&self) -> f64 {
        match *self {
            Self::InvalidHeaderPrefix(p)
            | Self::MismatchedQualityLength(p)
            | Self::MissingPlusLine(p)
            | Self::TruncateMidRecord(p) => p,
        }
    }
}

/// Unified generator and configuration for synthetic FASTA/FASTQ data.
#[derive(Clone, Debug)]
pub struct FastxGen {
    pub seed: u64,
    pub format: Format,
    pub num_records: usize,
    pub mean_len: usize,
    /// Relative standard deviation as a fraction of mean length (e.g. 0.1 for 10% variance, 0.0 for fixed).
    pub std_dev: f64,
    /// Emit Windows carriage returns (`\r\n`) instead of Unix (`\n`).
    pub crlf: bool,
    /// Wrap sequences across multiple lines (e.g. 60 or 80 cols). `None` for single-line.
    pub wrap_width: Option<usize>,
    pub alphabet: Alphabet,
    /// Proportion of non-ACTG bases (`alphabet[4..]`) injected into generated sequences (`[0.0, 1.0]`).
    /// Ignored for `Alphabet::Dna`.
    pub non_actg_prob: f64,
    /// Minimum Phred ASCII quality character (default `!` = 33).
    pub qual_min: u8,
    /// Maximum Phred ASCII quality character (default `~` = 126).
    pub qual_max: u8,
    pub header_prefix: String,
    /// Optional corruption mode and occurrence probability
    pub corrupt: Option<Corrupt>,
}

impl Default for FastxGen {
    fn default() -> Self {
        Self {
            seed: 42,
            format: Format::Fastq,
            num_records: 100,
            mean_len: 150,
            std_dev: 0.0,
            crlf: false,
            wrap_width: None,
            alphabet: Alphabet::Dna,
            non_actg_prob: 0.0,
            qual_min: b'!',
            qual_max: b'~',
            header_prefix: "record_".to_string(),
            corrupt: None,
        }
    }
}

impl FastxGen {
    /// Preset: Unwrapped FASTA sequences with fixed length.
    pub fn fasta(num_records: usize, len: usize) -> Self {
        Self {
            format: Format::Fasta,
            num_records,
            mean_len: len,
            ..Default::default()
        }
    }

    /// Preset: FASTQ records with fixed length.
    pub fn fastq(num_records: usize, len: usize) -> Self {
        Self {
            format: Format::Fastq,
            num_records,
            mean_len: len,
            ..Default::default()
        }
    }

    /// Preset: Illumina short reads (150bp fixed, FASTQ, Dna).
    pub fn illumina(num_records: usize) -> Self {
        Self::fastq(num_records, 150)
    }

    /// Preset: PacBio HiFi long reads (mean ~15kb, 10% variance, FASTQ, Dna).
    pub fn pacbio(num_records: usize) -> Self {
        Self::fastq(num_records, 15_000).with_std_dev(0.1)
    }

    /// Preset: Multi-line wrapped FASTA (e.g. chromosome/contig sequences wrapped at `wrap_width`).
    pub fn wrapped_fasta(num_records: usize, len: usize, wrap_width: usize) -> Self {
        Self::fasta(num_records, len).with_wrap(wrap_width)
    }

    /// Preset: Protein FASTA sequences.
    pub fn protein(num_records: usize, len: usize) -> Self {
        Self::fasta(num_records, len).with_alphabet(Alphabet::Protein)
    }

    /// Fluent builder method to override the random seed.
    pub const fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Fluent builder method to set relative length standard deviation (e.g. 0.1 for 10%).
    pub const fn with_std_dev(mut self, std_dev: f64) -> Self {
        self.std_dev = std_dev;
        self
    }

    /// Fluent builder method to set the mean sequence length.
    pub const fn with_mean_len(mut self, mean_len: usize) -> Self {
        self.mean_len = mean_len;
        self
    }

    /// Fluent builder method to set the sequence alphabet.
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        if self.non_actg_prob == 0.0 {
            match alphabet {
                Alphabet::DnaWithN => self.non_actg_prob = 0.05,
                Alphabet::DnaIupac => self.non_actg_prob = 0.05,
                _ => {}
            }
        }
        self.alphabet = alphabet;
        self
    }

    /// Fluent builder method to inject non-ACTG bases (`alphabet[4..]`) with a given probability.
    /// Ignored for `Alphabet::Dna`.
    pub const fn with_non_actg_prob(mut self, prob: f64) -> Self {
        self.non_actg_prob = prob;
        self
    }

    /// Fluent builder method to inject corrupted records for negative/resilience testing.
    pub const fn with_corrupt(mut self, corrupt: Corrupt) -> Self {
        self.corrupt = Some(corrupt);
        self
    }

    /// Fluent builder method to set multi-line sequence wrap width.
    pub const fn with_wrap(mut self, wrap_width: usize) -> Self {
        self.wrap_width = Some(wrap_width);
        self
    }

    /// Fluent builder method to configure the quality score range.
    pub const fn with_qual_range(mut self, min: u8, max: u8) -> Self {
        self.qual_min = min;
        self.qual_max = max;
        self
    }

    /// Fluent builder method to enable or disable Windows CRLF line endings.
    pub const fn with_crlf(mut self, crlf: bool) -> Self {
        self.crlf = crlf;
        self
    }

    /// Sample length for a single record using mean_len and relative std_dev.
    pub fn sample_len(&self, rng: &mut Xoshiro256PlusPlus) -> usize {
        if self.std_dev > 0.0 {
            let sigma = self.mean_len as f64 * self.std_dev;
            let n = sample_normal(rng, self.mean_len as f64, sigma);
            if n <= 1.0 { 1 } else { n.round() as usize }
        } else {
            self.mean_len.max(1)
        }
    }

    /// Generate a single ground-truth record.
    pub fn generate_record(&self, index: usize, rng: &mut Xoshiro256PlusPlus) -> GeneratedRecord {
        let id = format!("{}{}", self.header_prefix, index);
        let len = self.sample_len(rng);
        let mut seq = vec![0u8; len];
        fill_sequence(rng, &self.alphabet, &mut seq);

        let non_actg = self.alphabet.non_actg_bytes();
        if !non_actg.is_empty() && self.non_actg_prob > 0.0 {
            for byte in seq.iter_mut() {
                if rng.random_bool(self.non_actg_prob) {
                    *byte = non_actg[rng.random_range(0..non_actg.len())];
                }
            }
        }

        let qual = if self.format == Format::Fastq {
            let mut q = vec![0u8; len];
            fill_qualities(rng, self.qual_min, self.qual_max, &mut q);
            Some(q)
        } else {
            None
        };

        GeneratedRecord { id, seq, qual }
    }

    /// Iterator yielding ground-truth [`GeneratedRecord`] instances without file formatting.
    pub fn records(&self) -> RecordIterator {
        RecordIterator {
            generator: self.clone(),
            rng: Xoshiro256PlusPlus::seed_from_u64(self.seed),
            index: 0,
        }
    }

    /// Write formatted FASTA or FASTQ output to any writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(self.seed);
        let le: &[u8] = if self.crlf { b"\r\n" } else { b"\n" };
        let wrap = self.wrap_width.filter(|&w| w > 0).unwrap_or(usize::MAX);

        for i in 0..self.num_records {
            let GeneratedRecord { id, seq, qual } = self.generate_record(i, &mut rng);

            let active_corrupt = match self.corrupt {
                Some(c) if rng.random_bool(c.prob()) => Some(c),
                _ => None,
            };

            // Header
            match active_corrupt {
                Some(Corrupt::InvalidHeaderPrefix(_)) => writer.write_all(b"!")?,
                _ => match self.format {
                    Format::Fasta => writer.write_all(b">")?,
                    Format::Fastq => writer.write_all(b"@")?,
                },
            }
            writer.write_all(id.as_bytes())?;
            writer.write_all(le)?;

            // Sequence (wrapped or unwrapped)
            for chunk in seq.chunks(wrap) {
                writer.write_all(chunk)?;
                writer.write_all(le)?;
            }

            // FASTQ '+' and quality
            if self.format == Format::Fastq {
                if !matches!(active_corrupt, Some(Corrupt::MissingPlusLine(_))) {
                    writer.write_all(b"+")?;
                    writer.write_all(le)?;
                }

                if let Some(mut q) = qual {
                    if matches!(active_corrupt, Some(Corrupt::MismatchedQualityLength(_))) {
                        if q.len() > 1 {
                            q.pop();
                        } else {
                            q.push(b'I');
                        }
                    }
                    for chunk in q.chunks(wrap) {
                        writer.write_all(chunk)?;
                        writer.write_all(le)?;
                    }
                }
            }

            if matches!(active_corrupt, Some(Corrupt::TruncateMidRecord(_))) {
                return Ok(());
            }
        }

        Ok(())
    }

    /// Generate complete dataset directly into a `Vec<u8>`.
    pub fn generate_bytes(&self) -> Vec<u8> {
        let estimated_size = self.num_records * (self.mean_len * 2 + 64);
        let mut buf = Vec::with_capacity(estimated_size);
        self.write_to(&mut buf).expect("In-memory write failed");
        buf
    }

    /// Generate complete dataset as a `String`.
    pub fn generate_string(&self) -> String {
        String::from_utf8(self.generate_bytes()).expect("Generated non-UTF8 text")
    }
}

/// A single ground-truth record generated by [`FastxGen`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneratedRecord {
    /// Record header identifier (excluding leading `>` or `@`).
    pub id: String,
    /// Raw un-wrapped sequence bytes.
    pub seq: Vec<u8>,
    /// Raw un-wrapped quality bytes (if FASTQ).
    pub qual: Option<Vec<u8>>,
}

/// Iterator over ground-truth records.
pub struct RecordIterator {
    generator: FastxGen,
    rng: Xoshiro256PlusPlus,
    index: usize,
}

impl Iterator for RecordIterator {
    type Item = GeneratedRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.generator.num_records {
            return None;
        }
        let rec = self.generator.generate_record(self.index, &mut self.rng);
        self.index += 1;
        Some(rec)
    }
}

/// Fast batch table unpacker extracting `BITS` bits per letter from 64-bit random words.
#[inline(always)]
fn unpack_lut<const BITS: usize, const N: usize>(
    rng: &mut Xoshiro256PlusPlus,
    lut: &[u8],
    out: &mut [u8],
) {
    let mask = (1u64 << BITS) - 1;
    for chunk in out.chunks_mut(N) {
        let mut bits = rng.random::<u64>();
        for byte in chunk {
            *byte = lut[(bits & mask) as usize % lut.len()];
            bits >>= BITS;
        }
    }
}

/// High-performance sequence buffer filler.
#[inline]
pub fn fill_sequence(rng: &mut Xoshiro256PlusPlus, alphabet: &Alphabet, out: &mut [u8]) {
    match alphabet {
        Alphabet::Dna | Alphabet::DnaWithN | Alphabet::DnaIupac => {
            unpack_lut::<2, 32>(rng, Alphabet::DNA, out)
        }
        Alphabet::DnaMixed => unpack_lut::<3, 21>(rng, Alphabet::DNA_MIXED, out),
        Alphabet::Protein => unpack_lut::<5, 12>(rng, Alphabet::PROTEIN, out),
        Alphabet::Custom(bytes) => unpack_lut::<5, 12>(rng, bytes, out),
    }
}

/// High-performance FASTQ quality buffer filler unpacking 8 quality scores per 64-bit word.
#[inline]
pub fn fill_qualities(rng: &mut Xoshiro256PlusPlus, qual_min: u8, qual_max: u8, out: &mut [u8]) {
    let range = (qual_max - qual_min + 1) as u16;
    for chunk in out.chunks_mut(8) {
        let bytes = rng.random::<u64>().to_le_bytes();
        for (byte, b) in chunk.iter_mut().zip(bytes) {
            *byte = qual_min + (((b as u16 * range) >> 8) as u8);
        }
    }
}

/// Sample from a normal distribution using the Box-Muller transform.
fn sample_normal(rng: &mut Xoshiro256PlusPlus, mean: f64, std_dev: f64) -> f64 {
    if std_dev <= 0.0 {
        return mean;
    }
    let u1: f64 = rng.random_range(1e-7..1.0);
    let u2: f64 = rng.random_range(1e-7..1.0);
    let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * core::f64::consts::PI * u2).cos();
    z0.mul_add(std_dev, mean)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_deterministic() {
        let g1 = FastxGen::fasta(5, 20).with_seed(12345).generate_bytes();
        let g2 = FastxGen::fasta(5, 20).with_seed(12345).generate_bytes();
        assert_eq!(g1, g2);
        assert!(g1.starts_with(b">record_0\n"));
    }

    #[test]
    fn test_fastq_deterministic() {
        let g1 = FastxGen::illumina(3).with_seed(999).generate_bytes();
        let g2 = FastxGen::illumina(3).with_seed(999).generate_bytes();
        assert_eq!(g1, g2);
        assert!(g1.starts_with(b"@record_0\n"));
    }

    #[test]
    fn test_record_iterator_matches_content() {
        let generator = FastxGen::illumina(10).with_seed(42);
        let records: Vec<_> = generator.records().collect();
        assert_eq!(records.len(), 10);
        assert_eq!(records[0].id, "record_0");
        assert_eq!(
            records[0].seq.len(),
            records[0].qual.as_ref().unwrap().len()
        );
    }

    #[test]
    fn test_fasta_wrapped() {
        let generator = FastxGen::wrapped_fasta(1, 100, 40).with_seed(777);
        let bytes = generator.generate_bytes();
        let lines: Vec<_> = bytes.split(|&b| b == b'\n').collect();
        assert_eq!(lines[0], b">record_0");
        assert_eq!(lines[1].len(), 40);
        assert_eq!(lines[2].len(), 40);
        assert_eq!(lines[3].len(), 20);
    }

    #[test]
    fn test_crlf_line_endings() {
        let mut generator = FastxGen::fasta(2, 50).with_seed(101);
        generator.crlf = true;
        let bytes = generator.generate_bytes();
        assert!(bytes.contains(&b'\r'));
        assert!(bytes.windows(2).any(|w| w == b"\r\n"));
    }

    #[test]
    fn test_protein_alphabet() {
        let generator = FastxGen::protein(5, 200).with_seed(123);
        let valid_aa = b"ACDEFGHIKLMNPQRSTVWY";
        for rec in generator.records() {
            assert_eq!(rec.seq.len(), 200);
            assert!(rec.seq.iter().all(|b| valid_aa.contains(b)));
        }
    }

    #[test]
    fn test_dna_alphabet_pure() {
        let generator = FastxGen::fasta(5, 500).with_seed(456);
        let valid_dna = b"ACGT";
        for rec in generator.records() {
            assert_eq!(rec.seq.len(), 500);
            assert!(rec.seq.iter().all(|b| valid_dna.contains(b)));
        }
    }

    #[test]
    fn test_dna_alphabet_ignores_non_actg_prob() {
        // Even when non_actg_prob is high, Alphabet::Dna remains 100% pure ACGT
        let generator = FastxGen::fasta(10, 1000)
            .with_non_actg_prob(0.5)
            .with_seed(123);
        for rec in generator.records() {
            assert!(rec.seq.iter().all(|b| b"ACGT".contains(b)));
        }
    }

    #[test]
    fn test_non_actg_prob_dna_with_n() {
        let generator = FastxGen::fasta(10, 1000)
            .with_alphabet(Alphabet::DnaWithN)
            .with_non_actg_prob(0.1) // 10% Ns
            .with_seed(333);
        let records: Vec<_> = generator.records().collect();
        let total_len: usize = records.iter().map(|r| r.seq.len()).sum();
        let n_count: usize = records
            .iter()
            .flat_map(|r| &r.seq)
            .filter(|&&b| b == b'N')
            .count();
        let n_ratio = n_count as f64 / total_len as f64;
        assert!(
            n_ratio > 0.08 && n_ratio < 0.12,
            "Expected ~10% Ns, got {:.2}%",
            n_ratio * 100.0
        );
        for rec in records {
            assert!(rec.seq.iter().all(|b| b"ACGTN".contains(b)));
        }
    }

    #[test]
    fn test_non_actg_prob_iupac_draws_from_alphabet_slice() {
        let generator = FastxGen::fasta(20, 1000)
            .with_alphabet(Alphabet::DnaIupac)
            .with_non_actg_prob(0.2) // 20% degenerate bases
            .with_seed(777);
        let records: Vec<_> = generator.records().collect();
        let non_actg_letters: std::collections::HashSet<u8> = records
            .iter()
            .flat_map(|r| &r.seq)
            .copied()
            .filter(|b| !b"ACGT".contains(b))
            .collect();

        // Must draw across all degenerate letters in alphabet[4..], not just 'N'
        let deg_bases = Alphabet::DnaIupac.non_actg_bytes();
        for &deg in deg_bases {
            assert!(
                non_actg_letters.contains(&deg),
                "Missing expected degenerate base '{}' from IUPAC drawing",
                deg as char
            );
        }
    }

    #[test]
    fn test_sample_len() {
        let fixed = FastxGen::fasta(10, 200);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        assert_eq!(fixed.sample_len(&mut rng), 200);

        let var = FastxGen::fasta(100, 1000).with_std_dev(0.1);
        let lens: Vec<usize> = (0..100).map(|_| var.sample_len(&mut rng)).collect();
        let avg = lens.iter().sum::<usize>() as f64 / lens.len() as f64;
        assert!(avg > 950.0 && avg < 1050.0);
        assert!(lens.iter().any(|&l| l != 1000));
    }
}
