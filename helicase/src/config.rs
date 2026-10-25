//! Compile-time configuration of the parser.

/// Configuration for the parser, represented using bitflags.
pub type Config = u64;

pub mod advanced {
    //! Bitflags used for the configuration.

    use super::*;

    /// Returns `true` if any bit of `flag` is set in `config`.
    #[inline(always)]
    pub const fn flag_is_set(config: Config, flag: Config) -> bool {
        config & flag != 0
    }

    /// Returns `true` if no bit of `flag` is set in `config`.
    #[inline(always)]
    pub const fn flag_is_not_set(config: Config, flag: Config) -> bool {
        config & flag == 0
    }

    /// Configuration of [`ParserOptions::default`]: headers and DNA strings.
    pub const DEFAULT_CONFIG: Config = COMPUTE_HEADER | COMPUTE_DNA_STRING | RETURN_RECORD;

    /// See [`ParserOptions::compute_headers`].
    pub const COMPUTE_HEADER: Config = 1 << 0;
    /// See [`ParserOptions::and_dna_string`].
    pub const COMPUTE_DNA_STRING: Config = 1 << 1;
    /// See [`ParserOptions::and_dna_columnar`].
    pub const COMPUTE_DNA_COLUMNAR: Config = 1 << 2;
    /// See [`ParserOptions::and_dna_packed`].
    pub const COMPUTE_DNA_PACKED: Config = 1 << 3;
    /// See [`ParserOptions::compute_dna_len`].
    pub const COMPUTE_DNA_LEN: Config = 1 << 4;
    /// See [`ParserOptions::compute_quality`].
    pub const COMPUTE_QUALITY: Config = 1 << 5;
    /// See [`ParserOptions::split_non_actg`].
    pub const SPLIT_NON_ACTG: Config = 1 << 6;
    /// See [`ParserOptions::return_record`].
    pub const RETURN_RECORD: Config = 1 << 7;
    /// See [`ParserOptions::return_dna_chunk`].
    pub const RETURN_DNA_CHUNK: Config = 1 << 8;
    /// Append the chunks of a record instead of replacing them, see
    /// [`clear_chunk`](crate::HelicaseParser::clear_chunk).
    pub const MERGE_DNA_CHUNKS: Config = 1 << 9;
    /// See [`ParserOptions::merge_records`]. Implies [`MERGE_DNA_CHUNKS`].
    pub const MERGE_RECORDS: Config = 1 << 10;
    /// See [`ParserOptions::compute_mask_non_actg`].
    pub const COMPUTE_MASK_NON_ACTG: Config = 1 << 11;
    /// Compute a [`BitMask`](crate::dna_format::BitMask) of N bases, see
    /// [`get_mask_n`](crate::HelicaseParser::get_mask_n).
    pub const COMPUTE_MASK_N: Config = 1 << 12;
    /// See [`ParserOptions::reject_non_actg`].
    pub const REJECT_NON_ACTG: Config = 1 << 13;

    /// `true` if bitpacking uses PDEP (BMI2 target without the `no-pdep` feature).
    pub const PDEP_ENABLED: bool = cfg!(all(target_feature = "bmi2", not(feature = "no-pdep")));
}

use advanced::*;

/// Compile-time builder for the configuration of the parser.
///
/// Options are set using a builder pattern and must be finalized with
/// [`config`](ParserOptions::config) to produce a [`Config`] usable as a
/// const generic.
///
/// # Examples
///
/// Default configuration - compute headers and DNA as ASCII-encoded bytes, stopping after each record:
///
/// ```rust
/// use helicase::*;
///
/// const CONFIG: Config = ParserOptions::default().config();
/// ```
///
/// Ignore headers and compute [`PackedDNA`](crate::dna_format::PackedDNA), splitting non-ACTG characters and stopping after each chunk by default:
///
/// ```rust
/// use helicase::*;
///
/// const CONFIG: Config = ParserOptions::default()
///     .ignore_headers()
///     .dna_packed()
///     .config();
/// ```
///
/// Compute both a DNA string and [`PackedDNA`](crate::dna_format::PackedDNA), splitting non-ACTG characters and stopping after each chunk by default:
///
/// ```rust
/// use helicase::*;
///
/// const CONFIG: Config = ParserOptions::default()
///     .dna_string()
///     .and_dna_packed()
///     .config();
/// ```
///
/// Compute [`ColumnarDNA`](crate::dna_format::ColumnarDNA), lossily encode non-ACTG bases and produce a [`BitMask`](crate::dna_format::BitMask) marking their positions:
///
/// ```rust
/// use helicase::*;
///
/// const CONFIG: Config = ParserOptions::default()
///     .dna_columnar()
///     .keep_non_actg()
///     .compute_mask_non_actg()
///     .config();
/// ```
#[derive(Clone, Copy)]
pub struct ParserOptions(Config);

impl ParserOptions {
    /// Creates a default configuration, which computes headers and DNA as bytes.
    #[inline(always)]
    pub const fn default() -> Self {
        Self(DEFAULT_CONFIG)
    }

    /// Load an existing configuration.
    #[inline(always)]
    pub const fn from_config(config: Config) -> Self {
        Self(config)
    }

    /// Build the configuration of the parser.
    #[inline(always)]
    pub const fn config(self) -> Config {
        self.0
    }

    /// Enable a flag, typically one of [`advanced`] such as [`COMPUTE_MASK_N`].
    /// Call it after the other options, some of them reset flags.
    #[inline(always)]
    pub const fn with_flag(self, flag: Config) -> Self {
        Self(self.0 | flag)
    }

    /// Disable a flag, call it after the other options like [`with_flag`](Self::with_flag).
    #[inline(always)]
    pub const fn without_flag(self, flag: Config) -> Self {
        Self(self.0 & !flag)
    }

    /// Enable the computation of headers (default).
    #[inline(always)]
    pub const fn compute_headers(self) -> Self {
        self.with_flag(COMPUTE_HEADER)
    }

    /// Disable the computation of headers.
    #[inline(always)]
    pub const fn ignore_headers(self) -> Self {
        self.without_flag(COMPUTE_HEADER)
    }

    /// Enable the computation of quality.
    #[inline(always)]
    pub const fn compute_quality(self) -> Self {
        self.with_flag(COMPUTE_QUALITY)
    }

    /// Disable the computation of quality (default).
    #[inline(always)]
    pub const fn ignore_quality(self) -> Self {
        self.without_flag(COMPUTE_QUALITY)
    }

    /// Enable the computation DNA length.
    #[inline(always)]
    pub const fn compute_dna_len(self) -> Self {
        self.with_flag(COMPUTE_DNA_LEN)
    }

    /// Disable the computation DNA length (default).
    #[inline(always)]
    pub const fn ignore_dna_len(self) -> Self {
        self.without_flag(COMPUTE_DNA_LEN)
    }

    /// Disable the computation of DNA.
    #[inline(always)]
    pub const fn ignore_dna(self) -> Self {
        self.without_flag(
            COMPUTE_DNA_STRING
                | COMPUTE_DNA_COLUMNAR
                | COMPUTE_DNA_PACKED
                | SPLIT_NON_ACTG
                | REJECT_NON_ACTG
                | RETURN_DNA_CHUNK,
        )
        .return_record(true)
    }

    /// Set the DNA format to bytes (default).
    #[inline(always)]
    pub const fn dna_string(self) -> Self {
        self.ignore_dna().and_dna_string()
    }

    /// Set the DNA format to [`PackedDNA`](crate::dna_format::PackedDNA).
    ///
    /// By default, this option splits the sequence at non-ACTG bases, yielding one [`Event::DnaChunk`](crate::parser::Event) per contiguous ACTG run.
    /// Other behaviors can be set using [`keep_non_actg`](ParserOptions::keep_non_actg), [`skip_non_actg`](ParserOptions::skip_non_actg)
    /// or [`reject_non_actg`](ParserOptions::reject_non_actg).
    ///
    /// Note that the default [`split_non_actg`](ParserOptions::split_non_actg) behavior disables [`Event::Record`](crate::parser::Event) to avoid duplicate processing.
    /// This can be turned back on using [`return_record`](ParserOptions::return_record).
    #[inline(always)]
    pub const fn dna_packed(self) -> Self {
        self.ignore_dna().and_dna_packed()
    }

    /// Set the DNA format to [`ColumnarDNA`](crate::dna_format::ColumnarDNA).
    ///
    /// By default, this option splits the sequence at non-ACTG bases, yielding one [`Event::DnaChunk`](crate::parser::Event) per contiguous ACTG run.
    /// Other behaviors can be set using [`keep_non_actg`](ParserOptions::keep_non_actg), [`skip_non_actg`](ParserOptions::skip_non_actg)
    /// or [`reject_non_actg`](ParserOptions::reject_non_actg).
    ///
    /// Note that the default [`split_non_actg`](ParserOptions::split_non_actg) behavior disables [`Event::Record`](crate::parser::Event) to avoid duplicate processing.
    /// This can be turned back on using [`return_record`](ParserOptions::return_record).
    #[inline(always)]
    pub const fn dna_columnar(self) -> Self {
        self.ignore_dna().and_dna_columnar()
    }

    /// Also compute DNA as bytes.
    #[inline(always)]
    pub const fn and_dna_string(self) -> Self {
        self.with_flag(COMPUTE_DNA_STRING)
    }

    /// Also compute DNA as [`PackedDNA`](crate::dna_format::PackedDNA).
    ///
    /// This calls [`split_non_actg`](ParserOptions::split_non_actg) by default,
    /// which disables [`Event::Record`](crate::parser::Event),
    /// even if [`dna_string`](ParserOptions::dna_string) or [`and_dna_string`](ParserOptions::and_dna_string) was used.
    #[inline(always)]
    pub const fn and_dna_packed(self) -> Self {
        self.with_flag(COMPUTE_DNA_PACKED).split_non_actg()
    }

    /// Also compute DNA as [`ColumnarDNA`](crate::dna_format::ColumnarDNA).
    ///
    /// This calls [`split_non_actg`](ParserOptions::split_non_actg) by default,
    /// which disables [`Event::Record`](crate::parser::Event),
    /// even if [`dna_string`](ParserOptions::dna_string) or [`and_dna_string`](ParserOptions::and_dna_string) was used.
    #[inline(always)]
    pub const fn and_dna_columnar(self) -> Self {
        self.with_flag(COMPUTE_DNA_COLUMNAR).split_non_actg()
    }

    /// Compute a [`BitMask`](crate::dna_format::BitMask) indicating non-ACTG bases.
    /// This is only relevant when using [`keep_non_actg`](ParserOptions::keep_non_actg).
    #[inline(always)]
    pub const fn compute_mask_non_actg(self) -> Self {
        self.with_flag(COMPUTE_MASK_NON_ACTG)
    }

    /// Keep the non-ACTG bases in the sequence, even if their encoding is lossy.
    ///
    /// This is the default behaviour with [`dna_string`](ParserOptions::dna_string).
    /// Combine with [`compute_mask_non_actg`](ParserOptions::compute_mask_non_actg) to
    /// identify which positions were non-ACTG.
    ///
    /// # Example
    ///
    /// ```rust
    /// use helicase::*;
    ///
    /// const CONFIG: Config = ParserOptions::default()
    ///     .dna_string()
    ///     .keep_non_actg()
    ///     .compute_mask_non_actg()
    ///     .config();
    /// ```
    #[inline(always)]
    pub const fn keep_non_actg(self) -> Self {
        self.without_flag(SPLIT_NON_ACTG | RETURN_DNA_CHUNK | MERGE_DNA_CHUNKS | REJECT_NON_ACTG)
            .return_record(true)
    }

    /// Split the sequence at non-ACTG bases, yielding one [`Event::DnaChunk`](crate::parser::Event) per contiguous ACTG run.
    ///
    /// This is the default behaviour with [`dna_packed`](ParserOptions::dna_packed) and [`dna_columnar`](ParserOptions::dna_columnar).
    /// This also disables [`Event::Record`](crate::parser::Event) to avoid processing each sequence twice.
    /// Use [`return_record`](ParserOptions::return_record) to get `Record` events back.
    #[inline(always)]
    pub const fn split_non_actg(self) -> Self {
        self.without_flag(MERGE_DNA_CHUNKS | REJECT_NON_ACTG)
            .with_flag(SPLIT_NON_ACTG | RETURN_DNA_CHUNK)
            .return_record(false)
    }

    /// Skip non-ACTG bases, merging the remaining ACTG runs into a single chunk per record.
    #[inline(always)]
    pub const fn skip_non_actg(self) -> Self {
        self.without_flag(RETURN_DNA_CHUNK | REJECT_NON_ACTG)
            .with_flag(SPLIT_NON_ACTG | MERGE_DNA_CHUNKS)
            .return_record(true)
    }

    /// Return an [`ErrorKind::InvalidBase`](crate::parser::ErrorKind::InvalidBase) at the first non-ACTG base,
    /// yielding one [`Event::Record`](crate::parser::Event) per valid record.
    #[inline(always)]
    pub const fn reject_non_actg(self) -> Self {
        self.without_flag(RETURN_DNA_CHUNK | MERGE_DNA_CHUNKS)
            .with_flag(SPLIT_NON_ACTG | REJECT_NON_ACTG)
            .return_record(true)
    }

    /// Append the DNA of each record (or chunk) to the previous ones instead of replacing it, mostly useful with bitpacked formats.
    /// Header and quality stay per record.
    ///
    /// Memory grows with the input until [`clear_record`](crate::HelicaseParser::clear_record),
    /// `get_*_owned` or `swap_*` is called.
    #[inline(always)]
    pub const fn merge_records(self) -> Self {
        self.with_flag(MERGE_RECORDS)
    }

    /// Stop the parser iterator after each record (`true` by default).
    ///
    /// Set to `false` by default after [`split_non_actg`](ParserOptions::split_non_actg)
    /// (and anything that calls it, such as [`dna_packed`](ParserOptions::dna_packed) or [`and_dna_columnar`](ParserOptions::and_dna_columnar))
    /// to avoid processing each sequence twice.
    #[inline(always)]
    pub const fn return_record(self, enable: bool) -> Self {
        if enable {
            self.with_flag(RETURN_RECORD)
        } else {
            self.without_flag(RETURN_RECORD)
        }
    }

    /// Stop the parser iterator after each DNA chunk.
    ///
    /// Disabled by default with [`dna_string`](ParserOptions::dna_string),
    /// enabled by default with [`dna_packed`](ParserOptions::dna_packed) and [`dna_columnar`](ParserOptions::dna_columnar).
    #[inline(always)]
    pub const fn return_dna_chunk(self, enable: bool) -> Self {
        if enable {
            self.with_flag(RETURN_DNA_CHUNK)
        } else {
            self.without_flag(RETURN_DNA_CHUNK)
        }
    }
}
