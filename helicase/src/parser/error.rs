//! Parser error types and implementations.

use alloc::boxed::Box;
use alloc::string::String;
use core::fmt;

/// An error encountered while parsing FASTA/FASTQ data or reading inputs.
///
/// Guaranteed pointer-sized (8 bytes) via boxed representation to ensure that
/// `Result<Event, Error>` does not spill onto the stack in the hot parsing loop.
#[derive(Debug)]
pub struct Error {
    repr: Box<ErrorKind>,
}

/// The specific category of parser error.
#[derive(Debug)]
pub enum ErrorKind {
    /// The input data was empty.
    EmptyInput,

    /// The first record does not begin with an expected format marker
    /// (`>` for FASTA, `@` for FASTQ, or either for auto-detection).
    InvalidRecordStart { expected: &'static str, found: u8 },

    /// The input stream or buffer ended unexpectedly in the middle of a record.
    TruncatedRecord,

    /// Decompression error occurred while reading a compressed stream (gzip/BGZF/zstd).
    Decompression(String),

    /// An underlying I/O error occurred while reading from a file or stream.
    Io(std::io::Error),

    /// Generic or custom error message.
    Custom(String),
}

impl Error {
    /// Returns a reference to the inner [`ErrorKind`].
    #[inline(always)]
    pub fn kind(&self) -> &ErrorKind {
        &self.repr
    }

    /// Creates a new `Error` wrapping the specified [`ErrorKind`].
    #[cold]
    #[inline(never)]
    pub fn new(kind: ErrorKind) -> Self {
        Self {
            repr: Box::new(kind),
        }
    }

    /// Constructs an [`ErrorKind::EmptyInput`] error.
    #[cold]
    #[inline(never)]
    pub fn empty_input() -> Self {
        Self::new(ErrorKind::EmptyInput)
    }

    /// Constructs an [`ErrorKind::InvalidRecordStart`] error.
    #[cold]
    #[inline(never)]
    pub fn invalid_record_start(expected: &'static str, found: u8) -> Self {
        Self::new(ErrorKind::InvalidRecordStart { expected, found })
    }

    /// Constructs an [`ErrorKind::TruncatedRecord`] error.
    #[cold]
    #[inline(never)]
    pub fn truncated_record() -> Self {
        Self::new(ErrorKind::TruncatedRecord)
    }

    /// Constructs an [`ErrorKind::Decompression`] error.
    #[cold]
    #[inline(never)]
    pub fn decompression(msg: impl Into<String>) -> Self {
        Self::new(ErrorKind::Decompression(msg.into()))
    }

    /// Constructs an [`ErrorKind::Custom`] error.
    #[cold]
    #[inline(never)]
    pub fn custom(msg: impl Into<String>) -> Self {
        Self::new(ErrorKind::Custom(msg.into()))
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "input data is empty"),
            Self::InvalidRecordStart { expected, found } => {
                write!(
                    f,
                    "invalid record start: expected {expected}, found byte {found:#04x} ('{}')",
                    *found as char
                )
            }
            Self::TruncatedRecord => write!(f, "unexpected end of file mid-record"),
            Self::Decompression(msg) => write!(f, "decompression error: {msg}"),
            Self::Io(err) => write!(f, "I/O error: {err}"),
            Self::Custom(msg) => write!(f, "{msg}"),
        }
    }
}

impl fmt::Display for Error {
    #[inline(always)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&*self.repr, f)
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &*self.repr {
            ErrorKind::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    #[cold]
    fn from(err: std::io::Error) -> Self {
        if let Some(inner) = err.get_ref()
            && inner.is::<Self>()
        {
            if let Some(e) = err.into_inner()
                && let Ok(downcasted) = e.downcast::<Self>()
            {
                return *downcasted;
            }
            return Self::custom("failed to downcast Error from std::io::Error");
        }
        Self::new(ErrorKind::Io(err))
    }
}

impl From<Error> for std::io::Error {
    #[cold]
    fn from(err: Error) -> Self {
        match *err.repr {
            ErrorKind::Io(e) => e,
            _ => Self::new(std::io::ErrorKind::InvalidData, err),
        }
    }
}

impl From<&str> for Error {
    #[cold]
    fn from(s: &str) -> Self {
        Self::custom(s)
    }
}

impl From<String> for Error {
    #[cold]
    fn from(s: String) -> Self {
        Self::custom(s)
    }
}
