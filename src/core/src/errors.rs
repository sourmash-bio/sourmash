use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SourmashError {
    /// Raised for internal errors in the libraries.  Should not happen.
    #[error("internal error: {message:?}")]
    Internal { message: String },

    #[error("new scaled smaller than previous; cannot upsample")]
    CannotUpsampleScaled,

    #[error("must have same num: {n1} != {n2}")]
    MismatchNum { n1: u32, n2: u32 },

    #[error("different ksizes cannot be compared")]
    MismatchKSizes,

    #[error("DNA/prot minhashes cannot be compared")]
    MismatchDNAProt,

    #[error("mismatch in scaled; comparison fail")]
    MismatchScaled,

    #[error("mismatch in seed; comparison fail")]
    MismatchSeed,

    #[error("different signatures cannot be compared")]
    MismatchSignatureType,

    #[error("sketch needs abundance for this operation")]
    NeedsAbundanceTracking,

    #[error("Expected a MinHash sketch in this signature")]
    NoMinHashFound,

    #[error("Empty signature")]
    EmptySignature,

    #[error("Multiple sketches found, expected one")]
    MultipleSketchesFound,

    #[error("Invalid hash function: {function:?}")]
    InvalidHashFunction { function: String },

    #[error("Can only set {message:?} if the MinHash is empty")]
    NonEmptyMinHash { message: String },

    #[error("invalid DNA character in input k-mer: {message}")]
    InvalidDNA { message: String },

    #[error("invalid protein character in input: {message}")]
    InvalidProt { message: String },

    #[error("Codon is invalid length: {message}")]
    InvalidCodonLength { message: String },

    #[error("Skipmer ksize must be >= n ({n}), but got ksize: {ksize}")]
    InvalidSkipmerSize { ksize: usize, n: usize },

    #[error("Skipmer frame number must be < n ({n}), but got start: {start}")]
    InvalidSkipmerFrame { start: usize, n: usize },

    #[error("Frame number must be 0, 1, or 2, but got {frame_number}")]
    InvalidTranslateFrame { frame_number: usize },

    #[error("Set error rate to a value smaller than 0.367696 and larger than 0.00203125")]
    HLLPrecisionBounds,

    #[error("error while calculating ANI confidence intervals: {message}")]
    ANIEstimationError { message: String },

    #[error(transparent)]
    ReadDataError(#[from] ReadDataError),

    #[error(transparent)]
    StorageError(#[from] crate::storage::StorageError),

    #[error(transparent)]
    SerdeError(#[from] serde_json::error::Error),

    #[error(transparent)]
    NifflerError(#[from] niffler::Error),

    #[error(transparent)]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error(transparent)]
    IOError(#[from] std::io::Error),

    #[error(transparent)]
    CsvError(#[from] csv::Error),

    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    #[error(transparent)]
    Panic(#[from] crate::ffi::utils::Panic),

    #[cfg(not(target_arch = "wasm32"))]
    #[cfg(feature = "branchwater")]
    #[error(transparent)]
    RocksDBError(#[from] rocksdb::Error),

    #[error(transparent)]
    ZipError(#[from] piz::result::ZipError),
}

#[derive(Debug, Error)]
pub enum ReadDataError {
    #[error("Could not load data")]
    LoadError,
}

#[repr(u32)]
pub enum SourmashErrorCode {
    // no error
    NoError = 0,
    // panics and internals
    Panic = 1,
    Internal = 2,
    Msg = 3,
    Unknown = 4,
    // Compatibility errors
    MismatchKSizes = 1_01,
    MismatchDNAProt = 1_02,
    MismatchScaled = 1_03,
    MismatchSeed = 1_04,
    MismatchSignatureType = 1_05,
    NonEmptyMinHash = 1_06,
    MismatchNum = 1_07,
    NeedsAbundanceTracking = 1_08,
    CannotUpsampleScaled = 1_09,
    NoMinHashFound = 1_10,
    EmptySignature = 1_11,
    MultipleSketchesFound = 1_12,
    // Input sequence errors
    InvalidDNA = 11_01,
    InvalidProt = 11_02,
    InvalidCodonLength = 11_03,
    InvalidHashFunction = 11_04,
    InvalidSkipmerFrame = 11_05,
    InvalidSkipmerSize = 11_06,
    InvalidTranslateFrame = 11_07,
    // index-related errors
    ReadData = 12_01,
    Storage = 12_02,
    // HLL errors
    HLLPrecisionBounds = 13_01,
    // ANI errors
    ANIEstimationError = 14_01,
    // external errors
    Io = 100_001,
    Utf8Error = 100_002,
    ParseInt = 100_003,
    SerdeError = 100_004,
    NifflerError = 100_005,
    CsvError = 100_006,
    RocksDBError = 100_007,
    ZipError = 100_008,
}

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
impl SourmashErrorCode {
    pub fn from_error(error: &SourmashError) -> SourmashErrorCode {
        match error {
            SourmashError::Internal { .. } => SourmashErrorCode::Internal,
            SourmashError::Panic { .. } => SourmashErrorCode::Panic,
            SourmashError::CannotUpsampleScaled => SourmashErrorCode::CannotUpsampleScaled,
            SourmashError::MismatchNum { .. } => SourmashErrorCode::MismatchNum,
            SourmashError::NeedsAbundanceTracking => SourmashErrorCode::NeedsAbundanceTracking,
            SourmashError::MismatchKSizes => SourmashErrorCode::MismatchKSizes,
            SourmashError::MismatchDNAProt => SourmashErrorCode::MismatchDNAProt,
            SourmashError::MismatchScaled => SourmashErrorCode::MismatchScaled,
            SourmashError::MismatchSeed => SourmashErrorCode::MismatchSeed,
            SourmashError::MismatchSignatureType => SourmashErrorCode::MismatchSignatureType,
            SourmashError::NonEmptyMinHash { .. } => SourmashErrorCode::NonEmptyMinHash,
            SourmashError::NoMinHashFound => SourmashErrorCode::NoMinHashFound,
            SourmashError::EmptySignature => SourmashErrorCode::EmptySignature,
            SourmashError::MultipleSketchesFound => SourmashErrorCode::MultipleSketchesFound,
            SourmashError::InvalidDNA { .. } => SourmashErrorCode::InvalidDNA,
            SourmashError::InvalidProt { .. } => SourmashErrorCode::InvalidProt,
            SourmashError::InvalidCodonLength { .. } => SourmashErrorCode::InvalidCodonLength,
            SourmashError::InvalidHashFunction { .. } => SourmashErrorCode::InvalidHashFunction,
            SourmashError::InvalidSkipmerFrame { .. } => SourmashErrorCode::InvalidSkipmerFrame,
            SourmashError::InvalidSkipmerSize { .. } => SourmashErrorCode::InvalidSkipmerSize,
            SourmashError::InvalidTranslateFrame { .. } => SourmashErrorCode::InvalidTranslateFrame,
            SourmashError::ReadDataError { .. } => SourmashErrorCode::ReadData,
            SourmashError::StorageError { .. } => SourmashErrorCode::Storage,
            SourmashError::HLLPrecisionBounds => SourmashErrorCode::HLLPrecisionBounds,
            SourmashError::ANIEstimationError { .. } => SourmashErrorCode::ANIEstimationError,
            SourmashError::SerdeError { .. } => SourmashErrorCode::SerdeError,
            SourmashError::IOError { .. } => SourmashErrorCode::Io,
            SourmashError::NifflerError { .. } => SourmashErrorCode::NifflerError,
            SourmashError::Utf8Error { .. } => SourmashErrorCode::Utf8Error,
            SourmashError::CsvError { .. } => SourmashErrorCode::CsvError,

            #[cfg(not(target_arch = "wasm32"))]
            #[cfg(feature = "branchwater")]
            SourmashError::RocksDBError { .. } => SourmashErrorCode::RocksDBError,

            SourmashError::ZipError { .. } => SourmashErrorCode::ZipError,
        }
    }
}
