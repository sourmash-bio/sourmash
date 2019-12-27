use failure::Fail;

#[derive(Debug, Fail)]
pub enum SourmashError {
    /// Raised for internal errors in the libraries.  Should not happen.
    #[fail(display = "internal error: {}", message)]
    Internal { message: String },

    #[fail(display = "must have same num: {} != {}", n1, n2)]
    MismatchNum { n1: u32, n2: u32 },

    #[fail(display = "different ksizes cannot be compared")]
    MismatchKSizes,

    #[fail(display = "DNA/prot minhashes cannot be compared")]
    MismatchDNAProt,

    #[fail(display = "mismatch in max_hash; comparison fail")]
    MismatchMaxHash,

    #[fail(display = "mismatch in seed; comparison fail")]
    MismatchSeed,

    #[fail(display = "different signatures cannot be compared")]
    MismatchSignatureType,

    #[fail(display = "Can only set {} if the MinHash is empty", message)]
    NonEmptyMinHash { message: String },

    #[fail(display = "invalid DNA character in input k-mer: {}", message)]
    InvalidDNA { message: String },

    #[fail(display = "invalid protein character in input: {}", message)]
    InvalidProt { message: String },

    #[fail(display = "Codon is invalid length: {}", message)]
    InvalidCodonLength { message: String },

    #[fail(display = "Error from deserialization")]
    SerdeError,
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
    MismatchMaxHash = 1_03,
    MismatchSeed = 1_04,
    MismatchSignatureType = 1_05,
    NonEmptyMinHash = 1_06,
    MismatchNum = 1_07,
    // Input sequence errors
    InvalidDNA = 11_01,
    InvalidProt = 11_02,
    InvalidCodonLength = 11_03,
    // external errors
    Io = 100_001,
    Utf8Error = 100_002,
    ParseInt = 100_003,
    SerdeError = 100_004,
}

#[cfg(not(all(target_arch = "wasm32", target_vendor = "unknown")))]
impl SourmashErrorCode {
    pub fn from_error(error: &failure::Error) -> SourmashErrorCode {
        for cause in error.iter_chain() {
            use crate::ffi::utils::Panic;
            if cause.downcast_ref::<Panic>().is_some() {
                return SourmashErrorCode::Panic;
            }

            if let Some(err) = cause.downcast_ref::<SourmashError>() {
                return match err {
                    SourmashError::Internal { .. } => SourmashErrorCode::Internal,
                    SourmashError::MismatchNum { .. } => SourmashErrorCode::MismatchNum,
                    SourmashError::MismatchKSizes => SourmashErrorCode::MismatchKSizes,
                    SourmashError::MismatchDNAProt => SourmashErrorCode::MismatchDNAProt,
                    SourmashError::MismatchMaxHash => SourmashErrorCode::MismatchMaxHash,
                    SourmashError::MismatchSeed => SourmashErrorCode::MismatchSeed,
                    SourmashError::MismatchSignatureType => {
                        SourmashErrorCode::MismatchSignatureType
                    }
                    SourmashError::NonEmptyMinHash { .. } => SourmashErrorCode::NonEmptyMinHash,
                    SourmashError::InvalidDNA { .. } => SourmashErrorCode::InvalidDNA,
                    SourmashError::InvalidProt { .. } => SourmashErrorCode::InvalidProt,
                    SourmashError::InvalidCodonLength { .. } => {
                        SourmashErrorCode::InvalidCodonLength
                    }
                    SourmashError::SerdeError => SourmashErrorCode::SerdeError,
                };
            }
        }
        SourmashErrorCode::Unknown
    }
}
