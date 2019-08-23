use failure::{Error, Fail};

#[derive(Debug, Fail)]
pub enum SourmashError {
    /// Raised for internal errors in the libraries.  Should not happen.
    #[fail(display = "internal error: {}", message)]
    Internal { message: String },

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

    #[fail(display = "invalid DNA character in input k-mer: {}", message)]
    InvalidDNA { message: String },

    #[fail(display = "invalid protein character in input: {}", message)]
    InvalidProt { message: String },

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
    // Input sequence errors
    InvalidDNA = 11_01,
    InvalidProt = 11_02,
    // external errors
    Io = 100_001,
    Utf8Error = 100_002,
    ParseInt = 100_003,
    SerdeError = 100_004,
}

impl SourmashErrorCode {
    pub fn from_error(error: &Error) -> SourmashErrorCode {
        for cause in error.iter_chain() {
            use crate::utils::Panic;
            if cause.downcast_ref::<Panic>().is_some() {
                return SourmashErrorCode::Panic;
            }

            if let Some(err) = cause.downcast_ref::<SourmashError>() {
                return match err {
                    SourmashError::Internal { .. } => SourmashErrorCode::Internal,
                    SourmashError::MismatchKSizes => SourmashErrorCode::MismatchKSizes,
                    SourmashError::MismatchDNAProt => SourmashErrorCode::MismatchDNAProt,
                    SourmashError::MismatchMaxHash => SourmashErrorCode::MismatchMaxHash,
                    SourmashError::MismatchSeed => SourmashErrorCode::MismatchSeed,
                    SourmashError::MismatchSignatureType => {
                        SourmashErrorCode::MismatchSignatureType
                    }
                    SourmashError::InvalidDNA { .. } => SourmashErrorCode::InvalidDNA,
                    SourmashError::InvalidProt { .. } => SourmashErrorCode::InvalidProt,
                    SourmashError::SerdeError => SourmashErrorCode::SerdeError,
                };
            }
        }
        SourmashErrorCode::Unknown
    }
}
