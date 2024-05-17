pub mod hyperloglog;
pub mod minhash;

pub mod nodegraph;

pub mod cbl;

use serde::{Deserialize, Serialize};

use crate::sketch::cbl::CBL;
use crate::sketch::hyperloglog::HyperLogLog;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Serialize, rkyv::Deserialize, rkyv::Archive)
)]
pub enum Sketch {
    MinHash(KmerMinHash),
    LargeMinHash(KmerMinHashBTree),
    HyperLogLog(HyperLogLog),
    LowScaled(CBL),
}
