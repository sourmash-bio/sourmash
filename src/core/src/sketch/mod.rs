pub mod draff;
pub mod hyperloglog;
pub mod minhash;
pub mod nodegraph;

use serde::{Deserialize, Serialize};

use crate::sketch::draff::FlatUKHS as Flat;
use crate::sketch::hyperloglog::HyperLogLog;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
#[non_exhaustive]
pub enum Sketch {
    MinHash(KmerMinHash),
    LargeMinHash(KmerMinHashBTree),
    Draff(Flat),
    HyperLogLog(HyperLogLog),
}
