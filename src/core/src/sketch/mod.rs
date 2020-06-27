pub mod minhash;
pub mod nodegraph;

pub mod ukhs;

use serde_derive::{Deserialize, Serialize};

use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::sketch::ukhs::FlatUKHS;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Sketch {
    MinHash(KmerMinHash),
    LargeMinHash(KmerMinHashBTree),
    UKHS(FlatUKHS), // FIXME
}
