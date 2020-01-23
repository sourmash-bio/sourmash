use failure::Error;
use serde_derive::{Deserialize, Serialize};

use crate::signature::SigsTrait;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlatUKHS {}

impl FlatUKHS {
    pub fn md5sum(&self) -> String {
        unimplemented!()
    }
}

impl SigsTrait for FlatUKHS {
    fn size(&self) -> usize {
        unimplemented!()
    }

    fn to_vec(&self) -> Vec<u64> {
        unimplemented!()
    }

    fn ksize(&self) -> usize {
        unimplemented!()
    }

    fn check_compatible(&self, _other: &Self) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, _seq: &[u8], _force: bool) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_protein(&mut self, _seq: &[u8]) -> Result<(), Error> {
        unimplemented!()
    }
}
