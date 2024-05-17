use crate::signature::SigsTrait;

#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct CBL {
    inner: cbl::CBL<21, u64>,
}

impl std::fmt::Debug for CBL {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CBL<21, u64>")
    }
}

// FIXME: don't do this
unsafe impl std::marker::Send for CBL {}
// FIXME: don't do this either
unsafe impl std::marker::Sync for CBL {}

impl SigsTrait for CBL {
    fn size(&self) -> usize {
        unimplemented!()
    }

    fn seed(&self) -> u64 {
        unimplemented!()
    }

    fn ksize(&self) -> usize {
        unimplemented!()
    }

    fn to_vec(&self) -> Vec<u64> {
        unimplemented!()
    }

    fn add_hash(&mut self, hash: crate::HashIntoType) {
        unimplemented!()
    }

    fn add_protein(&mut self, seq: &[u8]) -> Result<(), crate::Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), crate::Error> {
        unimplemented!()
    }

    fn hash_function(&self) -> crate::encodings::HashFunctions {
        unimplemented!()
    }

    fn check_compatible(&self, other: &Self) -> Result<(), crate::Error> {
        unimplemented!()
    }
}
