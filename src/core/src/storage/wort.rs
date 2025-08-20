use crate::storage::{Storage, StorageArgs};
use crate::Result;

/// Load data from wort (https://wort.sourmash.bio)
///
/// This is read-only, no support for writing data to wort.
#[derive(Debug, Clone)]
pub struct WortStorage {
    // TODO: save reqwest blocking client here?
    client: reqwest::blocking::Client,

    // TODO: use a default URL pointing to wort.sourmash.bio, and allow mirrors?
}

impl WortStorage {
    pub fn new() -> Self {
        Self {
          client: reqwest::blocking::Client::new(),
        }
    }
}

impl Storage for WortStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String> {
        unimplemented!()
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let resp = self.client.get(format!("/{}", path)).send()?;
        Ok(resp.bytes()?.into())
    }

    fn args(&self) -> StorageArgs {
        unimplemented!()
    }

    fn spec(&self) -> String {
        "wort://".into()
    }
}
