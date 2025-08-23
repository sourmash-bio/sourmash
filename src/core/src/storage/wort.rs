use crate::storage::{Storage, StorageArgs};
use crate::Result;

#[cfg(not(target_arch = "wasm32"))]
use reqwest::blocking::Client;

#[cfg(target_arch = "wasm32")]
use reqwest::Client;

/// Load data from wort (https://wort.sourmash.bio)
///
/// This is read-only, no support for writing data to wort.
#[derive(Debug, Clone)]
pub struct WortStorage {
    // Save a reqwest client here to avoid initialization on every download
    client: Client,

    // Base URL for the wort API, by default https://wort.sourmash.bio/v1/view
    base_url: String,
}

impl Default for WortStorage {
    fn default() -> Self {
        Self::new()
    }
}

impl WortStorage {
    pub fn new() -> Self {
        Self {
            client: Client::new(),
            base_url: "https://wort.sourmash.bio/v1/view".to_string(),
        }
    }
}

impl Storage for WortStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String> {
        unimplemented!()
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let resp = self
            .client
            .get(format!("{}/{}", self.base_url, path))
            .send();

        #[cfg(target_arch = "wasm32")]
        let data = resp.await?.bytes().await?;

        #[cfg(not(target_arch = "wasm32"))]
        let data = resp?.bytes()?;

        Ok(data.into())
    }

    fn args(&self) -> StorageArgs {
        unimplemented!()
    }

    fn spec(&self) -> String {
        "wort://".into()
    }
}
