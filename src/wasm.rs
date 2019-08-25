use wasm_bindgen::prelude::*;

use serde_json;

use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;

#[wasm_bindgen]
impl KmerMinHash {
    #[wasm_bindgen(constructor)]
    pub fn new_with_scaled(
        num: u32,
        ksize: u32,
        is_protein: bool,
        seed: u32,
        scaled: u32,
        track_abundance: bool,
    ) -> KmerMinHash {
        let max_hash = if num != 0 {
            0
        } else if scaled == 0 {
            u64::max_value()
        } else {
            u64::max_value() / scaled as u64
        };

        KmerMinHash::new(
            num,
            ksize,
            is_protein,
            seed as u64,
            max_hash,
            track_abundance,
        )
    }

    #[wasm_bindgen]
    pub fn add_sequence_js(&mut self, buf: &str) {
        self.add_sequence(buf.as_bytes(), true);
    }

    #[wasm_bindgen]
    pub fn to_json(&mut self) -> String {
        serde_json::to_string(self).unwrap()
    }
}
