// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

use wasm_bindgen::prelude::*;

use crate::cmd::ComputeParameters as _ComputeParameters;
use crate::encodings::HashFunctions;
use crate::signature::Signature as _Signature;
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash as _KmerMinHash;

#[wasm_bindgen]
pub struct KmerMinHash(_KmerMinHash);

#[wasm_bindgen]
pub struct Signature(_Signature);

#[wasm_bindgen]
pub struct ComputeParameters(_ComputeParameters);

#[wasm_bindgen]
impl KmerMinHash {
    #[wasm_bindgen(constructor)]
    pub fn new_with_scaled(
        num: u32,
        ksize: u32,
        is_protein: bool,
        dayhoff: bool,
        hp: bool,
        seed: u32,
        scaled: u32,
        track_abundance: bool,
    ) -> KmerMinHash {
        // TODO: at most one of (prot, dayhoff, hp) should be true

        let hash_function = if dayhoff {
            HashFunctions::murmur64_dayhoff
        } else if hp {
            HashFunctions::murmur64_hp
        } else if is_protein {
            HashFunctions::murmur64_protein
        } else {
            HashFunctions::murmur64_DNA
        };

        KmerMinHash(_KmerMinHash::new(
            scaled as u64,
            ksize,
            hash_function,
            seed as u64,
            track_abundance,
            num,
        ))
    }

    #[wasm_bindgen]
    pub fn add_sequence_js(&mut self, buf: &str) {
        self.0
            .add_sequence(buf.as_bytes(), true)
            .expect("Error adding sequence");
    }

    #[wasm_bindgen]
    pub fn to_json(&mut self) -> String {
        serde_json::to_string(&self.0).unwrap()
    }
}

#[wasm_bindgen]
impl ComputeParameters {
    #[wasm_bindgen(constructor)]
    pub fn new_with_params() -> ComputeParameters {
        let params = _ComputeParameters::default();
        ComputeParameters(params)
    }

    #[wasm_bindgen]
    pub fn set_ksizes(&mut self, ksizes: Vec<u32>) {
        self.0.set_ksizes(ksizes);
    }
}

#[wasm_bindgen]
impl Signature {
    #[wasm_bindgen(constructor)]
    pub fn new_from_params(params: &ComputeParameters) -> Signature {
        //let params = ComputeParameters::default();

        Signature(_Signature::from_params(&params.0))
    }

    #[wasm_bindgen]
    pub fn add_sequence_js(&mut self, buf: &str) {
        self.0
            .add_sequence(buf.as_bytes(), true)
            .expect("Error adding sequence");
    }

    #[wasm_bindgen]
    pub fn add_from_file(&mut self, fp: web_sys::File) {
        unimplemented!()
    }

    #[wasm_bindgen]
    pub fn to_json(&mut self) -> String {
        serde_json::to_string(&self.0).unwrap()
    }

    pub fn size(&self) -> usize {
        self.0.size()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use wasm_bindgen_test::*;

    #[wasm_bindgen_test]
    fn wasm_test() {
        let mut params = ComputeParameters::new_with_params();
        params.set_ksizes(vec![19, 29, 49]);
        let sig = Signature::new_from_params(&params);
        assert_eq!(sig.size(), 3);
    }
}
