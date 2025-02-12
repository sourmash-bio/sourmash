use needletail::parse_fastx_reader;
use wasm_bindgen::prelude::*;

use crate::cmd::ComputeParameters as _ComputeParameters;
use crate::encodings::HashFunctions;
use crate::prelude::ToWriter;
use crate::signature::Signature as _Signature;
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash as _KmerMinHash;
use crate::ScaledType;

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
        scaled: ScaledType,
        track_abundance: bool,
    ) -> KmerMinHash {
        // TODO: at most one of (prot, dayhoff, hp) should be true

        let hash_function = if dayhoff {
            HashFunctions::Murmur64Dayhoff
        } else if hp {
            HashFunctions::Murmur64Hp
        } else if is_protein {
            HashFunctions::Murmur64Protein
        } else {
            HashFunctions::Murmur64Dna
        };

        KmerMinHash(_KmerMinHash::new(
            scaled,
            ksize,
            hash_function,
            seed as u64,
            track_abundance,
            num,
        ))
    }

    #[wasm_bindgen]
    pub fn add_sequence_js(&mut self, buf: &str) -> Result<(), JsErrors> {
        self.0.add_sequence(buf.as_bytes(), true)?;
        Ok(())
    }

    #[wasm_bindgen]
    pub fn to_json(&mut self) -> Result<String, JsErrors> {
        let mut st: Vec<u8> = vec![];
        self.0.to_writer(&mut st)?;
        Ok(unsafe { String::from_utf8_unchecked(st) })
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

    #[wasm_bindgen]
    pub fn set_scaled(&mut self, scaled: ScaledType) {
        self.0.set_scaled(scaled);
    }

    #[wasm_bindgen]
    pub fn set_num(&mut self, num: u32) {
        self.0.set_num_hashes(num);
    }

    #[wasm_bindgen]
    pub fn set_protein(&mut self, is_protein: bool) {
        self.0.set_protein(is_protein);
    }

    #[wasm_bindgen]
    pub fn set_dayhoff(&mut self, dayhoff: bool) {
        self.0.set_dayhoff(dayhoff);
    }

    #[wasm_bindgen]
    pub fn set_hp(&mut self, hp: bool) {
        self.0.set_hp(hp);
    }

    #[wasm_bindgen]
    pub fn set_track_abundance(&mut self, track: bool) {
        self.0.set_track_abundance(track);
    }
    #[wasm_bindgen]
    pub fn set_seed(&mut self, seed: u32) {
        self.0.set_seed(seed.into());
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
    pub fn add_sequence_js(&mut self, buf: &str) -> Result<(), JsErrors> {
        self.0.add_sequence(buf.as_bytes(), true)?;

        Ok(())
    }

    #[wasm_bindgen]
    pub fn add_from_file(
        &mut self,
        fp: web_sys::File,
        callback: Option<js_sys::Function>,
    ) -> Result<(), JsErrors> {
        let wf = SyncFile::new(fp, callback);

        let (rdr, _format) = niffler::send::get_reader(Box::new(wf))?;

        let mut parser = parse_fastx_reader(std::io::BufReader::with_capacity(
            1024 << 14, // 16 MiB
            rdr,
        ))?;

        while let Some(record) = parser.next() {
            let record = record?;
            self.0.add_sequence(&record.seq(), true)?;
        }

        Ok(())
    }

    #[wasm_bindgen]
    pub fn to_json(&mut self) -> Result<String, JsErrors> {
        let mut st: Vec<u8> = vec![];
        self.0.to_writer(&mut st)?;
        Ok(unsafe { String::from_utf8_unchecked(st) })
    }

    pub fn size(&self) -> usize {
        self.0.size()
    }
}

#[derive(thiserror::Error, Debug)]
pub enum JsErrors {
    #[error(transparent)]
    SourmashError(#[from] crate::Error),

    #[error(transparent)]
    SerdeError(#[from] serde_json::error::Error),

    #[error(transparent)]
    NifflerError(#[from] niffler::Error),

    #[error(transparent)]
    NeedletailError(#[from] needletail::errors::ParseError),
}

impl Into<JsValue> for JsErrors {
    fn into(self) -> JsValue {
        let error = js_sys::Error::new(&self.to_string());
        error.into()
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

// ==============================

use js_sys::Number;
use js_sys::Uint8Array;
use once_cell::sync::Lazy;
use web_sys::FileReaderSync;

thread_local! {
    static FILE_READER_SYNC: Lazy<FileReaderSync> = Lazy::new(|| {
      FileReaderSync::new().expect("Failed to create FileReaderSync. Is it running in a web worker context?")
    });
}

/// Wrapper around a `web_sys::File` that implements `Read` and `Seek`.
pub struct SyncFile {
    file: web_sys::File,
    pos: u64,
    cb: Option<js_sys::Function>,
}

/// Because this needs to be initialized in a Web Worker, it is safe to make it Send.
/// (hopefully. I don't think they can be sent across Web Workers, nor accessed from other WW)
unsafe impl Send for SyncFile {}

impl SyncFile {
    pub fn new(file: web_sys::File, cb: Option<js_sys::Function>) -> Self {
        Self { file, pos: 0, cb }
    }

    /// File size in bytes.
    pub fn size(&self) -> u64 {
        let size = self.file.size();
        if size <= Number::MAX_SAFE_INTEGER {
            return size as u64;
        } else {
            panic!("size is not safe to convert to integer from float")
        }
    }

    fn set_pos(&mut self, pos: u64) {
        self.pos = pos;
        self.cb.as_ref().map(|f| {
            let arr = js_sys::Array::new_with_length(1);
            arr.set(0, self.progress().into());
            f.apply(&JsValue::null(), &arr)
                .expect("Error calling progress callback");
        });
    }

    /// Current progress on the file
    pub fn progress(&self) -> f64 {
        self.pos as f64 / self.file.size()
    }
}

impl std::io::Read for SyncFile {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, std::io::Error> {
        let current_offset = self.pos;
        let new_offset_f64 = current_offset as f64;
        let new_offset_end_f64 = current_offset.saturating_add(
            u64::try_from(buf.len()).map_err(|_| std::io::Error::other("Can't convert to u64"))?,
        ) as f64;

        let blob = self
            .file
            .slice_with_f64_and_f64(new_offset_f64, new_offset_end_f64)
            .map_err(|_| std::io::Error::other("failed to slice file"))?;
        let array_buffer = FILE_READER_SYNC
            .with(|frs| frs.read_as_array_buffer(&blob))
            .map_err(|_| std::io::Error::other("failed to read as array buffer"))?;

        let array = Uint8Array::new(&array_buffer);
        let read_bytes = usize::try_from(array.byte_length())
            .map_err(|_| std::io::Error::other("read too many bytes at once"))?;

        // Copy to output buffer
        array.copy_to(&mut buf[..read_bytes]);

        // Update position
        self.set_pos(
            current_offset
                .checked_add(read_bytes as u64)
                .ok_or_else(|| std::io::Error::other("new position too large"))?,
        );

        Ok(read_bytes)
    }
}
