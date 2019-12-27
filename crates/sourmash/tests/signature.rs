use serde_json;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::signature::Signature;

#[test]
fn load_signature() {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/genome-s10+s11.sig");

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

    assert_eq!(sigs.len(), 1);

    let sig = sigs.get(0).unwrap();
    assert_eq!(sig.class, "sourmash_signature");
    assert_eq!(sig.email, "");
    if let Some(ref filename) = sig.filename {
        assert_eq!(filename, "-");
    }
    assert_eq!(sig.hash_function, "0.murmur64");
    if let Some(ref name) = sig.name {
        assert_eq!(name, "s10+s11");
    }
    assert_eq!(sig.signatures.len(), 4);
}
