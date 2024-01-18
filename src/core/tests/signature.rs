use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::signature::Signature;
use sourmash::sketch::Sketch;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[test]
fn selection_with_downsample() {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/47+63-multisig.sig");

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

    // create Selection object
    let mut selection = Selection::default();
    selection.set_scaled(2000);
    // iterate and check scaled
    for sig in &sigs {
        let modified_sig = sig.clone().select(&selection).unwrap();
        for sketch in modified_sig.sketches() {
            if let Sketch::MinHash(mh) = sketch {
                eprintln!("scaled: {:?}", mh.scaled());
                assert_eq!(mh.scaled(), 2000);
            }
        }
    }
}
