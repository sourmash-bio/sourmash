use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::result::Result;

use mqf::MQF;

use sourmash::signature::Signature;
use sourmash::sketch::Sketch;

fn main() {
    let mh_paths: HashMap<u8, String> = {
        [
            (6, "6d6e87e1154e95b279e5e7db414bc37b".into()),
            (7, "60f7e23c24a8d94791cc7a8680c493f9".into()),
            (8, "0107d767a345eff67ecdaed2ee5cd7ba".into()),
            (9, "f71e78178af9e45e6f1d87a0c53c465c".into()),
            (10, "f0c834bc306651d2b9321fb21d3e8d8f".into()),
            (11, "4e94e60265e04f0763142e20b52c0da1".into()),
            (12, "b59473c94ff2889eca5d7165936e64b3".into()),
        ]
        .into_iter()
        .cloned()
        .collect()
    };

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/.sbt.v5_mhmt/");

    let mh_sigs: HashMap<u8, Signature> = mh_paths
        .into_iter()
        .map(|(k, v)| {
            let mut mhpath = filename.clone();
            mhpath.push(v);

            let file = File::open(mhpath).unwrap();
            let reader = BufReader::new(file);
            let sig: Result<Vec<Signature>, _> = serde_json::from_reader(reader);
            (k, sig.unwrap()[0].clone())
        })
        .collect();

    let mut mqfs: HashMap<u8, MQF> = HashMap::default();

    for i in (0..=5).rev() {
        println!("Creating MQF {}", i);
        let mut mqf = MQF::new(1, 18);

        let left = i * 2 + 1;
        for child in left..=left + 1 {
            if mh_sigs.contains_key(&child) {
                println!("Loading values from MH {}", child);
                if let Sketch::MinHash(sig) = &mh_sigs[&child].signatures[0] {
                    sig.mins()
                        .iter()
                        .map(|h| {
                            dbg!(*h % u64::pow(2, 26));
                            mqf.insert(*h % u64::pow(2, 26), 1)
                            //mqf.insert(*h, 1)
                        })
                        .count();
                };
            } else if mqfs.contains_key(&child) {
                let mut cmqf = mqfs.get_mut(&child).unwrap();
                mqf.merge(&mut cmqf).expect("Error merging");
            } else {
                // TODO: shouldn't happen...
                unimplemented!()
            }
        }

        println!("Save MQF {} for later", i);
        mqfs.insert(i, mqf);

        println!("Saving MQFs to disk");
        let mut internal = filename.clone();
        internal.push(format!("internal.{}", i));
        mqfs[&i].serialize(internal).unwrap();
    }
}
