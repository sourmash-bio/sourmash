use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::index::lca_db::{LcaDB, LineagePair};
use sourmash::sketch::Sketch;
use sourmash::sketch::minhash::{KmerMinHash, HashFunctions};
use std::collections::HashMap;

use proptest::collection::vec;
use proptest::num::u64;
use proptest::proptest;

#[test]
fn build_default_struct() {
    let lca_db = LcaDB::new(31, 1, b"", b"DNA");

    assert!(lca_db.ksize() == 31);
    assert!(lca_db.scaled() == 1);
    assert!(lca_db.filename() == "");
    assert!(lca_db.moltype() == "DNA");
    assert!(lca_db._next_index() as u32 == 0);
    assert!(lca_db._next_lid() as u32 == 0);
    assert!(lca_db.ident_to_name() == HashMap::new());
    assert!(lca_db.ident_to_idx() == HashMap::new());
    assert!(lca_db.idx_to_lid() == HashMap::new());
    assert!(lca_db.lineage_to_lid() == HashMap::new());
    assert!(lca_db.lid_to_lineage() == HashMap::new());
    assert!(lca_db.hashval_to_idx() == HashMap::new());
}

#[test]
fn test_insert() {
    let mut lca_db = LcaDB::new(32, 1, b"", b"DNA");
    let mh = KmerMinHash::new(0, 32, HashFunctions::murmur64_protein, 42, 10000, false);
    let sig = Signature::default();
    let lineage: Vec<LineagePair> = vec![LineagePair::new("name1".to_string(), "rank1".to_string()), 
                                        LineagePair::new("name2".to_string(), "rank2".to_string())];

    unsafe { lca_db.insert(&mh, &sig, b"erik", lineage) };
                                        
    assert!(false);
    assert!(lca_db.ksize() == 32);
    assert!(lca_db.scaled() == 1);
    assert!(lca_db.filename() == "");
    assert!(lca_db.moltype() == "DNA");
    assert!(lca_db._next_index() as u32 == 0);
    assert!(lca_db._next_lid() as u32 == 0);
    assert!(lca_db.ident_to_name() == HashMap::new());
    assert!(lca_db.ident_to_idx() == HashMap::new());
    assert!(lca_db.idx_to_lid() == HashMap::new());
    assert!(lca_db.lineage_to_lid() == HashMap::new());
    assert!(lca_db.lid_to_lineage() == HashMap::new());
    assert!(lca_db.hashval_to_idx() == HashMap::new());
}

#[test]
fn test_save() {
    let mut lca_db2 = LcaDB::new(32, 1, b"", b"DNA");
    let mh = KmerMinHash::new(0, 32, HashFunctions::murmur64_protein, 42, 10000, false);
    let sig = Signature::default();
    let lineage: Vec<LineagePair> = vec![LineagePair::new("name1".to_string(), "rank1".to_string()), 
                                        LineagePair::new("name2".to_string(), "rank2".to_string())];
                                        
    let lca_db = lca_db2.save(&mh, b"eriksjson");

                                        
    assert!(false);
    assert!(lca_db.ksize() == 32);
    assert!(lca_db.scaled() == 1);
    assert!(lca_db.filename() == "");
    assert!(lca_db.moltype() == "DNA");
    assert!(lca_db._next_index() as u32 == 0);
    assert!(lca_db._next_lid() as u32 == 0);
    assert!(lca_db.ident_to_name() == HashMap::new());
    assert!(lca_db.ident_to_idx() == HashMap::new());
    assert!(lca_db.idx_to_lid() == HashMap::new());
    assert!(lca_db.lineage_to_lid() == HashMap::new());
    assert!(lca_db.lid_to_lineage() == HashMap::new());
    assert!(lca_db.hashval_to_idx() == HashMap::new());
}