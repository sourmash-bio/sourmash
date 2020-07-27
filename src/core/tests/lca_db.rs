use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::lca_db::LcaDB;
use sourmash::sketch::Sketch;

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
    assert!(lca_db.ident_to_name() == None);
    assert!(lca_db.ident_to_idx() == None);
    assert!(lca_db.idx_to_lid() == None);
    assert!(lca_db.lineage_to_lid() == None);
    assert!(lca_db.lid_to_lineage() == None);
    assert!(lca_db.hashval_to_idx() == None);
}