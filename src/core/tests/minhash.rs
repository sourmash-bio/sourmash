use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{
    max_hash_for_scaled, HashFunctions, KmerMinHash, KmerMinHashBTree,
};
use sourmash::sketch::Sketch;

use proptest::collection::vec;
use proptest::num::u64;
use proptest::proptest;

#[test]
fn throws_error() {
    let mut mh = KmerMinHash::new(1, 4, HashFunctions::murmur64_DNA, 42, 0, false);

    match mh.add_sequence(b"ATGR", false) {
        Ok(_) => assert!(false, "R is not a valid DNA character"),
        Err(_) => assert!(true),
    }
}

#[test]
fn merge() {
    let mut a = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, false);
    let mut b = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, false);

    a.add_sequence(b"TGCCGCCCAGCA", false).unwrap();
    b.add_sequence(b"TGCCGCCCAGCA", false).unwrap();

    a.add_sequence(b"GTCCGCCCAGTGA", false).unwrap();
    b.add_sequence(b"GTCCGCCCAGTGG", false).unwrap();

    a.merge(&b).unwrap();
    assert_eq!(
        a.to_vec(),
        vec![
            2996412506971915891,
            4448613756639084635,
            8373222269469409550,
            9390240264282449587,
            11085758717695534616,
            11668188995231815419,
            11760449009842383350,
            14682565545778736889,
        ]
    );
}

#[test]
fn invalid_dna() {
    let mut a = KmerMinHash::new(20, 3, HashFunctions::murmur64_DNA, 42, 0, false);

    a.add_sequence(b"AAANNCCCTN", true).unwrap();
    assert_eq!(a.mins().len(), 3);

    let mut b = KmerMinHash::new(20, 3, HashFunctions::murmur64_DNA, 42, 0, false);
    b.add_sequence(b"NAAA", true).unwrap();
    assert_eq!(b.mins().len(), 1);
}

#[test]
fn similarity() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerMinHash::new(5, 20, HashFunctions::murmur64_hp, 42, 0, true);
    let mut b = KmerMinHash::new(5, 20, HashFunctions::murmur64_hp, 42, 0, true);

    a.add_hash(1);
    b.add_hash(1);
    b.add_hash(2);

    assert!((a.similarity(&a, false, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false, false)? - 0.5).abs() < 0.001);

    Ok(())
}

#[test]
fn similarity_2() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerMinHash::new(5, 5, HashFunctions::murmur64_DNA, 42, 0, true);
    let mut b = KmerMinHash::new(5, 5, HashFunctions::murmur64_DNA, 42, 0, true);

    a.add_sequence(b"ATGGA", false)?;
    a.add_sequence(b"GGACA", false)?;

    a.add_sequence(b"ATGGA", false)?;
    b.add_sequence(b"ATGGA", false)?;

    assert!(
        (a.similarity(&b, false, false)? - 0.705).abs() < 0.001,
        "{}",
        a.similarity(&b, false, false)?
    );

    Ok(())
}

#[test]
fn similarity_3() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerMinHash::new(5, 20, HashFunctions::murmur64_dayhoff, 42, 0, true);
    let mut b = KmerMinHash::new(5, 20, HashFunctions::murmur64_dayhoff, 42, 0, true);

    a.add_hash(1);
    a.add_hash(1);
    a.add_hash(5);
    a.add_hash(5);

    b.add_hash(1);
    b.add_hash(2);
    b.add_hash(3);
    b.add_hash(4);

    assert!((a.similarity(&a, false, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false, false)? - 0.23).abs() < 0.001);

    assert!((a.similarity(&a, true, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, true, false)? - 0.2).abs() < 0.001);

    Ok(())
}

#[test]
fn dayhoff() {
    let mut a = KmerMinHash::new(10, 6, HashFunctions::murmur64_dayhoff, 42, 0, false);
    let mut b = KmerMinHash::new(10, 6, HashFunctions::murmur64_protein, 42, 0, false);

    a.add_sequence(b"ACTGAC", false).unwrap();
    b.add_sequence(b"ACTGAC", false).unwrap();

    assert_eq!(a.size(), 2);
    assert_eq!(b.size(), 2);
}

#[test]
fn hp() {
    let mut a = KmerMinHash::new(10, 6, HashFunctions::murmur64_hp, 42, 0, false);
    let mut b = KmerMinHash::new(10, 6, HashFunctions::murmur64_protein, 42, 0, false);

    a.add_sequence(b"ACTGAC", false).unwrap();
    b.add_sequence(b"ACTGAC", false).unwrap();

    assert_eq!(a.size(), 2);
    assert_eq!(b.size(), 2);
}

#[test]
fn max_for_scaled() {
    assert_eq!(max_hash_for_scaled(100), Some(184467440737095520));
}

proptest! {
#[test]
fn oracle_mins(hashes in vec(u64::ANY, 1..10000)) {
    let mut a = KmerMinHash::new(10, 6, HashFunctions::murmur64_DNA, 42, 0, true);
    let mut b = KmerMinHashBTree::new(10, 6, HashFunctions::murmur64_DNA, 42, 0, true);

    let mut c = KmerMinHash::new(10, 6, HashFunctions::murmur64_DNA, 42, 0, true);
    let mut d = KmerMinHashBTree::new(10, 6, HashFunctions::murmur64_DNA, 42, 0, true);

    for hash in hashes {
        a.add_hash(hash);
        b.add_hash(hash);

        if hash % 2 == 0 {
          c.add_hash(hash);
          d.add_hash(hash);
        }
    }
    assert_eq!(a.mins(), b.mins());
    assert_eq!(c.mins(), d.mins());

    assert_eq!(a.abunds(), b.abunds());
    assert_eq!(c.abunds(), d.abunds());

    assert_eq!(a.similarity(&c, false, false).unwrap(), b.similarity(&d, false, false).unwrap());
}
}

proptest! {
#[test]
fn oracle_mins_scaled(hashes in vec(u64::ANY, 1..10000)) {
    let max_hash = max_hash_for_scaled(100).unwrap();
    let mut a = KmerMinHash::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);
    let mut b = KmerMinHashBTree::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);

    let mut c = KmerMinHash::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);
    let mut d = KmerMinHashBTree::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);

    for hash in hashes {
        a.add_hash(hash);
        b.add_hash(hash);

        if hash % 2 == 0 {
          c.add_hash(hash);
          d.add_hash(hash);
        }
    }

    assert_eq!(a.mins(), b.mins());
    assert_eq!(c.mins(), d.mins());

    assert_eq!(a.md5sum(), b.md5sum());
    assert_eq!(c.md5sum(), d.md5sum());

    assert_eq!(a.is_protein(), b.is_protein());
    assert_eq!(a.num(), b.num());
    assert_eq!(a.seed(), b.seed());
    assert_eq!(a.ksize(), b.ksize());
    assert_eq!(a.max_hash(), b.max_hash());
    assert_eq!(a.track_abundance(), b.track_abundance());
    assert_eq!(a.hash_function(), b.hash_function());

    assert_eq!(a.abunds(), b.abunds());
    assert_eq!(c.abunds(), d.abunds());

    assert_eq!(a.similarity(&c, false, false).unwrap(), b.similarity(&d, false, false).unwrap());
    assert_eq!(a.similarity(&c, true, false).unwrap(), b.similarity(&d, true, false).unwrap());
}
}

proptest! {
#[test]
fn prop_merge(seq1 in "[ACGT]{6,100}", seq2 in "[ACGT]{6,200}") {
    let max_hash = max_hash_for_scaled(10).unwrap();
    let mut a = KmerMinHash::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);
    let mut b = KmerMinHashBTree::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);

    let mut c = KmerMinHash::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);
    let mut d = KmerMinHashBTree::new(0, 6, HashFunctions::murmur64_DNA, 42, max_hash, true);

    a.add_sequence(seq1.as_bytes(), false).unwrap();
    b.add_sequence(seq1.as_bytes(), false).unwrap();

    c.add_sequence(seq2.as_bytes(), false).unwrap();
    d.add_sequence(seq2.as_bytes(), false).unwrap();

    a.merge(&c).unwrap();
    b.merge(&d).unwrap();

    assert_eq!(a.mins(), b.mins());
    assert_eq!(c.mins(), d.mins());

    assert_eq!(a.abunds(), b.abunds());
    assert_eq!(c.abunds(), d.abunds());

    assert_eq!(a.similarity(&c, false, false).unwrap(), b.similarity(&d, false, false).unwrap());
}
}

#[test]
fn load_save_minhash_sketches() {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/genome-s10+s11.sig");

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

    let sig = sigs.get(0).unwrap();
    let sketches = sig.sketches();
    let mut buffer = vec![];

    if let Sketch::MinHash(mh) = &sketches[0] {
        let bmh: KmerMinHashBTree = mh.clone().into();
        {
            serde_json::to_writer(&mut buffer, &bmh).unwrap();
        }

        let new_mh: KmerMinHash = serde_json::from_reader(&buffer[..]).unwrap();
        let new_bmh: KmerMinHashBTree = serde_json::from_reader(&buffer[..]).unwrap();

        assert_eq!(mh.md5sum(), new_mh.md5sum());
        assert_eq!(bmh.md5sum(), new_bmh.md5sum());
        assert_eq!(bmh.md5sum(), new_mh.md5sum());
        assert_eq!(mh.md5sum(), new_bmh.md5sum());

        assert_eq!(mh.mins(), new_mh.mins());
        assert_eq!(bmh.mins(), new_bmh.mins());
        assert_eq!(bmh.mins(), new_mh.mins());
        assert_eq!(mh.mins(), new_bmh.mins());

        assert_eq!(mh.abunds(), new_mh.abunds());
        assert_eq!(bmh.abunds(), new_bmh.abunds());
        assert_eq!(bmh.abunds(), new_mh.abunds());
        assert_eq!(mh.abunds(), new_bmh.abunds());

        assert_eq!(
            mh.similarity(&new_mh, false, false).unwrap(),
            bmh.similarity(&new_bmh, false, false).unwrap()
        );

        assert_eq!(
            mh.similarity(&new_mh, true, false).unwrap(),
            bmh.similarity(&new_bmh, true, false).unwrap()
        );

        buffer.clear();
        let imh: KmerMinHash = bmh.clone().into();
        {
            serde_json::to_writer(&mut buffer, &imh).unwrap();
        }

        let new_mh: KmerMinHash = serde_json::from_reader(&buffer[..]).unwrap();
        let new_bmh: KmerMinHashBTree = serde_json::from_reader(&buffer[..]).unwrap();

        assert_eq!(mh.md5sum(), new_mh.md5sum());
        assert_eq!(bmh.md5sum(), new_bmh.md5sum());
        assert_eq!(bmh.md5sum(), new_mh.md5sum());
        assert_eq!(mh.md5sum(), new_bmh.md5sum());

        assert_eq!(mh.mins(), new_mh.mins());
        assert_eq!(bmh.mins(), new_bmh.mins());
        assert_eq!(bmh.mins(), new_mh.mins());
        assert_eq!(mh.mins(), new_bmh.mins());

        assert_eq!(mh.abunds(), new_mh.abunds());
        assert_eq!(bmh.abunds(), new_bmh.abunds());
        assert_eq!(bmh.abunds(), new_mh.abunds());
        assert_eq!(mh.abunds(), new_bmh.abunds());

        assert_eq!(
            mh.similarity(&new_mh, false, false).unwrap(),
            bmh.similarity(&new_bmh, false, false).unwrap()
        );

        assert_eq!(
            mh.similarity(&new_mh, true, false).unwrap(),
            bmh.similarity(&new_bmh, true, false).unwrap()
        );
    }
}

// TODO: test btree ::default()
// TODO: check .clear(), .is_empty()
// TODO: set_hash_function?
// TODO: enable/disable abundance
// TODO: downsample_max_hash
// TODO: remove, remove_many
// TODO: add_from, add_many, add_many_with_abund
// TODO: count_common, intersection, intersection_size
