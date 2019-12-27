use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::{HashFunctions, KmerMinHash};

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
fn compare() {
    let mut a = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, false);
    let mut b = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, false);

    a.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    //    assert_eq!(b.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);
    //    assert_eq!(a.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    //    assert_eq!(b.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);
    //    assert_eq!(a.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT", false)
        .unwrap();
    assert!(a.compare(&b).unwrap() >= 0.3);
    assert!(b.compare(&a).unwrap() >= 0.3);
}

#[test]
fn similarity() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerMinHash::new(5, 20, HashFunctions::murmur64_hp, 42, 0, true);
    let mut b = KmerMinHash::new(5, 20, HashFunctions::murmur64_hp, 42, 0, true);

    a.add_hash(1);
    b.add_hash(1);
    b.add_hash(2);

    assert!((a.similarity(&a, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false)? - 0.5).abs() < 0.001);

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
        (a.similarity(&b, false)? - 0.705).abs() < 0.001,
        "{}",
        a.similarity(&b, false)?
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

    assert!((a.similarity(&a, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false)? - 0.23).abs() < 0.001);

    assert!((a.similarity(&a, true)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, true)? - 0.2).abs() < 0.001);

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
