use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;

#[test]
fn throws_error() {
    let mut mh = KmerMinHash::new(1, 4, false, 42, 0, false);

    match mh.add_sequence(b"ATGR", false) {
        Ok(_) => assert!(false, "R is not a valid DNA character"),
        Err(_) => assert!(true),
    }
}

#[test]
fn merge() {
    let mut a = KmerMinHash::new(20, 10, false, 42, 0, false);
    let mut b = KmerMinHash::new(20, 10, false, 42, 0, false);

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
    let mut a = KmerMinHash::new(20, 10, false, 42, 0, false);
    let mut b = KmerMinHash::new(20, 10, false, 42, 0, false);

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
