#[cfg(feature = "from-finch")]
extern crate finch;

#[cfg(feature = "from-finch")]
extern crate needletail;

use std::collections::HashMap;
use std::collections::HashSet;
use std::iter::FromIterator;

use sourmash::KmerMinHash;

#[cfg(feature = "from-finch")]
use finch::minhashes::MinHashKmers;

#[cfg(feature = "from-finch")]
use needletail::kmer::canonical;

#[cfg(feature = "from-finch")]
#[test]
fn finch_behavior() {
    let mut a = KmerMinHash::new(20, 10, false, 42, 0, true);
    let mut b = MinHashKmers::new(20, 42);

    let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";

    a.add_sequence(seq, false);

    for kmer in seq.windows(10) {
        b.push(&canonical(kmer), 0);
    }

    let b_hashes = b.into_vec();

    let s1: HashSet<_> = HashSet::from_iter(a.mins.iter().map(|x| *x));
    let s2: HashSet<_> = HashSet::from_iter(b_hashes.iter().map(|x| x.hash as u64));
    let i1 = &s1 & &s2;

    assert!(i1.len() == a.mins.len());
    assert!(i1.len() == b_hashes.len());

    if let Some(abunds) = a.abunds {
        let smap: HashMap<_, _> = HashMap::from_iter(a.mins.iter().zip(abunds.iter()));
        println!("{:?}", smap);
        for item in b_hashes.iter() {
            assert!(smap.contains_key(&(item.hash as u64)));
            assert!(
                **smap.get(&(item.hash as u64)).unwrap()
                    == ((item.count + item.extra_count) as u64)
            );
        }
    }
}

#[cfg(feature = "from-finch")]
#[test]
fn from_finch() {
    let mut a = KmerMinHash::new(20, 10, false, 42, 0, true);
    let mut b = MinHashKmers::new(20, 42);

    let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";

    a.add_sequence(seq, false);

    for kmer in seq.windows(10) {
        b.push(&canonical(kmer), 0);
    }

    let c = KmerMinHash::from(b);

    let s1: HashSet<_> = HashSet::from_iter(a.mins.iter().map(|x| *x));
    let s2: HashSet<_> = HashSet::from_iter(c.mins.iter().map(|x| *x));
    let i1 = &s1 & &s2;

    assert!(i1.len() == a.mins.len());
    assert!(i1.len() == c.mins.len());

    if let Some(a_abunds) = a.abunds {
        if let Some(c_abunds) = c.abunds {
            let a_smap: HashMap<_, _> = HashMap::from_iter(a.mins.iter().zip(a_abunds.iter()));
            let c_smap: HashMap<_, _> = HashMap::from_iter(c.mins.iter().zip(c_abunds.iter()));
            for item in a_smap.iter() {
                assert!(c_smap.contains_key(*item.0));
                assert!(c_smap.get(*item.0).unwrap() == item.1);
            }
        }
    }
}
