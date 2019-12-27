use finch::sketch_schemes::mash::MashSketcher;
use finch::sketch_schemes::SketchScheme;

use crate::sketch::minhash::{HashFunctions, KmerMinHash};

/*
 TODO:
  - also convert scaled sketches
  - sourmash Signature equivalent is the finch Sketch, write conversions for that too
*/

impl From<MashSketcher> for KmerMinHash {
    fn from(other: MashSketcher) -> KmerMinHash {
        let values = other.to_vec();

        let mut new_mh = KmerMinHash::new(
            values.len() as u32,
            values.get(0).unwrap().kmer.len() as u32,
            HashFunctions::murmur64_DNA,
            42,
            0,
            true,
        );

        let hash_with_abunds: Vec<(u64, u64)> = values
            .iter()
            .map(|x| (x.hash as u64, x.count as u64))
            .collect();

        new_mh
            .add_many_with_abund(&hash_with_abunds)
            .expect("Error adding hashes with abund");

        new_mh
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;
    use std::collections::HashSet;
    use std::iter::FromIterator;

    use crate::signature::SigsTrait;
    use crate::sketch::minhash::{HashFunctions, KmerMinHash};

    use finch::sketch_schemes::mash::MashSketcher;
    use needletail::kmer::CanonicalKmers;
    use needletail::Sequence;

    use super::*;

    #[test]
    fn finch_behavior() {
        let mut a = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, true);
        let mut b = MashSketcher::new(20, 10, 42);

        let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";
        let rc = seq.reverse_complement();

        a.add_sequence(seq, false).unwrap();

        for (_, kmer, _) in CanonicalKmers::new(seq, &rc, 10) {
            b.push(&kmer, 0);
        }

        let b_hashes = b.to_vec();

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

    #[test]
    fn from_finch() {
        let mut a = KmerMinHash::new(20, 10, HashFunctions::murmur64_DNA, 42, 0, true);
        let mut b = MashSketcher::new(20, 10, 42);

        let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";
        let rc = seq.reverse_complement();

        a.add_sequence(seq, false).unwrap();

        for (_, kmer, _) in CanonicalKmers::new(seq, &rc, 10) {
            b.push(&kmer, 0);
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
}
