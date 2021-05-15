use finch::sketch_schemes::mash::MashSketcher;
use finch::sketch_schemes::SketchScheme;

use crate::encodings::HashFunctions;
use crate::sketch::minhash::KmerMinHash;

/*
 TODO:
  - also convert scaled sketches
  - sourmash Signature equivalent is the finch Sketch, write conversions for that too
*/

impl From<MashSketcher> for KmerMinHash {
    fn from(other: MashSketcher) -> KmerMinHash {
        let values = other.to_vec();

        let mut new_mh = KmerMinHash::new(
            0,
            values.get(0).unwrap().kmer.len() as u32,
            HashFunctions::murmur64_DNA,
            42,
            true,
            values.len() as u32,
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

    use crate::encodings::HashFunctions;
    use crate::signature::SigsTrait;
    use crate::sketch::minhash::KmerMinHash;

    use finch::sketch_schemes::mash::MashSketcher;
    use needletail::kmer::CanonicalKmers;
    use needletail::Sequence;

    use super::*;

    #[test]
    fn finch_behavior() {
        let mut a = KmerMinHash::new(0, 10, HashFunctions::murmur64_DNA, 42, true, 20);
        let mut b = MashSketcher::new(20, 10, 42);

        let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";
        let rc = seq.reverse_complement();

        a.add_sequence(seq, false).unwrap();

        for (_, kmer, _) in CanonicalKmers::new(seq, &rc, 10) {
            b.push(&kmer, 0);
        }

        let b_hashes = b.to_vec();

        let s1: HashSet<_> = a.mins().into_iter().collect();
        let s2: HashSet<_> = b_hashes.iter().map(|x| x.hash as u64).collect();
        let i1 = &s1 & &s2;

        assert!(i1.len() == a.size());
        assert!(i1.len() == b_hashes.len());

        if let Some(abunds) = a.abunds() {
            let mins = a.mins();
            let smap: HashMap<_, _> = mins.iter().zip(abunds.iter()).collect();
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
        let mut a = KmerMinHash::new(0, 10, HashFunctions::murmur64_DNA, 42, true, 20);
        let mut b = MashSketcher::new(20, 10, 42);

        let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";
        let rc = seq.reverse_complement();

        a.add_sequence(seq, false).unwrap();

        for (_, kmer, _) in CanonicalKmers::new(seq, &rc, 10) {
            b.push(&kmer, 0);
        }

        let c = KmerMinHash::from(b);

        let s1: HashSet<_> = a.mins().into_iter().collect();
        let s2: HashSet<_> = c.mins().into_iter().collect();
        let i1 = &s1 & &s2;

        assert!(i1.len() == a.mins().len());
        assert!(i1.len() == c.mins().len());

        if let Some(a_abunds) = a.abunds() {
            if let Some(c_abunds) = c.abunds() {
                let a_mins = a.mins();
                let a_smap: HashMap<_, _> = a_mins.iter().zip(a_abunds.iter()).collect();
                let c_mins = c.mins();
                let c_smap: HashMap<_, _> = c_mins.iter().zip(c_abunds.iter()).collect();
                for item in a_smap.iter() {
                    assert!(c_smap.contains_key(*item.0));
                    assert!(c_smap.get(*item.0).unwrap() == item.1);
                }
            }
        }
    }
}
