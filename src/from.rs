use finch::minhashes::MinHashKmers;

use KmerMinHash;

impl From<MinHashKmers> for KmerMinHash {
    fn from(other: MinHashKmers) -> KmerMinHash {
        let values = other.into_vec();

        let mut new_mh = KmerMinHash::new(
            values.len() as u32,
            values.get(0).unwrap().kmer.len() as u32,
            false,
            42,
            0,
            true,
        );

        let hash_with_abunds = values
            .iter()
            .map(|x| (x.hash as u64, x.count as u64))
            .collect();

        new_mh.add_many_with_abund(hash_with_abunds);

        new_mh
    }
}
