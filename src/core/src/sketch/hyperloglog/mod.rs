/*
Based on the HyperLogLog implementations in khmer
  https://github.com/dib-lab/khmer/blob/fb65d21eaedf0d397d49ae3debc578897f9d6eb4/src/oxli/hllcounter.cc
and pdatastructs
  https://github.com/crepererum/pdatastructs.rs/blob/0254dbce33c0404f444f29e0cca8f73002d6e5e6/src/hyperloglog.rs
*/

use std::cmp;

use serde::{Deserialize, Serialize};

use crate::encodings::HashFunctions;
use crate::signature::SigsTrait;
use crate::Error;
use crate::HashIntoType;

mod hyperloglog_data;
use hyperloglog_data::{
    BIAS_DATA_OFFSET, BIAS_DATA_VEC, RAW_ESTIMATE_DATA_OFFSET, RAW_ESTIMATE_DATA_VEC,
    THRESHOLD_DATA_OFFSET, THRESHOLD_DATA_VEC,
};

type CounterType = u8;

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct HyperLogLog {
    registers: Vec<CounterType>,
    p: usize,
    ksize: usize,
}

impl HyperLogLog {
    pub fn with_error_rate(error_rate: f64, ksize: usize) -> Result<HyperLogLog, Error> {
        let p = f64::ceil(f64::log2(f64::powi(1.04 / error_rate, 2)));
        HyperLogLog::new(p as usize, ksize)
    }

    pub fn new(p: usize, ksize: usize) -> Result<HyperLogLog, Error> {
        if p < 4 || p > 18 {
            return Err(Error::HLLPrecisionBounds);
        }

        let size = (1 as usize) << p;
        let registers = vec![0; size];

        Ok(HyperLogLog {
            registers,
            ksize,
            p,
        })
    }

    pub fn merge(&mut self, other: &HyperLogLog) -> Result<(), Error> {
        self.check_compatible(other)?;
        self.registers
            .iter_mut()
            .zip(other.registers.iter())
            .for_each(|(a, b)| *a = cmp::max(*a, *b));
        Ok(())
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = crate::_hash_murmur(word, 42); // TODO: decide on seed
        self.add_hash(hash);
    }

    pub fn add_many(&mut self, hashes: &[HashIntoType]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn cardinality(&self) -> usize {
        HyperLogLog::estimate_cardinality(&self.registers)
    }

    fn alpha(m: usize) -> f64 {
        match m {
            0..=15 => todo!("raise error_rate error, too small"),
            16..=32 => 0.673,
            33..=64 => 0.697,
            65..=128 => 0.709,
            262144..=usize::MAX => todo!("raise error_rate error, too big"),
            _ => 0.7213 / (1.0 + 1.079 / (m as f64)),
        }
    }

    pub fn estimate_cardinality(registers: &[CounterType]) -> usize {
        let m = registers.len() as f64;
        let v = bytecount::count(registers, 0) as f64;

        if v > 0. {
            // do linear counting
            let h = (m * (m / v).ln()) as usize;
            if h <= HyperLogLog::threshold(m as usize) {
                return h as usize;
            }
        }

        let z = 1.
            / registers
                .iter()
                .fold(0., |acc, &v| acc + 2f64.powi(-(i32::from(v))));

        let e = HyperLogLog::alpha(m as usize) * m * m * z;

        if e <= 5. * m {
            (e - HyperLogLog::estimate_bias(e, m as usize)) as usize
        } else {
            e as usize
        }
    }

    fn threshold(nregisters: usize) -> usize {
        let p = f64::log2(nregisters as f64) as usize;
        THRESHOLD_DATA_VEC[p - THRESHOLD_DATA_OFFSET]
    }

    fn neighbor_search_startpoints(lookup_array: &[f64], e: f64) -> (Option<usize>, Option<usize>) {
        // binary search first nearest neighbor
        match lookup_array.binary_search_by(|v| v.partial_cmp(&e).unwrap()) {
            Ok(i) => (Some(i), Some(i)),
            Err(i) => {
                if i == 0 {
                    // no left index
                    (None, Some(0))
                } else if i == lookup_array.len() {
                    // no right index
                    (Some(i - 1), None)
                } else {
                    (Some(i), Some(i + 1))
                }
            }
        }
    }

    fn estimate_bias(e: f64, nregisters: usize) -> f64 {
        let p = f64::log2(nregisters as f64) as usize;
        let lookup_array = RAW_ESTIMATE_DATA_VEC[p - RAW_ESTIMATE_DATA_OFFSET];
        let (mut idx_left, mut idx_right) = Self::neighbor_search_startpoints(lookup_array, e);

        // collect k nearest neighbors
        let k = 6;
        assert!(lookup_array.len() >= k);
        let mut neighbors = vec![];
        for _ in 0..k {
            let (right_instead_left, idx) = match (idx_left, idx_right) {
                (Some(i_left), Some(i_right)) => {
                    // 2 candidates, find better one
                    let delta_left = (lookup_array[i_left] - e).abs();
                    let delta_right = (lookup_array[i_right] - e).abs();
                    if delta_right < delta_left {
                        (true, i_right)
                    } else {
                        (false, i_left)
                    }
                }
                (Some(i_left), None) => {
                    // just left one is there, use it
                    (false, i_left)
                }
                (None, Some(i_right)) => {
                    // just right one is there, use it
                    (true, i_right)
                }
                _ => panic!("neighborhood search failed, this is bug!"),
            };
            neighbors.push(idx);
            if right_instead_left {
                idx_right = if idx < lookup_array.len() - 1 {
                    Some(idx + 1)
                } else {
                    None
                };
            } else {
                idx_left = if idx > 0 { Some(idx - 1) } else { None };
            }
        }

        // calculate mean of neighbors
        let bias_data = BIAS_DATA_VEC[p - BIAS_DATA_OFFSET];
        neighbors.iter().map(|&i| bias_data[i]).sum::<f64>() / (k as f64)
    }
}

impl SigsTrait for HyperLogLog {
    fn size(&self) -> usize {
        self.registers.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.registers.iter().map(|x| *x as u64).collect()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn seed(&self) -> u64 {
        // TODO: support other seeds
        42
    }

    fn hash_function(&self) -> HashFunctions {
        //TODO support other hash functions
        HashFunctions::murmur64_DNA
    }

    fn add_hash(&mut self, hash: HashIntoType) {
        let value = hash >> self.p;
        let index = (hash - (value << self.p)) as usize;

        let leftmost = value.leading_zeros() + 1 - (self.p as u32);

        let old_value = self.registers[index];
        self.registers[index] = cmp::max(old_value, leftmost as CounterType);
    }

    fn check_compatible(&self, other: &HyperLogLog) -> Result<(), Error> {
        if self.ksize() != other.ksize() {
            Err(Error::MismatchKSizes)
        } else if self.size() != other.size() {
            // TODO: create new error
            Err(Error::MismatchNum {
                n1: self.size() as u32,
                n2: other.size() as u32,
            })
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;
    use std::path::PathBuf;

    use crate::signature::SigsTrait;
    use needletail::{parse_fastx_file, Sequence};

    use super::HyperLogLog;

    // TODO: pull more tests from khmer HLL

    #[test]
    fn hll_add() {
        let ERR_RATE = 0.01;
        let N_UNIQUE = 3356;
        let KSIZE: u8 = 21;

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/ecoli.genes.fna");

        let mut hll = HyperLogLog::with_error_rate(ERR_RATE, KSIZE as usize).unwrap();
        let mut counter: HashSet<Vec<u8>> = HashSet::new();

        let mut parser = parse_fastx_file(filename).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            let norm_seq = record.normalize(false);
            let rc = norm_seq.reverse_complement();

            hll.add_sequence(&norm_seq, false).unwrap();
            for (_, kmer, _) in norm_seq.canonical_kmers(KSIZE, &rc) {
                counter.insert(kmer.into());
            }
        }

        assert_eq!(counter.len(), N_UNIQUE);
        let abs_error = (1. - (hll.cardinality() as f64 / N_UNIQUE as f64)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);
    }
}
