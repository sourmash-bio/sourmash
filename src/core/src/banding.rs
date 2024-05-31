use std::ops::Range;

use crate::HashIntoType;
use crate::sketch::minhash::max_hash_for_scaled;

pub struct FracBand(usize, Range<usize>);

pub struct Bands(Vec<FracBand>);

#[derive(PartialEq, Debug)]
pub struct BandIdx {
    scaled: usize,
    n_bands: usize,
    idx: usize,
}

impl Bands {
    fn new(scaled: usize, n_bands: usize) -> Self {
        Self(
            (0..n_bands)
                .into_iter()
                .map(|idx| {
                    FracBand::from(BandIdx {
                        scaled,
                        n_bands,
                        idx,
                    })
                })
                .collect(),
        )
    }

    pub fn contains(&self, n: &HashIntoType) -> bool {
        self.0.iter().any(|band| band.contains(n))
    }

    fn to_vec(self) -> Vec<FracBand> {
        self.0
    }
}

impl FracBand {
    pub fn contains(&self, n: &HashIntoType) -> bool {
        self.1.contains(&(*n as usize))
    }
}

impl From<BandIdx> for FracBand {
    fn from(bidx: BandIdx) -> FracBand {
        let band_size = max_hash_for_scaled((bidx.scaled * bidx.n_bands) as u64) as usize;
        let range = Range {
            start: bidx.idx * band_size,
            end: (bidx.idx + 1) * band_size,
        };
        Self(bidx.scaled, range)
    }
}

impl From<&FracBand> for BandIdx {
    fn from(band: &FracBand) -> BandIdx {
        let scaled = band.0;
        let max_hash = max_hash_for_scaled(scaled as u64);
        let band_size = band.1.len();
        let n_bands = max_hash as usize / band_size;

        let idx = band.1.start / band_size;

        Self {
            scaled,
            n_bands,
            idx,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::HashIntoType;

    use super::{Bands, FracBand, BandIdx};

    #[test]
    fn band_to_idx() {
        let bands: Vec<FracBand> = Bands::new(100, 10).to_vec();
        for (n, band) in bands.into_iter().enumerate() {
            assert_eq!(
                BandIdx::from(&band),
                BandIdx {
                    scaled: 100,
                    n_bands: 10,
                    idx: n,
                }
            );
        }
    }

    // Bands(100, 10, 0) == Bands(1000, 1, 0)
    #[test]
    fn contains_1000() {
        let band_100 = FracBand::from(BandIdx { scaled: 100, n_bands: 10, idx: 0 });
        let band_1000 = FracBand::from(BandIdx { scaled: 1000, n_bands: 1, idx: 0 });

        for n in band_100.1.clone() {
            assert!(band_1000.contains(&(n as HashIntoType)));
        }

        for n in band_1000.1 {
            assert!(band_100.contains(&(n as HashIntoType)));
        }
    }

    /*
    // Bands(100, 10, 0) == Bands(1000, 1, 0)
    #[test]
    fn bands_contains_1000() {
        let band_1000 = FracBand::from(BandIdx { scaled: 1000, n_bands: 1, idx: 0 });
        let bands_100 = Bands::new(100, 10);

        for n in band_1000.1 {
            assert!(bands_100.contains(&(n as HashIntoType)));
        }
    }
    */

    // Bands(100, 100, 0) == Bands(10000, 1, 0)
    #[test]
    fn contains_10000() {
        let band_100 = FracBand::from(BandIdx { scaled: 100, n_bands: 100, idx: 0 });
        let band_10000 = FracBand::from(BandIdx { scaled: 10000, n_bands: 1, idx: 0 });

        for n in band_100.1.clone() {
            assert!(band_10000.contains(&(n as HashIntoType)));
        }

        for n in band_10000.1 {
            assert!(band_100.contains(&(n as HashIntoType)));
        }
    }
}
