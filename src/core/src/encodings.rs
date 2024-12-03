use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault, Hash, Hasher};
use std::str;

use nohash_hasher::BuildNoHashHasher;
use once_cell::sync::Lazy;
use vec_collections::AbstractVecSet;

use crate::Error;

// To consider there: use a slab allocator for IdxTracker
// https://twitter.com/tomaka17/status/1391052081272967170
//   Pro-tip: you might be able to save a lot of hashmap lookups
//   if you replace a `HashMap<K, V>` with a `HashMap<K, usize>`
//   and a `Slab<V>`. This might be very useful if K is something
//   heavy such as a `String`.
pub type Color = u64;
pub type Idx = u32;
type IdxTracker = (vec_collections::VecSet<[Idx; 8]>, u64);
type ColorToIdx = HashMap<Color, IdxTracker, BuildNoHashHasher<Color>>;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Serialize, rkyv::Deserialize, rkyv::Archive)
)]
#[non_exhaustive]
pub enum HashFunctions {
    Murmur64Dna,
    Murmur64Protein,
    Murmur64Dayhoff,
    Murmur64Hp,
    Murmur64Skipm1n3,
    Murmur64Skipm2n3,
    Custom(String),
}

impl HashFunctions {
    pub fn dna(&self) -> bool {
        *self == HashFunctions::Murmur64Dna
    }

    pub fn protein(&self) -> bool {
        *self == HashFunctions::Murmur64Protein
    }

    pub fn dayhoff(&self) -> bool {
        *self == HashFunctions::Murmur64Dayhoff
    }

    pub fn hp(&self) -> bool {
        *self == HashFunctions::Murmur64Hp
    }
    pub fn skipm1n3(&self) -> bool {
        *self == HashFunctions::Murmur64Skipm1n3
    }
    pub fn skipm2n3(&self) -> bool {
        *self == HashFunctions::Murmur64Skipm2n3
    }
}

impl std::fmt::Display for HashFunctions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                HashFunctions::Murmur64Dna => "DNA",
                HashFunctions::Murmur64Protein => "protein",
                HashFunctions::Murmur64Dayhoff => "dayhoff",
                HashFunctions::Murmur64Hp => "hp",
                HashFunctions::Murmur64Skipm1n3 => "skipm1n3",
                HashFunctions::Murmur64Skipm2n3 => "skipm2n3",
                HashFunctions::Custom(v) => v,
            }
        )
    }
}

impl TryFrom<&str> for HashFunctions {
    type Error = Error;

    fn try_from(moltype: &str) -> Result<Self, Self::Error> {
        match moltype.to_lowercase().as_ref() {
            "dna" => Ok(HashFunctions::Murmur64Dna),
            "dayhoff" => Ok(HashFunctions::Murmur64Dayhoff),
            "hp" => Ok(HashFunctions::Murmur64Hp),
            "protein" => Ok(HashFunctions::Murmur64Protein),
            "skipm1n3" => Ok(HashFunctions::Murmur64Skipm1n3),
            "skipm2n3" => Ok(HashFunctions::Murmur64Skipm2n3),
            v => unimplemented!("{v}"),
        }
    }
}

#[derive(Debug)]
pub struct ReadingFrame {
    fw: Vec<u8>,         // Forward frame
    rc: Option<Vec<u8>>, // Reverse complement (optional, not used for protein input)
}

impl ReadingFrame {
    /// Create a k-mer iterator for this reading frame
    pub fn kmer_iter(&self, ksize: usize, seed: u64, force: bool) -> KmerIterator {
        KmerIterator::new(&self.fw, self.rc.as_deref(), ksize, seed, force)
    }
}

#[derive(Debug)]
pub struct ReadingFrames(Vec<ReadingFrame>);

impl ReadingFrames {
    /// Create ReadingFrames based on the input sequence, moltype, and protein flag
    pub fn new(sequence: &[u8], is_protein: bool, hash_function: &HashFunctions) -> Self {
        if is_protein {
            // for protein input, return one forward frame
            let frames = vec![ReadingFrame {
                fw: sequence.to_vec(),
                rc: None,
            }];
            Self(frames)
        } else if hash_function.dna() {
            // DNA: just forward + rc
            let dna_rc = revcomp(sequence);
            let frames = vec![ReadingFrame {
                fw: sequence.to_vec(),
                rc: Some(dna_rc),
            }];
            Self(frames)
        } else if hash_function.protein() || hash_function.dayhoff() || hash_function.hp() {
            // translation: build 6 frames
            let dna_rc = revcomp(sequence); // Compute reverse complement for translation
            let dayhoff = hash_function.dayhoff();
            let hp = hash_function.hp();
            Self::translate_frames(sequence, &dna_rc, dayhoff, hp)
        } else if hash_function.skipm1n3() || hash_function.skipm2n3() {
            // Skipmers: build 6 frames, following skip pattern
            let (m, n) = if hash_function.skipm1n3() {
                (1, 3)
            } else {
                (2, 3)
            };
            Self::skipmer_frames(sequence, n, m)
        } else {
            panic!("Unsupported moltype: {}", hash_function);
        }
    }

    /// Generate translated frames
    fn translate_frames(sequence: &[u8], dna_rc: &[u8], dayhoff: bool, hp: bool) -> Self {
        let frames: Vec<ReadingFrame> = (0..3)
            .map(|frame_number| ReadingFrame {
                fw: translated_frame(sequence, frame_number, dayhoff, hp),
                rc: Some(translated_frame(dna_rc, frame_number, dayhoff, hp)),
            })
            .collect();
        Self(frames)
    }

    /// Generate skipmer frames
    fn skipmer_frames(sequence: &[u8], n: usize, m: usize) -> Self {
        let frames: Vec<ReadingFrame> = (0..3)
            .map(|frame_number| {
                let fw = skipmer_frame(sequence, frame_number, n, m);
                ReadingFrame {
                    fw: fw.clone(),
                    rc: Some(revcomp(&fw)),
                }
            })
            .collect();
        Self(frames)
    }

    /// Access the frames
    pub fn frames(&self) -> &Vec<ReadingFrame> {
        &self.0
    }
}

pub struct KmerIterator<'a> {
    fw: &'a [u8],
    rc: Option<&'a [u8]>,
    ksize: usize,
    index: usize,
    len: usize,
    seed: u64,
    force: bool,
}

impl<'a> KmerIterator<'a> {
    pub fn new(fw: &'a [u8], rc: Option<&'a [u8]>, ksize: usize, seed: u64, force: bool) -> Self {
        Self {
            fw,
            rc,
            ksize,
            index: 0,
            len: fw.len(),
            seed,
            force,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.ksize > self.len {
            return None; // End of iteration
        }

        // Forward k-mer
        let kmer = &self.fw[self.index..self.index + self.ksize];

        // Validate the k-mer
        for j in self.index..self.index + self.ksize {
            if !VALID[self.fw[j] as usize] {
                if !self.force {
                    return Some(Err(Error::InvalidDNA {
                        message: String::from_utf8(kmer.to_vec()).unwrap(),
                    }));
                } else {
                    self.index += 1;
                    return Some(Ok(0)); // Skip invalid k-mer
                }
            }
        }

        // Reverse complement k-mer (if rc exists)

        // ... and then while moving the k-mer window forward for the sequence
        // we move another window backwards for the RC.
        //   For a ksize = 3, and a sequence AGTCGT (len = 6):
        //                   +-+---------+---------------+-------+
        //   seq      RC     |i|i + ksize|len - ksize - i|len - i|
        //  AGTCGT   ACGACT  +-+---------+---------------+-------+
        //  +->         +->  |0|    2    |       3       |   6   |
        //   +->       +->   |1|    3    |       2       |   5   |
        //    +->     +->    |2|    4    |       1       |   4   |
        //     +->   +->     |3|    5    |       0       |   3   |
        //                   +-+---------+---------------+-------+
        // (leaving this table here because I had to draw to
        //  get the indices correctly)
        let hash = if let Some(rc) = self.rc {
            let krc = &rc[self.len - self.ksize - self.index..self.len - self.index];
            crate::_hash_murmur(std::cmp::min(kmer, krc), self.seed)
        } else {
            crate::_hash_murmur(kmer, self.seed) // Use only forward k-mer if rc is None
        };

        self.index += 1;
        Some(Ok(hash))
    }
}

const COMPLEMENT: [u8; 256] = {
    let mut lookup = [0; 256];
    lookup[b'A' as usize] = b'T';
    lookup[b'C' as usize] = b'G';
    lookup[b'G' as usize] = b'C';
    lookup[b'T' as usize] = b'A';
    lookup[b'N' as usize] = b'N';
    lookup
};

#[inline]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|nt| COMPLEMENT[*nt as usize])
        .collect()
}

/// Generate a single translated frame from a DNA sequence
///
/// * `sequence`: The DNA sequence as a slice of bytes.
/// * `frame_number`: The frame to translate (0, 1, or 2).
/// * `dayhoff`: Whether to use the Dayhoff amino acid alphabet.
/// * `hp`: Whether to use the hydrophobic-polar amino acid alphabet.
///
/// Returns a translated frame as a `Vec<u8>`.
pub fn translated_frame(sequence: &[u8], frame_number: usize, dayhoff: bool, hp: bool) -> Vec<u8> {
    if frame_number > 2 {
        panic!("Frame number must be 0, 1, or 2");
    }

    sequence
        .iter()
        .cloned()
        .skip(frame_number) // Skip the initial bases for the frame
        .take(sequence.len() - frame_number) // Adjust length based on skipped bases
        .collect::<Vec<u8>>() // Collect the DNA subsequence
        .chunks(3) // Group into codons (triplets)
        .filter_map(|codon| to_aa(codon, dayhoff, hp).ok()) // Translate each codon
        .flatten() // Flatten the nested results into a single sequence
        .collect()
}

fn skipmer_frame(seq: &[u8], start: usize, n: usize, m: usize) -> Vec<u8> {
    seq.iter()
        .skip(start)
        .enumerate()
        .filter_map(|(i, &base)| if i % n < m { Some(base) } else { None })
        .collect()
}

static CODONTABLE: Lazy<HashMap<&'static str, u8>> = Lazy::new(|| {
    [
        // F
        ("TTT", b'F'),
        ("TTC", b'F'),
        // L
        ("TTA", b'L'),
        ("TTG", b'L'),
        // S
        ("TCT", b'S'),
        ("TCC", b'S'),
        ("TCA", b'S'),
        ("TCG", b'S'),
        ("TCN", b'S'),
        // Y
        ("TAT", b'Y'),
        ("TAC", b'Y'),
        // *
        ("TAA", b'*'),
        ("TAG", b'*'),
        // *
        ("TGA", b'*'),
        // C
        ("TGT", b'C'),
        ("TGC", b'C'),
        // W
        ("TGG", b'W'),
        // L
        ("CTT", b'L'),
        ("CTC", b'L'),
        ("CTA", b'L'),
        ("CTG", b'L'),
        ("CTN", b'L'),
        // P
        ("CCT", b'P'),
        ("CCC", b'P'),
        ("CCA", b'P'),
        ("CCG", b'P'),
        ("CCN", b'P'),
        // H
        ("CAT", b'H'),
        ("CAC", b'H'),
        // Q
        ("CAA", b'Q'),
        ("CAG", b'Q'),
        // R
        ("CGT", b'R'),
        ("CGC", b'R'),
        ("CGA", b'R'),
        ("CGG", b'R'),
        ("CGN", b'R'),
        // I
        ("ATT", b'I'),
        ("ATC", b'I'),
        ("ATA", b'I'),
        // M
        ("ATG", b'M'),
        // T
        ("ACT", b'T'),
        ("ACC", b'T'),
        ("ACA", b'T'),
        ("ACG", b'T'),
        ("ACN", b'T'),
        // N
        ("AAT", b'N'),
        ("AAC", b'N'),
        // K
        ("AAA", b'K'),
        ("AAG", b'K'),
        // S
        ("AGT", b'S'),
        ("AGC", b'S'),
        // R
        ("AGA", b'R'),
        ("AGG", b'R'),
        // V
        ("GTT", b'V'),
        ("GTC", b'V'),
        ("GTA", b'V'),
        ("GTG", b'V'),
        ("GTN", b'V'),
        // A
        ("GCT", b'A'),
        ("GCC", b'A'),
        ("GCA", b'A'),
        ("GCG", b'A'),
        ("GCN", b'A'),
        // D
        ("GAT", b'D'),
        ("GAC", b'D'),
        // E
        ("GAA", b'E'),
        ("GAG", b'E'),
        // G
        ("GGT", b'G'),
        ("GGC", b'G'),
        ("GGA", b'G'),
        ("GGG", b'G'),
        ("GGN", b'G'),
    ]
    .iter()
    .cloned()
    .collect()
});

// Dayhoff table from
// Peris, P., López, D., & Campos, M. (2008).
// IgTM: An algorithm to predict transmembrane domains and topology in
// proteins. BMC Bioinformatics, 9(1), 1029–11.
// http://doi.org/10.1186/1471-2105-9-367
//
// Original source:
// Dayhoff M. O., Schwartz R. M., Orcutt B. C. (1978).
// A model of evolutionary change in proteins,
// in Atlas of Protein Sequence and Structure,
// ed Dayhoff M. O., editor.
// (Washington, DC: National Biomedical Research Foundation; ), 345–352.
//
// | Amino acid    | Property              | Dayhoff |
// |---------------|-----------------------|---------|
// | C             | Sulfur polymerization | a       |
// | A, G, P, S, T | Small                 | b       |
// | D, E, N, Q    | Acid and amide        | c       |
// | H, K, R       | Basic                 | d       |
// | I, L, M, V    | Hydrophobic           | e       |
// | F, W, Y       | Aromatic              | f       |
static DAYHOFFTABLE: Lazy<HashMap<u8, u8>> = Lazy::new(|| {
    [
        // a
        (b'C', b'a'),
        // b
        (b'A', b'b'),
        (b'G', b'b'),
        (b'P', b'b'),
        (b'S', b'b'),
        (b'T', b'b'),
        // c
        (b'D', b'c'),
        (b'E', b'c'),
        (b'N', b'c'),
        (b'Q', b'c'),
        // d
        (b'H', b'd'),
        (b'K', b'd'),
        (b'R', b'd'),
        // e
        (b'I', b'e'),
        (b'L', b'e'),
        (b'M', b'e'),
        (b'V', b'e'),
        // e
        (b'F', b'f'),
        (b'W', b'f'),
        (b'Y', b'f'),
        // stop aa
        (b'*', b'*'),
    ]
    .iter()
    .cloned()
    .collect()
});

// HP Hydrophobic/hydrophilic mapping
// From: Phillips, R., Kondev, J., Theriot, J. (2008).
// Physical Biology of the Cell. New York: Garland Science, Taylor & Francis Group. ISBN: 978-0815341635

//
// | Amino acid                            | HP
// |---------------------------------------|---------|
// | A, F, G, I, L, M, P, V, W, Y          | h       |
// | N, C, S, T, D, E, R, H, K, Q          | p       |
static HPTABLE: Lazy<HashMap<u8, u8>> = Lazy::new(|| {
    [
        // h
        (b'A', b'h'),
        (b'F', b'h'),
        (b'G', b'h'),
        (b'I', b'h'),
        (b'L', b'h'),
        (b'M', b'h'),
        (b'P', b'h'),
        (b'V', b'h'),
        (b'W', b'h'),
        (b'Y', b'h'),
        // p
        (b'N', b'p'),
        (b'C', b'p'),
        (b'S', b'p'),
        (b'T', b'p'),
        (b'D', b'p'),
        (b'E', b'p'),
        (b'R', b'p'),
        (b'H', b'p'),
        (b'K', b'p'),
        (b'Q', b'p'),
        // stop aa
        (b'*', b'*'),
    ]
    .iter()
    .cloned()
    .collect()
});

#[inline]
pub fn translate_codon(codon: &[u8]) -> Result<u8, Error> {
    if codon.len() == 1 {
        return Ok(b'X');
    }

    if codon.len() == 2 {
        let mut v = codon.to_vec();
        v.push(b'N');
        match CODONTABLE.get(str::from_utf8(v.as_slice()).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    if codon.len() == 3 {
        match CODONTABLE.get(str::from_utf8(codon).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    Err(Error::InvalidCodonLength {
        message: format!("{}", codon.len()),
    })
}

#[inline]
pub fn aa_to_dayhoff(aa: u8) -> u8 {
    match DAYHOFFTABLE.get(&aa) {
        Some(letter) => *letter,
        None => b'X',
    }
}

pub fn aa_to_hp(aa: u8) -> u8 {
    match HPTABLE.get(&aa) {
        Some(letter) => *letter,
        None => b'X',
    }
}

#[inline]
pub fn to_aa(seq: &[u8], dayhoff: bool, hp: bool) -> Result<Vec<u8>, Error> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() < 3 {
            break;
        }

        let residue = translate_codon(chunk)?;
        if dayhoff {
            converted.push(aa_to_dayhoff(residue));
        } else if hp {
            converted.push(aa_to_hp(residue));
        } else {
            converted.push(residue);
        }
    }

    Ok(converted)
}

pub const VALID: [bool; 256] = {
    let mut lookup = [false; 256];
    lookup[b'A' as usize] = true;
    lookup[b'C' as usize] = true;
    lookup[b'G' as usize] = true;
    lookup[b'T' as usize] = true;
    lookup
};

#[derive(Serialize, Deserialize, Default)]
pub struct Colors {
    colors: ColorToIdx,
}

impl Colors {
    pub fn new() -> Colors {
        Default::default()
    }

    /// Given a color and a new idx, return an updated color
    ///
    /// This might create a new one, or find an already existing color
    /// that contains the new_idx
    ///
    /// Future optimization: store a count for each color, so we can track
    /// if there are extra colors that can be removed at the end.
    /// (the count is decreased whenever a new color has to be created)
    pub fn update<'a, I: IntoIterator<Item = &'a Idx>>(
        &mut self,
        current_color: Option<Color>,
        new_idxs: I,
    ) -> Result<Color, Error> {
        if let Some(color) = current_color {
            if let Some(idxs) = self.colors.get_mut(&color) {
                let idx_to_add: Vec<_> = new_idxs
                    .into_iter()
                    .filter(|new_idx| !idxs.0.contains(new_idx))
                    .collect();

                if idx_to_add.is_empty() {
                    // Easy case, it already has all the new_idxs, so just return this color
                    idxs.1 += 1;
                    Ok(color)
                } else {
                    // We need to either create a new color,
                    // or find an existing color that have the same idxs

                    let mut idxs = idxs.clone();
                    idxs.0.extend(idx_to_add.into_iter().cloned());
                    let new_color = Colors::compute_color(&idxs);

                    if new_color != color {
                        self.colors.get_mut(&color).unwrap().1 -= 1;
                        if self.colors[&color].1 == 0 {
                            self.colors.remove(&color);
                        };
                    };

                    self.colors
                        .entry(new_color)
                        .and_modify(|old_idxs| {
                            assert_eq!(old_idxs.0, idxs.0);
                            old_idxs.1 += 1;
                        })
                        .or_insert_with(|| (idxs.0, 1));
                    Ok(new_color)
                }
            } else {
                unimplemented!("throw error, current_color must exist in order to be updated. current_color: {:?}, colors: {:#?}", current_color, &self.colors);
            }
        } else {
            let mut idxs = IdxTracker::default();
            idxs.0.extend(new_idxs.into_iter().cloned());
            idxs.1 = 1;
            let new_color = Colors::compute_color(&idxs);
            self.colors
                .entry(new_color)
                .and_modify(|old_idxs| {
                    assert_eq!(old_idxs.0, idxs.0);
                    old_idxs.1 += 1;
                })
                .or_insert_with(|| (idxs.0, 1));
            Ok(new_color)
        }
    }

    fn compute_color(idxs: &IdxTracker) -> Color {
        let s = BuildHasherDefault::<twox_hash::Xxh3Hash128>::default();
        let mut hasher = s.build_hasher();
        idxs.0.hash(&mut hasher);
        hasher.finish()
    }

    pub fn len(&self) -> usize {
        self.colors.len()
    }

    pub fn is_empty(&self) -> bool {
        self.colors.is_empty()
    }

    pub fn contains(&self, color: Color, idx: Idx) -> bool {
        if let Some(idxs) = self.colors.get(&color) {
            idxs.0.contains(&idx)
        } else {
            false
        }
    }

    pub fn indices(&self, color: &Color) -> Indices {
        // TODO: what if color is not present?
        Indices {
            iter: self.colors.get(color).unwrap().0.iter(),
        }
    }

    pub fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&Color, &mut IdxTracker) -> bool,
    {
        self.colors.retain(f)
    }
}

pub struct Indices<'a> {
    iter: vec_collections::VecSetIter<core::slice::Iter<'a, Idx>>,
}

impl<'a> Iterator for Indices<'a> {
    type Item = &'a Idx;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn colors_update() {
        let mut colors = Colors::new();

        let color = colors.update(None, &[1_u32]).unwrap();
        assert_eq!(colors.len(), 1);

        dbg!("update");
        let new_color = colors.update(Some(color), &[1_u32]).unwrap();
        assert_eq!(colors.len(), 1);
        assert_eq!(color, new_color);

        dbg!("upgrade");
        let new_color = colors.update(Some(color), &[2_u32]).unwrap();
        assert_eq!(colors.len(), 2);
        assert_ne!(color, new_color);
    }

    #[test]
    fn colors_retain() {
        let mut colors = Colors::new();

        let color1 = colors.update(None, &[1_u32]).unwrap();
        assert_eq!(colors.len(), 1);
        // used_colors:
        //   color1: 1

        dbg!("update");
        let same_color = colors.update(Some(color1), &[1_u32]).unwrap();
        assert_eq!(colors.len(), 1);
        assert_eq!(color1, same_color);
        // used_colors:
        //   color1: 2

        dbg!("upgrade");
        let color2 = colors.update(Some(color1), &[2_u32]).unwrap();
        assert_eq!(colors.len(), 2);
        assert_ne!(color1, color2);
        // used_colors:
        //   color1: 1
        //   color2: 1

        dbg!("update");
        let same_color = colors.update(Some(color2), &[2_u32]).unwrap();
        assert_eq!(colors.len(), 2);
        assert_eq!(color2, same_color);
        // used_colors:
        //   color1: 1
        //   color1: 2

        dbg!("upgrade");
        let color3 = colors.update(Some(color1), &[3_u32]).unwrap();
        assert_ne!(color1, color3);
        assert_ne!(color2, color3);
        // used_colors:
        //   color1: 0
        //   color2: 2
        //   color3: 1

        // This is the pre color-count tracker, where it is needed
        // to call retain to maintain colors
        //assert_eq!(colors.len(), 3);
        //colors.retain(|c, _| [color2, color3].contains(c));

        assert_eq!(colors.len(), 2);
    }
}
