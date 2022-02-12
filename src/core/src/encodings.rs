use std::collections::HashMap;
use std::convert::TryFrom;
use std::iter::Iterator;
use std::str;

use once_cell::sync::Lazy;

use crate::Error;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(u32)]
pub enum HashFunctions {
    murmur64_DNA = 1,
    murmur64_protein = 2,
    murmur64_dayhoff = 3,
    murmur64_hp = 4,
}

impl HashFunctions {
    pub fn dna(&self) -> bool {
        *self == HashFunctions::murmur64_DNA
    }

    pub fn protein(&self) -> bool {
        *self == HashFunctions::murmur64_protein
    }

    pub fn dayhoff(&self) -> bool {
        *self == HashFunctions::murmur64_dayhoff
    }

    pub fn hp(&self) -> bool {
        *self == HashFunctions::murmur64_hp
    }
}

impl std::fmt::Display for HashFunctions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                HashFunctions::murmur64_DNA => "dna",
                HashFunctions::murmur64_protein => "protein",
                HashFunctions::murmur64_dayhoff => "dayhoff",
                HashFunctions::murmur64_hp => "hp",
            }
        )
    }
}

impl TryFrom<&str> for HashFunctions {
    type Error = Error;

    fn try_from(moltype: &str) -> Result<Self, Self::Error> {
        match moltype.to_lowercase().as_ref() {
            "dna" => Ok(HashFunctions::murmur64_DNA),
            "dayhoff" => Ok(HashFunctions::murmur64_dayhoff),
            "hp" => Ok(HashFunctions::murmur64_hp),
            "protein" => Ok(HashFunctions::murmur64_protein),
            _ => unimplemented!(),
        }
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
            converted.push(aa_to_dayhoff(residue) as u8);
        } else if hp {
            converted.push(aa_to_hp(residue) as u8);
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
