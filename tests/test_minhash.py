# This file is part of sourmash, https://github.com/dib-lab/sourmash/, and is
# Copyright (C) 2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: titus@idyll.org
# pylint: disable=missing-docstring,protected-access

import itertools
import pickle
import math
import numpy as np

import pytest

import screed

import sourmash
from sourmash.minhash import (
    MinHash,
    FrozenMinHash,
    hash_murmur,
    _get_scaled_for_max_hash,
    _get_max_hash_for_scaled,
    translate_codon
)
from sourmash import signature

import sourmash_tst_utils as utils

# add:
# * get default params from Python
# * keyword args for minhash constructor
# * trap error from handing protein/non-DNA to a DNA MH
# * fail on untagged/unloaded countgraph
# * nan on empty minhash
# * define equals

scaled50 = _get_scaled_for_max_hash(50)
scaled100 = _get_scaled_for_max_hash(100)
scaled5000 = _get_scaled_for_max_hash(5000)


def _kmers_from_all_coding_frames(sequence, ksize):
    """Mimic the internal rust code for translation of DNA into aa.

    For each frame, yield all fwd k-mers, then all reverse k-mers
    from that frame. Then do next frame.
    """
    seqrc = screed.rc(sequence)

    for frame in (0, 1, 2):
        # get forward k-mers
        for start in range(0, len(sequence) - ksize + 1 - frame, 3):
            kmer = sequence[start + frame:start + frame + ksize]
            yield kmer

        # get rc k-mers
        for start in range(0, len(seqrc) - ksize + 1 - frame, 3):
            kmer = seqrc[start + frame:start + frame + ksize]
            yield kmer


def _hash_fwd_only(mh_translate, seq):
    "Return the first hashval only, for coding frame +1."
    assert len(seq) == mh_translate.ksize*3
    xx = mh_translate.seq_to_hashes(seq)[0]
    return xx


def test_basic_dna(track_abundance):
    # verify that MHs of size 1 stay size 1, & act properly as bottom sketches.
    mh = MinHash(1, 4, track_abundance=track_abundance)
    assert mh.moltype == 'DNA'

    mh.add_sequence('ATGC')
    a = mh.hashes

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.hashes

    print(a, b)
    assert list(a) == list(b)
    assert len(b) == 1
    assert list(a)[0] == list(b)[0] == 12415348535738636339


def test_div_zero(track_abundance):
    # verify that empty MHs do not yield divide by zero errors for similarity
    mh = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = mh.copy_and_clear()

    mh.add_sequence('ATGC')
    assert mh.similarity(mh2) == 0
    assert mh2.similarity(mh) == 0


def test_div_zero_contained(track_abundance):
    # verify that empty MHs do not yield divide by zero errors for contained_by
    mh = MinHash(0, 4, scaled=1, track_abundance=track_abundance)
    mh2 = mh.copy_and_clear()

    mh.add_sequence('ATGC')
    assert mh.contained_by(mh2) == 0
    assert mh2.contained_by(mh) == 0


def test_contained_requires_scaled(track_abundance):
    # test that contained_by requires scaled signatures
    mh1 = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)

    mh1.add_sequence('ATGC')
    mh2.add_sequence('ATGC')

    with pytest.raises(TypeError):
        mh2.contained_by(mh1)

    with pytest.raises(TypeError):
        mh1.contained_by(mh2)


def test_contained_requires_scaled_2(track_abundance):
    # test that max_containment requires scaled signatures
    mh1 = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)

    mh1.add_sequence('ATGC')
    mh2.add_sequence('ATGC')

    with pytest.raises(TypeError):
        mh2.max_containment(mh1)

    with pytest.raises(TypeError):
        mh1.max_containment(mh2)


def test_contained_requires_scaled_3(track_abundance):
    # test that avg_containment requires scaled signatures
    mh1 = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)

    mh1.add_sequence('ATGC')
    mh2.add_sequence('ATGC')

    with pytest.raises(TypeError):
        mh2.avg_containment(mh1)

    with pytest.raises(TypeError):
        mh1.avg_containment(mh2)


def test_bytes_dna(track_abundance):
    mh = MinHash(1, 4, track_abundance=track_abundance)
    mh.add_sequence('ATGC')
    mh.add_sequence(b'ATGC')
    mh.add_sequence('ATGC')
    a = mh.hashes

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    mh.add_sequence(b'GCAT')             # this will not get added; hash > ATGC
    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.hashes

    print(a, b)
    assert list(a) == list(b)
    assert len(b) == 1

def test_add_long_seqs_force():
    # Test for (All kmers are invalid)

    mh = sourmash.minhash.MinHash(n = 0, ksize=21, scaled =10, seed = 42)
    seq = "ACGTN" * 100000
    hashes = mh.seq_to_hashes(seq, force = True)
    assert(len(mh.hashes) == 0)


def test_seq_to_hashes(track_abundance):
    mh = sourmash.minhash.MinHash(n=0, ksize=21, scaled=1, track_abundance=track_abundance)
    seq = "ATGAGAGACGATAGACAGATGACC"
    mh.add_sequence(seq)

    golden_hashes = mh.hashes
    
    # New seq to hashes without adding to the sketch
    new_hashes = mh.seq_to_hashes(seq)

    assert set(golden_hashes) == set(new_hashes)


def test_seq_to_hashes_protein_1(track_abundance, dayhoff):
    mh = MinHash(10, 2, is_protein=True, dayhoff=dayhoff, hp=False, track_abundance=track_abundance)
    prot_seq = "AGYYG"

    mh.add_protein(prot_seq)

    golden_hashes = mh.hashes

    # New seq to hashes without adding to the sketch
    new_hashes = mh.seq_to_hashes(prot_seq, is_protein = True)

    assert set(golden_hashes) == set(new_hashes)

def test_seq_to_hashes_protein_2(track_abundance):
    mh = sourmash.minhash.MinHash(n=0, ksize=21, scaled=1, track_abundance=track_abundance)
    seq = "ATGAGAGACGATAGACAGATGACC"

    with pytest.raises(ValueError):
        mh.seq_to_hashes(seq, is_protein = True)


def test_seq_to_hashes_translated(track_abundance):
    mh_protein = MinHash(10, 2, is_protein=True, track_abundance=track_abundance)
    seq = "ACTGAC"
    mh_protein.add_sequence(seq)

    golden_hashes = mh_protein.hashes

    # New seq to hashes without adding to the sketch
    new_hashes = mh_protein.seq_to_hashes(seq)

    assert set(golden_hashes) == set(new_hashes)


def test_seq_to_hashes_bad_kmers_as_zeroes_1():
    mh = sourmash.minhash.MinHash(n=0, ksize=21, scaled=1)
    seq = "ATGAGAGACGATAGACAGATGACN"
    
    # New seq to hashes without adding to the sketch
    hashes = mh.seq_to_hashes(seq, force=True, bad_kmers_as_zeroes=True)

    assert len(hashes) == len(seq) - 21 + 1


def test_seq_to_hashes_bad_kmers_as_zeroes_2():
    mh = sourmash.minhash.MinHash(n=0, ksize=21, scaled=1)
    seq = "ATGAGAGACGATAGACAGATGACN"
    
    with pytest.raises(ValueError):
        hashes = mh.seq_to_hashes(seq, bad_kmers_as_zeroes=True)


def test_seq_to_hashes_translated_short():
    mh = MinHash(0, 2, is_protein=True, dayhoff=True, hp=False, scaled = 1)
    hashes = mh.seq_to_hashes("ACTGA")

    assert(len(hashes) == 0)


def test_bytes_protein_dayhoff(track_abundance, dayhoff):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 2, is_protein=True, dayhoff=dayhoff, hp=False,
                 track_abundance=track_abundance)

    expected_moltype = 'protein'
    if dayhoff:
        expected_moltype = 'dayhoff'
    assert mh.moltype == expected_moltype

    mh.add_protein('AGYYG')
    mh.add_protein('AGYYG')
    mh.add_protein(b'AGYYG')

    assert len(mh.hashes) == 4


def test_protein_dayhoff(track_abundance, dayhoff):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 2, is_protein=True, dayhoff=dayhoff, hp=False, track_abundance=track_abundance)
    mh.add_protein('AGYYG')

    assert len(mh.hashes) == 4


def test_bytes_protein_hp(track_abundance, hp):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 2, is_protein=True, dayhoff=False, hp=hp, track_abundance=track_abundance)
    expected_moltype = 'protein'
    if hp:
        expected_moltype = 'hp'
    assert mh.moltype == expected_moltype

    mh.add_protein('AGYYG')
    mh.add_protein(u'AGYYG')
    mh.add_protein(b'AGYYG')

    if hp:
        assert len(mh.hashes) == 1
    else:
        assert len(mh.hashes) == 4


def test_protein_hp(track_abundance, hp):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 2, is_protein=True, dayhoff=False, hp=hp, track_abundance=track_abundance)
    mh.add_protein('AGYYG')

    if hp:
        assert len(mh.hashes) == 1
    else:
        assert len(mh.hashes) == 4


def test_module_translate_codon(track_abundance):
    # Ensure that translation occurs properly - module level function tests
    assert "S" == translate_codon('TCT')
    assert "S" == translate_codon('TC')
    assert "X" == translate_codon("T")

    with pytest.raises(ValueError):
        translate_codon("")
        translate_codon("TCTA")


def test_dayhoff(track_abundance):
    # verify that we can hash to dayhoff-encoded protein/aa sequences
    mh_dayhoff = MinHash(10, 2, is_protein=True,
                         dayhoff=True, hp=False, track_abundance=track_abundance)
    mh_dayhoff.add_sequence('ACTGAC')

    assert len(mh_dayhoff.hashes) == 2
    # verify that dayhoff-encoded hashes are different from protein/aa hashes
    mh_protein = MinHash(10, 2, is_protein=True, track_abundance=track_abundance)
    mh_protein.add_sequence('ACTGAC')

    assert len(mh_protein.hashes) == 2
    print(mh_protein.hashes)
    print(mh_dayhoff.hashes)
    assert mh_protein.hashes != mh_dayhoff.hashes


def test_dayhoff_2(track_abundance):
    mh = MinHash(0, 7, scaled=1, dayhoff=True, track_abundance=1)

    # first, check protein -> dayhoff hashes via minhash
    mh.add_protein('CADHIFC')
    assert len(mh) == 1
    hashval = list(mh.hashes)[0]
    assert hashval == hash_murmur('abcdefa')

    # also check seq_to_hashes
    hashes = list(mh.seq_to_hashes('CADHIFC', is_protein=True))
    assert hashval == hashes[0]

    # do we handle stop codons properly?
    mh = mh.copy_and_clear()
    mh.add_protein('CADHIF*')
    assert len(mh) == 1
    hashval = list(mh.hashes)[0]
    assert hashval == hash_murmur('abcdef*')

    # again, check seq_to_hashes
    hashes = list(mh.seq_to_hashes('CADHIF*', is_protein=True))
    assert hashval == hashes[0]


def test_hp(track_abundance):
    # verify that we can hash to hp-encoded protein/aa sequences
    mh_hp = MinHash(10, 2, is_protein=True,
                    dayhoff=False, hp=True, track_abundance=track_abundance)
    assert mh_hp.moltype == 'hp'

    mh_hp.add_sequence('ACTGAC')

    assert len(mh_hp.hashes) == 2
    # verify that hp-encoded hashes are different from protein/aa hashes
    mh_protein = MinHash(10, 2, is_protein=True, track_abundance=track_abundance)
    mh_protein.add_sequence('ACTGAC')

    assert len(mh_protein.hashes) == 2
    assert mh_protein.hashes != mh_hp.hashes


def test_hp_2(track_abundance):
    mh = MinHash(0, 3, scaled=1, hp=True, track_abundance=track_abundance)

    mh.add_protein('ANA')
    assert len(mh) == 1
    hashval = list(mh.hashes)[0]
    assert hashval == hash_murmur('hph')

    # also check seq_to_hashes
    hashes = list(mh.seq_to_hashes('ANA', is_protein=True))
    assert hashval == hashes[0]

    mh = mh.copy_and_clear()
    mh.add_protein('AN*')
    assert len(mh) == 1
    hashval = list(mh.hashes)[0]
    assert hashval == hash_murmur('hp*')

    # also check seq_to_hashes
    hashes = list(mh.seq_to_hashes('AN*', is_protein=True))
    assert hashval == hashes[0]


def test_protein_short(track_abundance):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 9, is_protein=True, track_abundance=track_abundance)
    mh.add_protein('AG')

    assert len(mh.hashes) == 0, mh.hashes


def test_size_limit(track_abundance):
    # test behavior with size limit of 3
    mh = MinHash(3, 4, track_abundance=track_abundance)
    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)
    assert list(sorted(mh.hashes)) == [10, 20, 30]
    mh.add_hash(5) # -> should push 30 off end
    assert list(sorted(mh.hashes)) == [5, 10, 20]


def test_scaled(track_abundance):
    # test behavior with scaled
    scaled = _get_scaled_for_max_hash(35)
    print('XX', scaled, _get_max_hash_for_scaled(scaled))
    mh = MinHash(0, 4, track_abundance=track_abundance, scaled=scaled)
    assert mh._max_hash == 35

    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)

    assert list(sorted(mh.hashes)) == [10, 20, 30]
    mh.add_hash(40)
    assert list(sorted(mh.hashes)) == [10, 20, 30]
    mh.add_hash(36)
    assert list(sorted(mh.hashes)) == [10, 20, 30]


def test_no_scaled(track_abundance):
    # no 'scaled', num=0 - should fail
    with pytest.raises(ValueError):
        mh = MinHash(0, 4, track_abundance=track_abundance)


def test_max_hash_conversion():
    SCALED=100000
    max_hash = _get_max_hash_for_scaled(SCALED)
    new_scaled = _get_scaled_for_max_hash(max_hash)
    assert new_scaled == SCALED


def test_max_hash_and_scaled_zero():
    max_hash = _get_max_hash_for_scaled(0)
    new_scaled = _get_scaled_for_max_hash(0)
    assert max_hash == new_scaled
    assert max_hash == 0


def test_max_hash_and_scaled_error(track_abundance):
    # test behavior when supplying both max_hash and scaled
    with pytest.raises(ValueError):
        mh = MinHash(0, 4, track_abundance=track_abundance, max_hash=35,
                     scaled=5)


def test_max_hash_cannot_limit(track_abundance):
    # make sure you can't set both n and scaled.
    with pytest.raises(ValueError):
        mh = MinHash(2, 4, track_abundance=track_abundance,
                     scaled=_get_scaled_for_max_hash(1))


def test_no_downsample_scaled_if_n(track_abundance):
    # make sure you can't set max_n and then downsample scaled
    mh = MinHash(2, 4, track_abundance=track_abundance)
    with pytest.raises(ValueError) as excinfo:
        mh.downsample(scaled=100000000)

    assert 'cannot downsample a num MinHash using scaled' in str(excinfo.value)


def test_scaled_num_both(track_abundance):
    # make sure you can't set both max_n and scaled.
    with pytest.raises(ValueError):
        mh = MinHash(2, 4, track_abundance=track_abundance, scaled=2)


def test_mh_jaccard_similarity():
    # check actual Jaccard value for a non-trivial case
    a = MinHash(0, 20, scaled=scaled50, track_abundance=False)
    b = MinHash(0, 20, scaled=scaled50, track_abundance=False)
    a.add_many([1, 3, 5, 8])
    b.add_many([1, 3, 5, 6, 8, 10])

    assert a.similarity(b) == 4. / 6.


def test_mh_similarity_downsample_jaccard_value():
    # check jaccard value after downsampling

    # max_hash = 50
    a = MinHash(0, 20, scaled=scaled50, track_abundance=False)
    # max_hash = 100
    b = MinHash(0, 20, scaled=scaled100, track_abundance=False)

    a.add_many([1, 3, 5, 8, 70])
    b.add_many([1, 3, 5, 6, 8, 10, 70 ])

    # the hash=70 will be truncated by downsampling
    assert a.similarity(b, downsample=True) == 4. / 6.


def test_mh_angular_similarity():
    # check actual angular similarity for a non-trivial case, taken from:
    # https://www.sciencedirect.com/topics/computer-science/cosine-similarity
    # note: angular similarity is 1 - 2*(acos(sim) / pi), when elements
    # are always positive (https://en.wikipedia.org/wiki/Cosine_similarity)
    a = MinHash(0, 20, scaled=scaled50, track_abundance=True)
    b = MinHash(0, 20, scaled=scaled50, track_abundance=True)
    a.set_abundances({ 1:5, 3:3, 5:2, 8:2})
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 })

    cos_sim = 0.9356
    angular_sim = 1 - 2*math.acos(cos_sim) / math.pi
    assert round(angular_sim, 4) == 0.7703

    assert round(a.similarity(b), 4) == round(angular_sim, 4)


def test_mh_angular_similarity_2():
    # check actual angular similarity for a second non-trivial case
    a = MinHash(0, 20, scaled=scaled100, track_abundance=True)
    b = MinHash(0, 20, scaled=scaled100, track_abundance=True)
    a.set_abundances({ 1:5, 3:3, 5:2, 8:2, 70:70 })
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1, 70:70 })

    assert round(a.similarity(b), 4) == 0.9728

    # ignore_abundance => jaccard
    assert a.similarity(b, ignore_abundance=True) == 5. / 7.


def test_mh_similarity_downsample_angular_value():
    # test downsample=True argument to MinHash.similarity

    # max_hash = 50
    a = MinHash(0, 20, scaled=scaled50, track_abundance=True)
    # max_hash = 100
    b = MinHash(0, 20, scaled=scaled100, track_abundance=True)

    a.set_abundances({ 1:5, 3:3, 5:2, 8:2, 70:70 })
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1, 70:70 })

    # the hash=70 will be truncated by downsampling
    sim = a.similarity(b, downsample=True)
    assert round(sim, 4) == 0.7703

    # with ignore_abundance, will be equal to jaccard
    jaccard = a.similarity(b, downsample=True, ignore_abundance=True)
    assert jaccard == 4. / 6.

def test_mh_angular_similarity_fail():
    # raise TypeError if calling angular_similarity directly and
    # one or both sketches do not have abundance info
    a = MinHash(0, 20, scaled=scaled50, track_abundance=True)
    b = MinHash(0, 20, scaled=scaled50, track_abundance=False)
    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    a.set_abundances(a_values)
    b.add_many(b_values.keys())

    # one sketch lacks track_abundance
    with pytest.raises(TypeError) as exc:
        a.angular_similarity(b)
    print(str(exc))
    assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
    # both sketches lack track abundance
    a = MinHash(0, 20, scaled=scaled50, track_abundance=False)
    a.add_many(a_values.keys())
    with pytest.raises(TypeError) as exc:
        a.angular_similarity(b)
    print(str(exc))
    assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)


def test_mh_similarity_downsample_true(track_abundance):
    # verify sim(a, b) == sim(b, a), with and without ignore_abundance

    # max_hash = 50
    a = MinHash(0, 20, scaled=scaled50, track_abundance=track_abundance)
    # max_hash = 100
    b = MinHash(0, 20, scaled=scaled100, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    # downsample=True => no error; values should match either way
    x = a.similarity(b, ignore_abundance=True, downsample=True)
    y = b.similarity(a, ignore_abundance=True, downsample=True)
    assert x == y

    # downsample=True => no error; values should match either way
    x = a.similarity(b, ignore_abundance=False, downsample=True)
    y = b.similarity(a, ignore_abundance=False, downsample=True)
    assert x == y


def test_mh_similarity_downsample_errors(track_abundance):
    # test downsample=False (default) argument to MinHash.similarity

    # max_hash = 50
    a = MinHash(0, 20, scaled=scaled50, track_abundance=track_abundance)
    # max_hash = 100
    b = MinHash(0, 20, scaled=scaled100, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    # error, incompatible max hash
    with pytest.raises(ValueError) as e:
        a.similarity(b, ignore_abundance=True)   # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        a.similarity(b, ignore_abundance=False)  # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        b.similarity(a, ignore_abundance=True)   # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        b.similarity(a, ignore_abundance=False)  # downsample=false
    assert 'mismatch in scaled; comparison fail' in str(e.value)


def test_basic_dna_bad(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(1, 4, track_abundance=track_abundance)

    with pytest.raises(ValueError) as e:
        mh.add_sequence('ATGR')
    print(e)

    assert 'invalid DNA character in input k-mer: ATGR' in str(e.value)


def test_basic_dna_bad_2(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(1, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        mh.add_protein('YYYY')


def test_basic_dna_bad_force(track_abundance):
    # test behavior on bad DNA; use 100 so multiple hashes get added.
    mh = MinHash(100, 4, track_abundance=track_abundance)
    assert len(mh.hashes) == 0
    mh.add_sequence('ATGN', True)     # ambiguous kmer skipped.
    assert len(mh.hashes) == 0
    mh.add_sequence('AATGN', True)    # but good k-mers still used.
    assert len(mh.hashes) == 1
    mh.add_sequence('AATG', True)     # checking that right kmer was added
    assert len(mh.hashes) == 1    # (only 1 hash <- this is a dup)


def test_basic_dna_bad_force_2(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(100, 4, track_abundance=track_abundance)
    assert len(mh.hashes) == 0
    mh.add_sequence('AAGNCGG', True)     # ambiguous kmers skipped.
    assert len(mh.hashes) == 0
    mh.add_sequence('AATGNGCGG', True)  # ambiguous kmers skipped.
    assert len(mh.hashes) == 2
    mh.add_sequence('AATG', True)        # checking that right kmers were added
    mh.add_sequence('GCGG', True)
    assert len(mh.hashes) == 2       # (only 2 hashes should be there)


def test_consume_lowercase(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA'.lower())
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert round(a.similarity(b), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0
    assert round(b.similarity(a), 3) == 1.0
    assert round(a.similarity(a), 3) == 1.0


def test_similarity_1(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert round(a.similarity(b), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0
    assert round(b.similarity(a), 3) == 1.0
    assert round(a.similarity(a), 3) == 1.0

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert round(a.similarity(b), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0
    assert round(b.similarity(a), 3) == 1.0
    assert round(a.similarity(a), 3) == 1.0


    b.add_sequence('GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT')
    x = a.similarity(b)
    assert x >= 0.3, x

    x = b.similarity(a)
    assert x >= 0.3, x
    assert round(a.similarity(a), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0


def test_copy(track_abundance):
    a = MinHash(20, 21, track_abundance=track_abundance)
    a.add_hash(5)
    b = a.copy()
    assert a == b
    a.add_hash(6)
    assert a != b


def test_frozen_copy(track_abundance):
    a = MinHash(20, 21, track_abundance=track_abundance)
    a.add_hash(5)
    b = a.copy()
    assert 5 in b.hashes
    a.add_hash(6)
    assert 6 not in b.hashes


def test_mh_copy(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b = a.__copy__()
    assert round(b.similarity(a), 3) == 1.0


def test_mh_len(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)

    assert len(a) == 0
    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert len(a) == 20


def test_mh_len_2(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    assert list(sorted(a.hashes)) == list(range(0, 40, 2))


def test_mh_unsigned_long_long(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    a.add_hash(9227159859419181011)        # too big for a C long int.
    assert 9227159859419181011 in a.hashes


def test_mh_count_common(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10


def test_mh_count_common_diff_protein(track_abundance):
    a = MinHash(20, 5, is_protein=False, track_abundance=track_abundance)
    b = MinHash(20, 5, is_protein=True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_maxhash(track_abundance):
    a = MinHash(0, 5, is_protein=False, track_abundance=track_abundance,
                scaled=_get_scaled_for_max_hash(1))
    b = MinHash(0, 5, is_protein=True, track_abundance=track_abundance,
                scaled=_get_scaled_for_max_hash(2))

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_seed(track_abundance):
    a = MinHash(20, 5, is_protein=False, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, is_protein=True, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_notmh(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = set()

    with pytest.raises(TypeError):
        a.count_common(b)


def test_mh_downsample_num_error(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    with pytest.raises(ValueError):
        a.downsample(num=30)


def test_mh_jaccard_asymmetric_num(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10

    # with 'jaccard', this will raise an error b/c different num
    with pytest.raises(TypeError):
        a.jaccard(b)

    a = a.downsample(num=10)
    # CTB note: this used to be 'compare', is now 'jaccard'
    assert a.jaccard(b) == 0.5
    assert b.jaccard(a) == 0.5


def test_mh_merge_typeerror(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    with pytest.raises(TypeError):
        a.merge(set())


def test_mh_merge(track_abundance):
    # test merging two identically configured minhashes
    a = MinHash(100, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(100, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c.merge(b)

    d = b.__copy__()
    d.merge(a)

    assert len(c) == len(d)
    assert list(sorted(c.hashes.items())) == list(sorted(d.hashes.items()))

    assert round(c.similarity(d), 3) == 1.0
    assert round(d.similarity(c), 3) == 1.0


def test_mh_merge_empty_num(track_abundance):
    # test merging two identically configured minhashes, one empty
    a = MinHash(100, 10, track_abundance=track_abundance)

    b = MinHash(100, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c.merge(b)

    d = b.__copy__()
    d.merge(a)

    assert len(c)
    assert len(c) == len(d)

    assert list(sorted(c.hashes.items())) == list(sorted(d.hashes.items()))
    assert round(c.similarity(d), 3) == 1.0
    assert round(d.similarity(c), 3) == 1.0


def test_mh_merge_empty_scaled(track_abundance):
    # test merging two identically configured minhashes, one empty
    a = MinHash(0, 10, scaled=1, track_abundance=track_abundance)

    b = MinHash(0, 10, scaled=1, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c.merge(b)

    d = b.__copy__()
    d.merge(a)

    assert len(c)
    assert len(c) == len(d)

    assert list(sorted(c.hashes.items())) == list(sorted(d.hashes.items()))
    assert round(c.similarity(d), 3) == 1.0
    assert round(d.similarity(c), 3) == 1.0


def test_mh_merge_check_length(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c.merge(b)
    assert len(c.hashes) == 20


def test_mh_merge_check_length2(track_abundance):
    # merged MH doesn't have full number of elements
    a = MinHash(4, 10, track_abundance=track_abundance)
    a.add_hash(3)
    a.add_hash(1)
    a.add_hash(4)

    b = MinHash(4, 10, track_abundance=track_abundance)
    b.add_hash(3)
    b.add_hash(1)
    b.add_hash(4)

    c = a.__copy__()
    c.merge(b)
    assert len(c.hashes) == 3

def test_mh_asymmetric_merge(track_abundance):
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c.merge(b)
    d = b.__copy__()
    d.merge(a)

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    # can't use jaccard on different nums without downsampling
    with pytest.raises(TypeError):
        d.jaccard(a)

    a = a.downsample(num=d.num)

    if track_abundance:
        assert round(d.similarity(a), 3) == 0.795
    else:
        assert round(d.similarity(a), 3) == 1.0

    c = c.downsample(num=b.num)
    if track_abundance:
        assert round(c.similarity(b), 3) == 0.436
    else:
        assert c.similarity(b) == 0.5


def test_mh_inplace_concat_asymmetric(track_abundance):
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b

    d = b.__copy__()
    d += a

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    try:
        d.similarity(a)
    except TypeError as exc:
        assert 'must have same num' in str(exc)

    a = a.downsample(num=d.num)
    if track_abundance:
        assert round(d.similarity(a), 3) == 0.795 # see: d += a, above.
    else:
        assert d.similarity(a) == 1.0 # see: d += a, above.

    c = c.downsample(num=b.num)
    if track_abundance:
        assert round(c.similarity(b), 3) == 0.436
    else:
        assert c.similarity(b) == 0.5


def test_mh_inplace_concat(track_abundance):
    # test merging two identically configured minhashes
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b
    d = b.__copy__()
    d += a

    assert len(c) == len(d)
    assert c.hashes == d.hashes
    assert round(c.similarity(d), 3) == 1.0
    assert round(d.similarity(c), 3) == 1.0


def test_mh_merge_diff_protein(track_abundance):
    a = MinHash(20, 5, is_protein=False, track_abundance=track_abundance)
    b = MinHash(20, 5, is_protein=True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.merge(b)


def test_mh_merge_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.merge(b)


def test_mh_similarity_diff_protein(track_abundance):
    a = MinHash(20, 5, is_protein=False, track_abundance=track_abundance)
    b = MinHash(20, 5, is_protein=True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.similarity(b)


def test_mh_similarity_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.similarity(b)


def test_mh_similarity_diff_seed(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a.similarity(b)


def test_mh_compare_diff_max_hash(track_abundance):
    a = MinHash(0, 5, track_abundance=track_abundance,
                scaled=_get_max_hash_for_scaled(5))

    b = MinHash(0, 5, track_abundance=track_abundance,
                scaled=_get_max_hash_for_scaled(10))

    with pytest.raises(ValueError):
        a.similarity(b)


def test_mh_concat_diff_protein(track_abundance):
    a = MinHash(20, 5, is_protein=False, track_abundance=track_abundance)
    b = MinHash(20, 5, is_protein=True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_max_hash(track_abundance):
    a = MinHash(0, 5, track_abundance=track_abundance,
                scaled=_get_max_hash_for_scaled(5))
    b = MinHash(0, 5, track_abundance=track_abundance,
                scaled=_get_max_hash_for_scaled(10))

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_seed(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a += b


def test_short_sequence(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    a.add_sequence('GGGG')
    # adding a short sequence should fail silently
    assert len(a.hashes) == 0


def test_bytes_murmur():
    x = hash_murmur("ACG")
    assert x == 1731421407650554201

    x = hash_murmur(b"ACG")
    assert x == 1731421407650554201

    x = hash_murmur(u"ACG")
    assert x == 1731421407650554201


def test_murmur():
    x = hash_murmur("ACG")
    assert x == 1731421407650554201

    try:
        x = hash_murmur()
        assert 0, "hash_murmur requires an argument"
    except TypeError:
        pass

    x = hash_murmur("ACG", 42)
    assert x == 1731421407650554201

    y = hash_murmur("ACG", 43)
    assert y != x


def test_abundance_simple():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_sequence('AAAAA')
    assert list(a.hashes) == [2110480117637990133]
    assert a.hashes == {2110480117637990133: 1}

    a.add_sequence('AAAAA')
    assert list(a.hashes) == [2110480117637990133]
    assert a.hashes == {2110480117637990133: 2}


def test_add_hash_with_abundance():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_hash_with_abundance(10, 1)
    assert a.hashes == {10: 1}

    a.add_hash_with_abundance(20, 2)
    assert a.hashes == {10: 1, 20: 2}

    a.add_hash_with_abundance(10, 2)
    assert a.hashes == {10: 3, 20: 2}


def test_add_hash_with_abundance_2():
    a = MinHash(20, 5, is_protein=False, track_abundance=False)

    with pytest.raises(RuntimeError) as e:
        a.add_hash_with_abundance(10, 1)

    assert "track_abundance=True when constructing" in e.value.args[0]


def test_clear():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_hash(10)
    assert a.hashes == {10: 1}

    a.clear()
    assert a.hashes == {}


def test_clear_2():
    a = MinHash(20, 5, is_protein=False, track_abundance=False)

    a.add_hash(10)
    assert list(a.hashes) == [10]

    a.clear()
    assert list(a.hashes) == []


def test_abundance_simple_2():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)
    b = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_sequence('AAAAA')
    assert list(a.hashes) == [2110480117637990133]
    assert a.hashes == {2110480117637990133: 1}

    a.add_sequence('AAAAA')
    assert list(a.hashes) == [2110480117637990133]
    assert a.hashes == {2110480117637990133: 2}

    b.add_sequence('AAAAA')
    assert a.count_common(b) == 1


def test_abundance_count_common():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)
    b = MinHash(20, 5, is_protein=False, track_abundance=False)

    a.add_sequence('AAAAA')
    a.add_sequence('AAAAA')
    assert list(a.hashes) == [2110480117637990133]
    assert a.hashes == {2110480117637990133: 2}

    b.add_sequence('AAAAA')
    b.add_sequence('GGGGG')
    assert a.count_common(b) == 1
    assert a.count_common(b) == b.count_common(a)

    assert list(sorted(b.hashes)) == [2110480117637990133, 10798773792509008305]


def test_abundance_similarity():
    a = MinHash(20, 10, track_abundance=True)
    b = MinHash(20, 10, track_abundance=False)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert round(a.similarity(b), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0
    assert round(b.similarity(a), 3) == 1.0
    assert round(a.similarity(a), 3) == 1.0

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert round(a.similarity(b), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0
    assert round(b.similarity(a), 3) == 1.0
    assert round(a.similarity(a), 3) == 1.0

    b.add_sequence('GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT')
    x = a.similarity(b)
    assert x >= 0.3, x

    x = b.similarity(a)
    assert x >= 0.3, x
    assert round(a.similarity(a), 3) == 1.0
    assert round(b.similarity(b), 3) == 1.0


def test_set_abundance():
    a = MinHash(20, 10, track_abundance=False)

    with pytest.raises(RuntimeError) as e:
        a.set_abundances({1: 3, 2: 4})

    assert "track_abundance=True when constructing" in e.value.args[0]


def test_set_abundance_2():
    datapath = utils.get_test_data("genome-s12.fa.gz.sig")
    sig = sourmash.load_one_signature(datapath,
                                      ksize=30,
                                      select_moltype='dna')
    new_mh = sig.minhash.copy_and_clear()
    mins = sig.minhash.hashes
    mins = {k: 1 for k in mins}
    new_mh.track_abundance = True
    new_mh.set_abundances(mins)

    assert set(new_mh.hashes) == set(mins)


def test_set_abundance_clear():
    # on empty minhash, clear should have no effect
    a = MinHash(20, 5, is_protein=False, track_abundance=True)
    b = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.set_abundances({1: 3, 2: 4}, clear=True)
    b.set_abundances({1: 3, 2: 4}, clear=False)

    assert list(sorted(a.hashes)) == list(sorted(b.hashes))


def test_set_abundance_clear_2():
    # default should be clear=True
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_hash(10)
    assert a.hashes == {10: 1}

    a.set_abundances({20: 2})
    assert a.hashes == {20: 2}


def test_set_abundance_clear_3():
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.add_hash(10)
    assert a.hashes == {10: 1}
    
    a.set_abundances({20: 1, 30: 4}, clear=False)
    assert a.hashes == {10: 1, 20: 1, 30: 4}


def test_set_abundance_clear_4():
    # setting the abundance of an already set hash should add
    # the abundances together
    a = MinHash(20, 5, is_protein=False, track_abundance=True)

    a.set_abundances({20: 2, 10: 1}, clear=False)   # should also sort the hashes
    assert a.hashes == {10: 1, 20: 2}

    a.set_abundances({20: 1, 10: 2}, clear=False)
    assert a.hashes == {10: 3, 20: 3}

def test_clear_abundance_on_zero():
    mh = sourmash.minhash.MinHash(n=0, ksize=31, scaled=1, track_abundance=True)
    mh.set_abundances({ 1: 5, 2: 3, 3 : 5 })
    mh.set_abundances({ 1: 0 }, clear=False)
    assert 1 not in dict(mh.hashes)
    assert dict(mh.hashes)[2] == 3
    assert dict(mh.hashes)[3] == 5
    assert len(mh) == 2

    with pytest.raises(ValueError):
        mh.set_abundances({ 2: -1 }) # Test on clear = True

    with pytest.raises(ValueError):
        mh.set_abundances({ 2: -1 }, clear=False)    
    
    assert len(mh) == 2 # Assert that nothing was affected

def test_reset_abundance_initialized():
    a = MinHash(1, 4, track_abundance=True)
    a.add_sequence('ATGC')

    # If we had a minhash with abundances and drop it, this shouldn't fail.
    # Convert from Abundance to Regular MinHash
    a.track_abundance = False

    assert list(a.hashes) == [12415348535738636339]


def test_set_abundance_initialized():
    a = MinHash(1, 4, track_abundance=False)
    a.add_sequence('ATGC')

    with pytest.raises(RuntimeError) as e:
        a.track_abundance = True

    assert "Can only set track_abundance=True if the MinHash is empty" in e.value.args[0]


def test_set_abundance_num():
    a = MinHash(2, 10, track_abundance=True)

    a.set_abundances({1: 3, 2: 4})

    assert a.hashes == {1: 3, 2: 4}


def test_mh_copy_and_clear(track_abundance):
    # test basic creation of new, empty MinHash
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = a.copy_and_clear()
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == 0
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == 0
    assert a.scaled == b.scaled
    assert b.scaled == 0


def test_mh_copy_and_clear_with_max_hash(track_abundance):
    # test basic creation of new, empty MinHash w/max_hash param set
    a = MinHash(0, 10, track_abundance=track_abundance,
                scaled=_get_scaled_for_max_hash(20))
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = a.copy_and_clear()
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == 0
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_scaled_property(track_abundance):
    scaled = 10000
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled)
    assert a.scaled == scaled


def test_pickle_protein(track_abundance):
    # check that protein/etc ksize is handled properly during serialization.
    a = MinHash(0, 10, track_abundance=track_abundance, is_protein=True,
                scaled=_get_scaled_for_max_hash(20))
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == a._max_hash
    assert b._max_hash == 20
    assert b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == len(a.hashes)
    assert len(b.hashes) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_pickle_dayhoff(track_abundance):
    # check that dayhoff ksize is handled properly during serialization.
    a = MinHash(0, 10, track_abundance=track_abundance, dayhoff=True,
                scaled=_get_scaled_for_max_hash(20))
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == a._max_hash
    assert b._max_hash == 20
    assert b.dayhoff
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == len(a.hashes)
    assert len(b.hashes) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_pickle_hp(track_abundance):
    # check that hp ksize is handled properly during serialization.
    a = MinHash(0, 10, track_abundance=track_abundance, hp=True,
                scaled=_get_scaled_for_max_hash(20))
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == a._max_hash
    assert b._max_hash == 20
    assert b.hp
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == len(a.hashes)
    assert len(b.hashes) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_pickle_max_hash(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance,
                scaled=_get_scaled_for_max_hash(20))
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == a._max_hash
    assert b._max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == len(a.hashes)
    assert len(b.hashes) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_pickle_scaled(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=922337203685477632)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b._max_hash == a._max_hash
    assert b._max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.hashes) == len(a.hashes)
    assert len(b.hashes) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_minhash_abund_add():
    # this targets part of bug #319, a segfault caused by invalidation of
    # std::vector iterators upon vector resizing - in this case, there
    # was also a bug in inserting into the middle of mins when scaled was set.

    a = MinHash(0, 10, track_abundance=True, scaled=scaled5000)

    n = 0
    for i in range(10, 0, -1):
        a.add_hash(i)
        n += 1
        assert len(a.hashes) == n
        print(len(a.hashes))


def test_minhash_abund_capacity_increase():
    # this targets bug #319, a segfault caused by invalidation of
    # std::vector iterators upon vector resizing.

    # this should set capacity to 1000 - see KmerMinHash constructor call
    # to 'reserve' when n > 0 for specific parameter.
    a = MinHash(0, 10, track_abundance=True, scaled=scaled5000)

    # 1001 is dependent on the value passed to reserve (currently 1000).
    for i in range(1001, 0, -1):
        a.add_hash(i)


def test_minhash_abund_merge_flat():
    # this targets a segfault caused by trying to compute similarity
    # of a signature with abundance and a signature without abundance.
    # the correct behavior for now is to calculate simple Jaccard,
    # i.e. 'flatten' both of them.
    a = MinHash(0, 10, track_abundance=True, scaled=scaled5000)
    b = MinHash(0, 10, scaled=scaled5000)

    for i in range(0, 10, 2):
        a.add_hash(i)

    for j in range(0, 10, 3):
        b.add_hash(i)

    # these crashed, previously.
    assert a.similarity(b) == 0.2
    assert b.similarity(a) == 0.2


def test_minhash_abund_merge_flat_2():
    # this targets a segfault caused by trying to merge
    # a signature with abundance and a signature without abundance.

    a = MinHash(0, 10, track_abundance=True, scaled=scaled5000)
    b = MinHash(0, 10, scaled=scaled5000)

    for i in range(0, 10, 2):
        a.add_hash(i)

    for j in range(0, 10, 3):
        b.add_hash(i)

    a.merge(b)


def test_distance_matrix(track_abundance):
    import numpy

    siglist = [next(signature.load_signatures(utils.get_test_data(f)))
               for f in utils.SIG_FILES]

    D1 = numpy.zeros([len(siglist), len(siglist)])
    D2 = numpy.zeros([len(siglist), len(siglist)])

    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            if i < j:
                continue
            similarity = E.similarity(E2, track_abundance)
            D2[i][j] = similarity
            D2[j][i] = similarity

    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            D1[i][j] = E.similarity(E2, track_abundance)

    assert numpy.array_equal(D1, D2)


def test_remove_many(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)

    a.add_many(list(range(0, 100, 2)))

    orig_sig = signature.SourmashSignature(a)
    orig_md5 = orig_sig.md5sum()

    a.remove_many(list(range(0, 100, 3)))
    new_sig = signature.SourmashSignature(a)
    new_md5 = new_sig.md5sum()

    assert orig_md5 == "f1cc295157374f5c07cfca5f867188a1"
    assert new_md5 == "dd93fa319ef57f4a019c59ee1a8c73e2"
    assert orig_md5 != new_md5

    assert len(a) == 33
    assert all(c % 6 != 0 for c in a.hashes)

def test_remove_minhash(track_abundance):
    original_mh = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)
    added_mh = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)
    tested_mh = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)

    original_mh.add_many(list(range(101)))
    added_mh.add_many(list(range(101,201))) # contains original in it
    tested_mh.add_many(list(range(201))) # original + added

    # Now we should expect tested_minhash == original_minhash
    # Note we are passing a MinHash object instead of an iterable object
    tested_mh.remove_many(added_mh)

    # Assertion
    original_sig = signature.SourmashSignature(original_mh)
    tested_sig = signature.SourmashSignature(tested_mh)

    # Should pass if the hashes list in the same order
    assert original_mh.hashes == tested_mh.hashes
    assert len(original_mh) == len(tested_mh)
    assert original_sig.md5sum() == tested_sig.md5sum()


def test_add_many(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)
    b = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)

    a.add_many(list(range(0, 100, 2)))
    a.add_many(list(range(0, 100, 2)))    # => abundance = 2

    assert len(a) == 50
    assert all(c % 2 == 0 for c in a.hashes)

    for h in range(0, 100, 2):
        b.add_hash(h)
        b.add_hash(h)

    assert len(b) == 50
    assert a == b


def test_set_abundances_huge():
    max_hash = 4000000
    a = MinHash(0, 10, track_abundance=True,
                scaled=_get_scaled_for_max_hash(max_hash))

    hashes = list(range(max_hash))
    abundances = itertools.repeat(2)

    a.set_abundances(dict(zip(hashes, abundances)))


def test_try_change_hashes(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)
    b = MinHash(0, 10, track_abundance=track_abundance, scaled=scaled5000)

    a.add_many(list(range(0, 100, 2)))

    h = a.hashes
    with pytest.raises(RuntimeError):
        h[5] = 10


def test_flatten():
    # test behavior with scaled
    scaled = _get_scaled_for_max_hash(35)
    mh = MinHash(0, 4, track_abundance=True, scaled=scaled)
    assert mh._max_hash == 35

    mh.add_hash(10)
    mh.add_hash(10)
    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(20)
    mh.add_hash(30)
    mh.add_hash(30)
    mh.add_hash(30)

    assert mh.hashes[10] == 3
    assert mh.hashes[20] == 2
    assert mh.hashes[30] == 3

    mh2 = mh.flatten()

    assert mh2.hashes[10] == 1
    assert mh2.hashes[20] == 1
    assert mh2.hashes[30] == 1
    assert len(mh2) == 3


def test_inflate():
    # test behavior of inflate function
    scaled = _get_scaled_for_max_hash(35)
    mh = MinHash(0, 4, track_abundance=False, scaled=scaled)
    mh2 = MinHash(0, 4, track_abundance=True, scaled=scaled)
    assert mh._max_hash == 35

    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)

    assert mh.hashes[10] == 1
    assert mh.hashes[20] == 1
    assert mh.hashes[30] == 1

    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(20)
    mh2.add_hash(20)
    mh2.add_hash(30)
    mh2.add_hash(30)
    mh2.add_hash(30)

    assert mh2.hashes[10] == 3
    assert mh2.hashes[20] == 2
    assert mh2.hashes[30] == 3

    mh3 = mh.inflate(mh2)

    assert mh3.hashes[10] == 3
    assert mh3.hashes[20] == 2
    assert mh3.hashes[30] == 3


def test_inflate_error():
    # test behavior of inflate function with 'self' as an abund sketch
    scaled = _get_scaled_for_max_hash(35)
    mh = MinHash(0, 4, track_abundance=True, scaled=scaled)
    mh2 = MinHash(0, 4, track_abundance=True, scaled=scaled)
    assert mh._max_hash == 35

    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)

    assert mh.hashes[10] == 1
    assert mh.hashes[20] == 1
    assert mh.hashes[30] == 1

    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(20)
    mh2.add_hash(20)
    mh2.add_hash(30)
    mh2.add_hash(30)
    mh2.add_hash(30)

    assert mh2.hashes[10] == 3
    assert mh2.hashes[20] == 2
    assert mh2.hashes[30] == 3

    with pytest.raises(ValueError) as exc:
        mh = mh.inflate(mh2)

    assert "inflate operates on a flat MinHash and takes a MinHash object with track_abundance=True" in str(exc.value)


def test_inflate_not_a_subset():
    # test behavior of inflate function when 'from_mh' is not a subset.
    scaled = _get_scaled_for_max_hash(35)
    mh = MinHash(0, 4, track_abundance=False, scaled=scaled)
    mh2 = MinHash(0, 4, track_abundance=True, scaled=scaled)
    assert mh._max_hash == 35

    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)

    assert mh.hashes[10] == 1
    assert mh.hashes[20] == 1
    assert mh.hashes[30] == 1

    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(10)
    mh2.add_hash(30)
    mh2.add_hash(30)
    mh2.add_hash(30)

    assert mh2.hashes[10] == 3
    assert 20 not in mh2.hashes
    assert mh2.hashes[30] == 3

    mh3 = mh.inflate(mh2)

    assert mh3.hashes[10] == 3
    assert 20 not in mh3.hashes # should intersect, in this case.
    assert mh3.hashes[30] == 3


def test_add_kmer(track_abundance):
    # test add_kmer method
    mh1 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)

    mh1.add_sequence('ATGCGTGC')
    a = mh1.hashes

    mh2.add_kmer('ATGC')
    mh2.add_kmer('TGCG')
    mh2.add_kmer('GCGT')
    mh2.add_kmer('CGTG')
    mh2.add_kmer('GTGC')
    b = mh2.hashes

    assert set(a.items()) == set(b.items())


def test_add_kmer_too_long(track_abundance):
    # test add_kmer method - should only take length k
    mh1 = MinHash(0, 4, scaled=1, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        mh1.add_kmer('ATGCGTGC')


def test_get_mins_deprecated(track_abundance):
    mh = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mins = (28945103950853965, 74690756200987412, 82962372765557409)

    mh.add_many(mins)
    mh.add_many(mins)
    mh.add_many(mins)
    mh.add_many(mins)

    with pytest.warns(DeprecationWarning):
        assert set(mh.get_mins()) == set(mins)
        if track_abundance:
            d = mh.get_mins(with_abundance=True)
            for k in mins:
                assert d[k] == 4
            assert len(d) == len(mins)


def test_get_hashes_deprecated(track_abundance):
    mh = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mins = (28945103950853965, 74690756200987412, 82962372765557409)

    mh.add_many(mins)
    mh.add_many(mins)
    mh.add_many(mins)
    mh.add_many(mins)

    with pytest.warns(DeprecationWarning):
        assert set(mh.get_hashes()) == set(mins)


def test_downsample_num(track_abundance):
    # test downsample(num=...) function
    mh = MinHash(10, 21, track_abundance=track_abundance)
    for i in range(20):
        mh.add_hash(i)

    assert mh.num == 10
    assert len(mh) == 10

    assert list(sorted(mh.hashes)) == list(range(10))

    mh2 = mh.downsample(num=5)
    assert mh2.num == 5
    assert len(mh2) == 5

    assert list(sorted(mh2.hashes)) == list(range(5))


def test_downsample_scaled(track_abundance):
    # test downsample(scaled...) method
    mh = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mins = (1, 2, 3,
            9223372036854775808 + 1, 9223372036854775808 + 2,
            9223372036854775808 + 3)
    mh.add_many(mins)

    assert len(mh) == 6
    assert list(sorted(mh.hashes)) == list(mins)

    mh2 = mh.downsample(scaled=2)
    print(mh._max_hash, mh2._max_hash)

    assert len(mh2) == 3
    assert list(sorted(mh2.hashes)) == list(mins[:3])


def test_is_molecule_type_1(track_abundance):
    mh = MinHash(1, 21, track_abundance=track_abundance)
    assert mh.moltype == 'DNA'
    assert mh.is_dna
    assert not mh.is_protein
    assert not mh.hp
    assert not mh.dayhoff


def test_is_molecule_type_2(track_abundance):
    mh = MinHash(1, 21, track_abundance=track_abundance, is_protein=True)
    assert mh.moltype == 'protein'
    assert not mh.is_dna
    assert mh.is_protein
    assert not mh.hp
    assert not mh.dayhoff


def test_is_molecule_type_3(track_abundance):
    mh = MinHash(1, 21, track_abundance=track_abundance, hp=True)
    assert mh.moltype == 'hp'
    assert not mh.is_dna
    assert not mh.is_protein
    assert mh.hp
    assert not mh.dayhoff



def test_is_molecule_type_4(track_abundance):
    mh = MinHash(1, 21, track_abundance=track_abundance, dayhoff=True)
    assert mh.moltype == 'dayhoff'
    assert not mh.is_dna
    assert not mh.is_protein
    assert not mh.hp
    assert mh.dayhoff


def test_addition_num_incompatible():
    mh1 = MinHash(10, 21)
    mh2 = MinHash(20, 21)

    mh1.add_hash(0)
    mh2.add_hash(1)

    with pytest.raises(TypeError) as exc:
        mh3 = mh1 + mh2

    assert "incompatible num values: self=10 other=20" in str(exc.value)


def test_addition_abund():
    mh1 = MinHash(10, 21, track_abundance=True)
    mh2 = MinHash(10, 21, track_abundance=True)

    mh1.set_abundances({ 0: 1 })
    mh2.set_abundances({ 0: 3 })

    mh3 = mh1 + mh2
    hashcounts = mh3.hashes
    assert len(hashcounts) == 1

    assert hashcounts[0] == 4


def test_addition_noabund():
    mh1 = MinHash(10, 21, track_abundance=False)
    mh2 = MinHash(10, 21, track_abundance=False)

    mh1.add_hash(0)
    mh2.add_hash(0)

    mh3 = mh1 + mh2
    hashcounts = mh3.hashes
    assert len(hashcounts) == 1
    assert hashcounts[0] == 1


def test_iaddition_abund():
    mh1 = MinHash(10, 21, track_abundance=True)
    mh2 = MinHash(10, 21, track_abundance=True)

    mh1.set_abundances({ 0: 1 })
    mh2.set_abundances({ 0: 3 })

    mh1 += mh2
    hashcounts = mh1.hashes
    assert len(hashcounts) == 1
    assert hashcounts[0] == 4

    hashcounts2 = mh2.hashes
    assert len(hashcounts2) == 1
    assert hashcounts2[0] == 3


def test_iaddition_noabund():
    mh1 = MinHash(10, 21, track_abundance=False)
    mh2 = MinHash(10, 21, track_abundance=False)

    mh1.add_hash(0)
    mh2.add_hash(0)

    mh1 += mh2
    hashcounts = mh1.hashes
    assert len(hashcounts) == 1
    assert hashcounts[0] == 1


def test_intersection_1_num():
    mh1 = MinHash(10, 21)
    mh2 = MinHash(10, 21)

    mh1.add_hash(0)
    mh1.add_hash(1)
    mh2.add_hash(0)
    mh2.add_hash(2)

    mh3 = mh1.intersection(mh2)
    print("mh.intersection INTERSECTION HASHES:",set(mh3.hashes))
    assert len(mh3) == 1
    assert 0 in mh3.hashes

def test_and_operator():
    mh1 = MinHash(20, 21)
    mh1.add_hash(5)
    mh1.add_hash(6)
    mh2 = MinHash(20, 21)
    mh2.add_hash(6)
    mh2.add_hash(7)

    print("\n \n mh1 EQUALS ", mh1.hashes, "\n mh2 EQUALS", mh2.hashes)

    mh3 = mh1.intersection(mh2)
    mh4 = mh1 & mh2

    print("\n Intersection hashes (mh3): ", mh3.hashes, "\n '&' hashes: (mh4)", mh4.hashes)

    assert mh3
    assert mh3 == mh4

def test_intersection_2_scaled():
    mh1 = MinHash(0, 21, scaled=1)
    mh2 = MinHash(0, 21, scaled=1)

    mh1.add_hash(0)
    mh1.add_hash(1)
    mh2.add_hash(0)
    mh2.add_hash(2)

    mh3 = mh1.intersection(mh2)
    print(set(mh3.hashes))
    assert len(mh3) == 1
    assert 0 in mh3.hashes


def test_intersection_3_abundance_error():
    # cannot intersect abundance MinHash
    mh1 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=True)

    with pytest.raises(TypeError) as exc:
        mh3 = mh1.intersection(mh2)

    assert str(exc.value) == "can only intersect flat MinHash objects"


def test_intersection_4_incompatible_ksize():
    # cannot intersect incompatible ksize etc
    mh1 = MinHash(500, 21)
    mh2 = MinHash(500, 31)

    with pytest.raises(ValueError) as exc:
        mh3 = mh1.intersection(mh2)

    assert str(exc.value) == "different ksizes cannot be compared"


def test_intersection_5_incompatible():
    # cannot intersect with non-MinHash objects
    mh1 = MinHash(0, 21, scaled=1)

    with pytest.raises(TypeError) as exc:
        mh3 = mh1.intersection(set())

    assert str(exc.value) == "can only intersect MinHash objects"


def test_intersection_6_full_num():
    # intersection of two "full" num objects is correct
    mh1 = MinHash(20, 21)
    mh2 = MinHash(20, 21)

    for i in range(100):
        mh1.add_hash(i)

    for i in range(0, 100, 2):
        mh2.add_hash(i)

    # they are both full:
    assert len(mh1) == 20
    assert len(mh2) == 20

    # intersection is symmetric:
    mh3 = mh1.intersection(mh2)
    mh4 = mh2.intersection(mh1)
    assert mh3 == mh4

    # everything in intersection is in both:
    for k in mh3.hashes:
        assert k in mh1.hashes
        assert k in mh2.hashes

    assert mh1.intersection_and_union_size(mh2) == (10, 20)

def test_intersection_7_full_scaled():
    # intersection of two scaled objects is correct
    mh1 = MinHash(0, 21, scaled=100)
    mh2 = MinHash(0, 21, scaled=100)

    for i in range(100):
        mh1.add_hash(i)

    for i in range(0, 200, 2):
        mh2.add_hash(i)

    # they both have everything:
    assert len(mh1) == 100
    assert len(mh2) == 100

    # intersection is symmetric:
    mh3 = mh1.intersection(mh2)
    mh4 = mh2.intersection(mh1)
    assert mh3 == mh4

    # everything in intersection is in both:
    for k in mh3.hashes:
        assert k in mh1.hashes
        assert k in mh2.hashes

    assert mh1.intersection_and_union_size(mh2) == (50, 150)


def test_intersection_and_union_8_incompatible_ksize():
    # cannot intersect different ksizes
    mh1 = MinHash(0, 21, scaled=1)
    mh2 = MinHash(0, 31, scaled=1)

    with pytest.raises(TypeError) as exc:
        mh1.intersection_and_union_size(mh2)
    assert "incompatible MinHash objects" in str(exc)


def test_merge_abund():
    mh1 = MinHash(10, 21, track_abundance=True)
    mh2 = MinHash(10, 21, track_abundance=True)

    mh1.set_abundances({ 0: 1 })
    mh2.set_abundances({ 0: 3 })

    ret = mh1.merge(mh2)
    assert ret is None

    print("MH1 EQUALS ", mh1.hashes)

    hashcounts = mh1.hashes
    assert len(hashcounts) == 1
    assert hashcounts[0] == 4

    hashcounts2 = mh2.hashes
    assert len(hashcounts2) == 1
    assert hashcounts2[0] == 3


def test_merge_noabund():
    mh1 = MinHash(10, 21, track_abundance=False)
    mh2 = MinHash(10, 21, track_abundance=False)

    mh1.add_hash(0)
    mh2.add_hash(0)

    ret = mh1.merge(mh2)
    assert ret is None

    hashcounts = mh1.hashes
    assert len(hashcounts) == 1
    assert hashcounts[0] == 1


def test_merge_full_num():
    # merge/union of two "full" num objects is correct
    mh1 = MinHash(20, 21)
    mh2 = MinHash(20, 21)

    for i in range(100):
        mh1.add_hash(i)

    for i in range(0, 100, 2):
        mh2.add_hash(i)

    # they are both full:
    assert len(mh1) == 20
    assert len(mh2) == 20

    # add is symmetric:
    mh3 = mh1 + mh2
    mh4 = mh2 + mh1
    assert mh3 == mh4

    # merge is full
    assert len(mh3) == 20

    # everything in union is in at least one
    for k in mh3.hashes:
        assert k in mh1.hashes or k in mh2.hashes


def test_merge_scaled():
    # merge/union of two reasonably full scaled objects is correct
    mh1 = MinHash(0, 21, scaled=100)
    mh2 = MinHash(0, 21, scaled=100)

    for i in range(100):
        mh1.add_hash(i)

    for i in range(0, 200, 2):
        mh2.add_hash(i)

    assert len(mh1) == 100
    assert len(mh2) == 100

    # merge contains all the things
    mh3 = mh1 + mh2
    assert len(mh3) == 150

    # everything in either one is in union
    for k in mh1.hashes:
        assert k in mh3.hashes
    for k in mh2.hashes:
        assert k in mh3.hashes

def test_add_is_symmetric():
    mh1 = MinHash(20, 21)
    mh1.add_hash(5)
    mh2 = MinHash(20, 21)
    mh2.add_hash(6)
    print("\n mh1 EQUALS ", mh1.hashes, "\n mh2 EQUALS", mh2.hashes)
    mh3 = mh1 + mh2
    mh4 = mh2 + mh1
    print("\n mh3 EQUALS ", mh3.hashes, "\n mh4 EQUALS", mh4.hashes)
    #if mh3 != 0, then it is "true", so it passes
    assert mh3
    assert mh3 == mh4

def test_or_equals_add():
    mh1 = MinHash(20, 21)
    mh1.add_hash(5)
    mh2 = MinHash(20, 21)
    mh2.add_hash(6)
    print("\n mh1 EQUALS ", mh1.hashes, "\n mh2 EQUALS", mh2.hashes)
    mh3 = mh1 + mh2
    mh4 = mh1 | mh2
    print("\n mh3 EQUALS ", mh3.hashes, "\n mh4 EQUALS", mh4.hashes)
    assert mh3
    assert mh3 == mh4

def test_max_containment():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 5))

    assert mh1.contained_by(mh2) == 1/4
    assert mh2.contained_by(mh1) == 1/2
    assert mh1.max_containment(mh2) == 1/2
    assert mh2.max_containment(mh1) == 1/2

def test_max_containment_debiased():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 5))

    assert mh1.contained_by_debiased(mh2) == 1/4
    assert mh2.contained_by_debiased(mh1) == 1/2
    assert mh1.max_containment_debiased(mh2) == 1/2
    assert mh2.max_containment_debiased(mh1) == 1/2


def test_max_containment_empty():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))

    assert mh1.contained_by(mh2) == 0
    assert mh2.contained_by(mh1) == 0
    assert mh1.max_containment(mh2) == 0
    assert mh2.max_containment(mh1) == 0


def test_max_containment_debiased_empty():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))

    assert mh1.contained_by_debiased(mh2) == 0
    assert mh2.contained_by_debiased(mh1) == 0
    assert mh1.max_containment_debiased(mh2) == 0
    assert mh2.max_containment_debiased(mh1) == 0


def test_max_containment_equal():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    assert mh1.contained_by(mh2) == 1
    assert mh2.contained_by(mh1) == 1
    assert mh1.max_containment(mh2) == 1
    assert mh2.max_containment(mh1) == 1


def test_max_containment_debiased_equal():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    assert mh1.contained_by_debiased(mh2) == 1
    assert mh2.contained_by_debiased(mh1) == 1
    assert mh1.max_containment_debiased(mh2) == 1
    assert mh2.max_containment_debiased(mh1) == 1


def test_avg_containment():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 5))

    assert mh1.contained_by(mh2) == 1/4
    assert mh2.contained_by(mh1) == 1/2
    assert mh1.avg_containment(mh2) == 0.375
    assert mh2.avg_containment(mh1) == 0.375


def test_avg_containment_debiased():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 5))

    assert mh1.contained_by_debiased(mh2) == 1/4
    assert mh2.contained_by_debiased(mh1) == 1/2
    assert mh1.avg_containment_debiased(mh2) == 0.375
    assert mh2.avg_containment_debiased(mh1) == 0.375


def test_avg_containment_empty():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))

    assert mh1.contained_by(mh2) == 0
    assert mh2.contained_by(mh1) == 0
    assert mh1.avg_containment(mh2) == 0
    assert mh2.avg_containment(mh1) == 0


def test_avg_containment_debiased_empty():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))

    assert mh1.contained_by_debiased(mh2) == 0
    assert mh2.contained_by_debiased(mh1) == 0
    assert mh1.avg_containment_debiased(mh2) == 0
    assert mh2.avg_containment_debiased(mh1) == 0


def test_avg_containment_equal():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    assert mh1.contained_by(mh2) == 1
    assert mh2.contained_by(mh1) == 1
    assert mh1.avg_containment(mh2) == 1
    assert mh2.avg_containment(mh1) == 1


def test_avg_containment_debiased_equal():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    assert mh1.contained_by_debiased(mh2) == 1
    assert mh2.contained_by_debiased(mh1) == 1
    assert mh1.avg_containment_debiased(mh2) == 1
    assert mh2.avg_containment_debiased(mh1) == 1


def test_frozen_and_mutable_1(track_abundance):
    # mutable minhashes -> mutable minhashes creates new copy
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = mh1.to_mutable()

    mh1.add_hash(10)
    assert 10 not in mh2.hashes


def test_frozen_and_mutable_2(track_abundance):
    # check that mutable -> frozen are separate
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh1.add_hash(10)

    mh2 = mh1.to_frozen()
    assert 10 in mh2.hashes
    mh1.add_hash(11)
    assert 11 not in mh2.hashes


def test_frozen_and_mutable_3(track_abundance):
    # check that mutable -> frozen -> mutable are all separate from each other
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh1.add_hash(10)

    mh2 = mh1.to_frozen()
    assert 10 in mh2.hashes
    mh1.add_hash(11)
    assert 11 not in mh2.hashes

    mh3 = mh2.to_mutable()
    mh3.add_hash(12)
    assert 12 not in mh2.hashes
    assert 12 not in mh1.hashes


def test_dna_kmers():
    # test seq_to_hashes for dna -> dna
    mh = MinHash(0, ksize=31, scaled=1) # DNA
    seq = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGAT"

    # first calculate seq to hashes
    hashes = mh.seq_to_hashes(seq)

    # then calculate all hashes for the sequence
    mh.add_sequence(seq)

    # identical?
    assert set(hashes) == set(mh.hashes)

    # k-mer by k-mer?
    for i in range(0, len(seq) - 31 + 1):
        # calculate each k-mer
        kmer = seq[i:i+31]

        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_sequence(kmer)
        assert len(single_mh) == 1

        # also calculate via seq_to_hashes
        hashvals = mh.seq_to_hashes(kmer)
        assert len(hashvals) == 1
        hashval = hashvals[0]

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]
        assert hashval == hashes[i]


def test_dna_kmers_2():
    # test kmers_and_hashes for dna -> dna
    mh = MinHash(0, ksize=31, scaled=1) # DNA
    seq = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGAT"

    # k-mer by k-mer?
    for kmer, hashval in mh.kmers_and_hashes(seq):
        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_sequence(kmer)
        assert len(single_mh) == 1

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]


def test_dna_kmers_3_bad_dna():
    # test kmers_and_hashes for dna -> dna, with some bad k-mers in there
    mh = MinHash(0, ksize=31, scaled=1) # DNA
    seq = "NTGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGAT"

    with pytest.raises(ValueError) as exc:
        list(mh.kmers_and_hashes(seq))

    assert "invalid DNA character in input k-mer: NTGCGAGTGT" in str(exc)


def test_dna_kmers_4_bad_dna():
    # test kmers_and_hashes for bad dna -> dna, using force
    mh = MinHash(0, ksize=31, scaled=1) # DNA
    seq = "NTGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGAT"

    # k-mer by k-mer?
    found_bad_kmer = False
    for kmer, hashval in mh.kmers_and_hashes(seq, force=True):
        # add to minhash obj
        single_mh = mh.copy_and_clear()

        if hashval == None:
            assert kmer == seq[:31] # first k-mer is baaaaad.
            found_bad_kmer = True
            continue

        # if the below code raises an exception, it's because the above
        # 'if' statement was not triggered (but should have been :)
        single_mh.add_sequence(kmer)
        assert len(single_mh) == 1
        assert hashval == list(single_mh.hashes)[0]

    assert found_bad_kmer, "there is one bad k-mer in here"


def test_protein_kmers():
    # test seq_to_hashes for protein -> protein
    mh = MinHash(0, ksize=7, is_protein=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # first calculate seq to hashes
    hashes = mh.seq_to_hashes(seq, is_protein=True)

    # then calculate all hashes for the sequence
    mh.add_protein(seq)

    # identical?
    assert set(hashes) == set(mh.hashes)

    # k-mer by k-mer?
    for i in range(0, len(seq) - 7 + 1):
        # calculate each k-mer
        kmer = seq[i:i+7]

        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # also calculate via seq_to_hashes
        hashvals = mh.seq_to_hashes(kmer, is_protein=True)
        assert len(hashvals) == 1
        hashval = hashvals[0]

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]
        assert hashval == hashes[i]


def test_protein_kmers_2():
    # test kmers_and_hashes for protein -> protein
    mh = MinHash(0, ksize=7, is_protein=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # k-mer by k-mer?
    for kmer, hashval in mh.kmers_and_hashes(seq, is_protein=True):
        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]


def test_dayhoff_kmers():
    # test seq_to_hashes for protein -> dayhoff
    mh = MinHash(0, ksize=7, dayhoff=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # first calculate seq to hashes
    hashes = mh.seq_to_hashes(seq, is_protein=True)

    # then calculate all hashes for the sequence
    mh.add_protein(seq)

    # identical?
    assert set(hashes) == set(mh.hashes)

    # k-mer by k-mer?
    for i in range(0, len(seq) - 7 + 1):
        # calculate each k-mer
        kmer = seq[i:i+7]

        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # also calculate via seq_to_hashes
        hashvals = mh.seq_to_hashes(kmer, is_protein=True)
        assert len(hashvals) == 1
        hashval = hashvals[0]

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]
        assert hashval == hashes[i]


def test_dayhoff_kmers_2():
    # test kmers_and_hashes for protein -> dayhoff
    mh = MinHash(0, ksize=7, dayhoff=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # k-mer by k-mer?
    for kmer, hashval in mh.kmers_and_hashes(seq, is_protein=True):
        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]


def test_hp_kmers():
    # test hashes_to_seq for protein -> hp
    mh = MinHash(0, ksize=7, hp=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # first calculate seq to hashes
    hashes = mh.seq_to_hashes(seq, is_protein=True)

    # then calculate all hashes for the sequence
    mh.add_protein(seq)

    # identical?
    assert set(hashes) == set(mh.hashes)

    # k-mer by k-mer?
    for i in range(0, len(seq) - 7 + 1):
        # calculate each k-mer
        kmer = seq[i:i+7]

        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # also calculate via seq_to_hashes
        hashvals = mh.seq_to_hashes(kmer, is_protein=True)
        assert len(hashvals) == 1
        hashval = hashvals[0]

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]
        assert hashval == hashes[i]


def test_hp_kmers_2():
    # test kmers_and_hashes for protein -> hp
    mh = MinHash(0, ksize=7, hp=True, scaled=1)
    seq = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    # k-mer by k-mer?
    for kmer, hashval in mh.kmers_and_hashes(seq, is_protein=True):
        # add to minhash obj
        single_mh = mh.copy_and_clear()
        single_mh.add_protein(kmer)
        assert len(single_mh) == 1

        # confirm it all matches
        assert hashval == list(single_mh.hashes)[0]


def test_translate_protein_hashes():
    # test seq_to_hashes for dna -> protein
    mh = MinHash(0, ksize=7, is_protein=True, scaled=1)
    mh_translate = mh.copy()
    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa"
    prot = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    hashes_translate = mh_translate.seq_to_hashes(dna)
    hashes_prot = mh_translate.seq_to_hashes(prot, is_protein=True)

    # one is a subset of other, b/c of six frame translation
    assert set(hashes_prot).issubset(set(hashes_translate))
    assert not set(hashes_translate).issubset(set(hashes_prot))


def test_translate_protein_hashes_2():
    # test kmers_and_hashes for dna -> protein
    mh_translate = MinHash(0, ksize=7, is_protein=True, scaled=1)

    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa".upper()

    # does everything match? check!
    k_and_h = list(mh_translate.kmers_and_hashes(dna))
    for idx, kmer in enumerate(_kmers_from_all_coding_frames(dna, 21)):
        k, h = k_and_h[idx]

        assert kmer == k
        assert _hash_fwd_only(mh_translate, kmer) == h


def test_translate_hp_hashes():
    # test seq_to_hashes for dna -> protein -> hp
    mh = MinHash(0, ksize=7, hp=True, scaled=1)
    mh_translate = mh.copy()
    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa"
    prot = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    hashes_translate = mh_translate.seq_to_hashes(dna)
    hashes_prot = mh_translate.seq_to_hashes(prot, is_protein=True)

    # one is a subset of other, b/c of six frame translation
    assert set(hashes_prot).issubset(set(hashes_translate))
    assert not set(hashes_translate).issubset(set(hashes_prot))


def test_translate_hp_hashes_2():
    # test kmers_and_hashes for dna -> protein -> hp
    mh_translate = MinHash(0, ksize=7, hp=True, scaled=1)
    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa"
    dna = dna.upper()

    # does everything match? check!
    k_and_h = list(mh_translate.kmers_and_hashes(dna))
    for idx, kmer in enumerate(_kmers_from_all_coding_frames(dna, 21)):
        k, h = k_and_h[idx]

        assert kmer == k
        assert _hash_fwd_only(mh_translate, kmer) == h


def test_translate_dayhoff_hashes():
    # test seq_to_hashes for dna -> protein -> dayhoff
    mh = MinHash(0, ksize=7, dayhoff=True, scaled=1)
    mh_translate = mh.copy()
    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa"
    prot = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"

    hashes_translate = mh_translate.seq_to_hashes(dna)
    hashes_prot = mh_translate.seq_to_hashes(prot, is_protein=True)

    # one is a subset of other, b/c of six frame translation
    assert set(hashes_prot).issubset(set(hashes_translate))
    assert not set(hashes_translate).issubset(set(hashes_prot))


def test_translate_dayhoff_hashes_2():
    # test kmers_and_hashes for dna -> protein -> dayhoff
    mh_translate = MinHash(0, ksize=7, dayhoff=True, scaled=1)
    dna = "atggttaaagtttatgccccggcttccagtgccaatatgagcgtcgggtttgatgtgctcggggcggcggtgacacctgttgatggtgcattgctcggagatgtagtcacggttgaggcggcagagacattcagtctcaacaacctcggacgctttgccgataagctgccgtcagaaccacgggaaaatatcgtttatcagtgctgggagcgtttttgccaggaactgggtaagcaaattccagtggcgatgaccctggaaaagaatatgccgatcggttcgggcttaggctccagtgcctgttcggtggtcgcggcgctgatggcgatgaatgaacactgcggcaagccgcttaatgacactcgtttgctggctttgatgggcgagctggaaggccgtatctccggcagcattcattacgacaacgtggcaccgtgttttctcggtggtatgcagttgatgatcgaagaaaacgacatcatcagccagcaagtgccagggtttgatgagtggctgtgggtgctggcgtatccggggattaaagtctcgacggcagaagccagggctattttaccggcgcagtatcgccgccaggattgcattgcgcacgggcgacatctggcaggcttcattcacgcctgctattcccgtcagcctgagcttgccgcgaagctgatgaaagatgttatcgctgaaccctaccgtgaacggttactgccaggcttccggcaggcgcggcaggcggtcgcggaaatcggcgcggtagcgagcggtatctccggctccggcccgaccttgttcgctctgtgtgacaagccggaaaccgcccagcgcgttgccgactggttgggtaagaactacctgcaaaatcaggaaggttttgttcatatttgccggctggatacggcgggcgcacgagtactggaaaactaa"
    dna = dna.upper()

    # does everything match? check!
    k_and_h = list(mh_translate.kmers_and_hashes(dna))
    for idx, kmer in enumerate(_kmers_from_all_coding_frames(dna, 21)):
        k, h = k_and_h[idx]

        assert kmer == k
        assert _hash_fwd_only(mh_translate, kmer) == h


def test_containment(track_abundance):
    "basic containment test. note: containment w/abundance ignores abundance."
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    assert mh1.contained_by(mh2) == 1/4
    assert mh2.contained_by(mh1) == 1/2


def test_containment_debiased(track_abundance):
    "basic containment test. note: containment w/abundance ignores abundance."
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    assert mh1.contained_by_debiased(mh2) == 1/4
    assert mh2.contained_by_debiased(mh1) == 1/2


def test_sum_abundances(track_abundance):
    "test sum_abundances"
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    if track_abundance:
        assert mh1.sum_abundances == 6
        assert mh2.sum_abundances == 6
    else:
        assert mh1.sum_abundances == None
        assert mh2.sum_abundances == None


def test_mean_abundance(track_abundance):
    "test mean_abundance"
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    if track_abundance:
        assert mh1.mean_abundance == 1.5
        assert mh2.mean_abundance == 3
    else:
        assert not mh1.mean_abundance
        assert not mh2.mean_abundance


def test_median_abundance(track_abundance):
    "test median_abundance"
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    if track_abundance:
        assert mh1.median_abundance == 1.5
        assert mh2.median_abundance == 3
    else:
        assert not mh1.median_abundance
        assert not mh2.median_abundance


def test_std_abundance(track_abundance):
    "test std_abundance"
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))
    mh2.add_many((1, 5))

    if track_abundance:
        assert mh1.std_abundance == 0.5
        assert mh2.std_abundance == 0.0
    else:
        assert not mh1.std_abundance
        assert not mh2.std_abundance


def test_unique_dataset_hashes(track_abundance):
    "test total_hashes approximation"
    mh1 = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    mh2 = MinHash(4, 21, track_abundance=track_abundance)

    mh1.add_many((1, 2, 3, 4))
    mh1.add_many((1, 2))
    mh2.add_many((1, 5))

    assert mh1.unique_dataset_hashes == 4
    with pytest.raises(TypeError) as exc:
        mh2.unique_dataset_hashes
    assert "can only approximate unique_dataset_hashes for scaled MinHashes" in str(exc)


def test_containment_ANI():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash

    m1_cont_m2 = mh1.containment_ani(mh2, estimate_ci =True)
    m2_cont_m1 = mh2.containment_ani(mh1, estimate_ci =True)
    print("\nmh1 contained by mh2", m1_cont_m2)
    print("mh2 contained by mh1", m2_cont_m1)

    assert (round(m1_cont_m2.ani,3), m1_cont_m2.ani_low, m1_cont_m2.ani_high) == (1.0, 1.0, 1.0)
    assert (round(m2_cont_m1.ani,3), round(m2_cont_m1.ani_low,3), round(m2_cont_m1.ani_high,3)) == (0.966, 0.965, 0.967)

    m1_mc_m2 = mh1.max_containment_ani(mh2, estimate_ci =True)
    m2_mc_m1 = mh2.max_containment_ani(mh1, estimate_ci =True)
    print("mh1 max containment", m1_mc_m2)
    print("mh2 max containment", m2_mc_m1)
    m1_mc_m2.size_is_inaccurate = False
    m2_mc_m1.size_is_inaccurate = False
    assert m1_mc_m2 == m2_mc_m1
    assert (round(m1_mc_m2.ani, 3), round(m1_mc_m2.ani_low, 3), round(m1_mc_m2.ani_high, 3)) == (1.0,1.0,1.0)


def test_containment_ANI_precalc_containment():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash
    # precalc containments and assert same results
    s1c = mh1.contained_by(mh2)
    s2c = mh2.contained_by(mh1)
    mc = max(s1c, s2c)

    assert mh1.containment_ani(mh2, estimate_ci=True) ==  mh1.containment_ani(mh2, containment=s1c, estimate_ci=True)
    assert mh2.containment_ani(mh1) ==  mh2.containment_ani(mh1, containment=s2c)
    assert mh1.max_containment_ani(mh2) ==  mh2.max_containment_ani(mh1)
    assert mh1.max_containment_ani(mh2) ==  mh1.max_containment_ani(mh2, max_containment=mc)
    assert mh1.max_containment_ani(mh2) ==  mh2.max_containment_ani(mh1, max_containment=mc)


def test_avg_containment_ani():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash
    # check average_containment_ani
    ac_m1 = mh1.avg_containment_ani(mh2)
    ac_m2 = mh2.avg_containment_ani(mh1)
    assert ac_m1 == ac_m2 == (mh1.containment_ani(mh2).ani + mh2.containment_ani(mh1).ani)/2 


def test_containment_ANI_downsample():
    f2 = utils.get_test_data('2+63.fa.sig')
    f3 = utils.get_test_data('47+63.fa.sig')
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash
    mh3 = sourmash.load_one_signature(f3, ksize=31).minhash
    # check that downsampling works properly
    print(mh2.scaled)
    mh2 = mh2.downsample(scaled=1100)
    assert mh2.scaled != mh3.scaled
    ds_s3c = mh2.containment_ani(mh3, downsample=True)
    ds_s4c = mh3.containment_ani(mh2, downsample=True)
    mc_w_ds_1 =  mh2.max_containment_ani(mh3, downsample=True)
    mc_w_ds_2 =  mh3.max_containment_ani(mh2, downsample=True)
    print(ds_s3c)
    with pytest.raises(ValueError) as e:
        mh2.containment_ani(mh3)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    with pytest.raises(ValueError) as e:
        mh2.max_containment_ani(mh3)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    mh3 = mh3.downsample(scaled=1100)
    assert mh2.scaled == mh3.scaled
    ds_s3c_manual = mh2.containment_ani(mh3)
    ds_s4c_manual = mh3.containment_ani(mh2)
    ds_mc_manual =  mh2.max_containment_ani(mh3)
    assert ds_s3c == ds_s3c_manual
    assert ds_s4c == ds_s4c_manual
    assert mc_w_ds_1 == mc_w_ds_2 == ds_mc_manual

    ac_m2 = mh2.avg_containment_ani(mh3)
    ac_m3 = mh3.avg_containment_ani(mh2)
    assert ac_m2 == ac_m3 == (ds_s3c.ani + ds_s4c.ani)/2


def test_jaccard_ANI():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash

    print("\nJACCARD_ANI", mh1.jaccard_ani(mh2))

    m1_jani_m2 = mh1.jaccard_ani(mh2)
    m2_jani_m1 = mh2.jaccard_ani(mh1)

    assert m1_jani_m2 == m2_jani_m1
    assert (m1_jani_m2.ani, m1_jani_m2.p_nothing_in_common, m1_jani_m2.jaccard_error) == (0.9783711630110239, 0.0, 3.891666770716877e-07)


def test_jaccard_ANI_untrustworthy():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash

    print("\nJACCARD_ANI", mh1.jaccard_ani(mh2))

    m1_jani_m2 = mh1.jaccard_ani(mh2, err_threshold=1e-7)

    # since size is inaccurate on 2.fa.sig, need to override to be able to get ani
    m1_jani_m2.size_is_inaccurate = False

    assert m1_jani_m2.ani == None
    assert m1_jani_m2.je_exceeds_threshold==True
    assert m1_jani_m2.je_threshold == 1e-7


def test_jaccard_ANI_precalc_jaccard():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash
    # precalc jaccard and assert same result
    jaccard = mh1.jaccard(mh2)
    print("\nJACCARD_ANI", mh1.jaccard_ani(mh2,jaccard=jaccard))

    assert mh1.jaccard_ani(mh2) == mh1.jaccard_ani(mh2, jaccard=jaccard) == mh2.jaccard_ani(mh1, jaccard=jaccard)
    wrong_jaccard = jaccard - 0.1
    assert mh1.jaccard_ani(mh2) != mh1.jaccard_ani(mh2, jaccard=wrong_jaccard)


def test_jaccard_ANI_downsample():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash

    print(mh1.scaled)
    mh1 = mh1.downsample(scaled=2000)
    assert mh1.scaled != mh2.scaled
    with pytest.raises(ValueError) as e:
        mh1.jaccard_ani(mh2)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    ds_s1c = mh1.jaccard_ani(mh2, downsample=True)
    ds_s2c = mh2.jaccard_ani(mh1, downsample=True)

    mh2 = mh2.downsample(scaled=2000)
    assert mh1.scaled == mh2.scaled
    ds_j_manual = mh1.jaccard_ani(mh2)
    assert ds_s1c == ds_s2c == ds_j_manual


def test_containment_ani_ci_tiny_testdata():
    """
    tiny test data to trigger the following:
    WARNING: Cannot estimate ANI confidence intervals from containment. Do your sketches contain enough hashes?
    Error: varN <0.0!
    """
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    m2_cani_m1 = mh2.containment_ani(mh1, estimate_ci=True)
    print(m2_cani_m1)
    # from the formula ANI = c^(1/k) for c=3/4 and k=21
    np.testing.assert_almost_equal(m2_cani_m1.ani, 0.986394259982259, decimal=3)
    m2_cani_m1.size_is_inaccurate = False
    assert m2_cani_m1.ani_low == None
    assert m2_cani_m1.ani_high == None


def test_containment_num_fail():
    f1 = utils.get_test_data('num/47.fa.sig')
    f2 = utils.get_test_data('num/63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash

    with pytest.raises(TypeError) as exc:
        mh1.contained_by(mh2)
    print(str(exc))
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh2.contained_by_debiased(mh1)
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.max_containment(mh2)
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.max_containment_debiased(mh2)
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.avg_containment(mh2)
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.avg_containment_debiased(mh2)
    assert "Error: can only calculate containment for scaled MinHashes" in str(exc)


def test_ANI_num_fail():
    f1 = utils.get_test_data('num/47.fa.sig')
    f2 = utils.get_test_data('num/63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2, ksize=31).minhash

    with pytest.raises(TypeError) as exc:
        mh1.containment_ani(mh2)
    print(str(exc))
    assert "Error: can only calculate ANI for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh2.containment_ani(mh1, estimate_ci =True)
    assert "Error: can only calculate ANI for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.max_containment_ani(mh2)
    assert "Error: can only calculate ANI for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.avg_containment_ani(mh2)
    assert "Error: can only calculate ANI for scaled MinHashes" in str(exc)
    with pytest.raises(TypeError) as exc:
        mh1.jaccard_ani(mh2)
    assert "Error: can only calculate ANI for scaled MinHashes" in str(exc)


def test_minhash_set_size_estimate_is_accurate():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash
    mh1_ds = mh1.downsample(scaled=100000)
    # check accuracy using default thresholds (rel_err= 0.2, confidence=0.95)
    assert mh1.size_is_accurate() == True
    assert mh1_ds.size_is_accurate() == False
    assert mh2.size_is_accurate() == True

    # change rel err
    assert mh1.size_is_accurate(relative_error=0.5) == True
    assert mh2.size_is_accurate(relative_error=0.0001) == False

    # change prob
    assert mh1.size_is_accurate(confidence=0.5) == True
    assert mh1.size_is_accurate(relative_error=0.001, confidence=1) == False

    # check that relative error and confidence must be between 0 and 1
    with pytest.raises(ValueError) as exc:
        mh2.size_is_accurate(relative_error=-1)
    assert "Error: relative error and confidence values must be between 0 and 1." in str(exc)

    with pytest.raises(ValueError) as exc:
        mh2.size_is_accurate(confidence=-1)
    assert "Error: relative error and confidence values must be between 0 and 1." in str(exc)

    with pytest.raises(ValueError) as exc:
        mh2.size_is_accurate(relative_error=-1, confidence=-1)
    assert "Error: relative error and confidence values must be between 0 and 1." in str(exc)


def test_minhash_ani_inaccurate_size_est():
    # TODO: It's actually really tricky to get the set size to be inaccurate. Eg. For a scale factor of 10000,
    # you would need
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash
    mh2 = sourmash.load_one_signature(f2).minhash
    # downsample
    mh1_ds = mh1.downsample(scaled=100000)
    mh2_ds = mh2.downsample(scaled=100000)
    assert mh1.size_is_accurate(relative_error=0.05, confidence=0.95) == True
    assert mh1.size_is_accurate() == True
    assert mh1_ds.size_is_accurate() == False
    assert mh2.size_is_accurate() == True

    assert round(mh1.jaccard_ani(mh2).ani, 3) == 0.978

    m2_ca_m1 = mh2.containment_ani(mh1)
    assert round(m2_ca_m1.ani, 3) == 0.966
    assert m2_ca_m1.size_is_inaccurate == False

    m1_ca_m2_ds = mh1_ds.containment_ani(mh2_ds)
    print(m1_ca_m2_ds)
    assert m1_ca_m2_ds.ani == None #0.987
    assert m1_ca_m2_ds.size_is_inaccurate == True


def test_size_num_fail():
    f1 = utils.get_test_data('num/47.fa.sig')
    mh1 = sourmash.load_one_signature(f1, ksize=31).minhash

    with pytest.raises(TypeError) as exc:
        mh1.size_is_accurate()
    print(str(exc))
    assert "Error: can only estimate dataset size for scaled MinHashes" in str(exc)
