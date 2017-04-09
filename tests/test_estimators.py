"""
Legacy tests from when there were Estimator objects and not just MinHash
objects.
"""

from __future__ import print_function, unicode_literals

import pytest
from sourmash_lib import MinHash

# below, 'track_abundance' is toggled to both True and False by py.test --
# see conftest.py.


def test_jaccard_1(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 5.0


def test_jaccard_2_difflen(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4]:
        E2.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 4.0


def test_common_1(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.add_hash(i)

    assert E1.count_common(E2) == 4
    assert E2.count_common(E1) == 4


def test_diff_seed(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance, seed=1)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance, seed=2)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.add_hash(i)

    with pytest.raises(ValueError):
        E1.count_common(E2)


def test_dna_mh(track_abundance):
    e1 = MinHash(n=5, ksize=4, track_abundance=track_abundance)
    e2 = MinHash(n=5, ksize=4, track_abundance=track_abundance)

    seq = 'ATGGCAGTGACGATGCCAG'
    e1.add_sequence(seq)
    for i in range(len(seq) - 3):
        e2.add(seq[i:i + 4])

    assert e1.get_mins() == e2.get_mins()
    print(e1.get_mins())
    assert 726311917625663847 in e1.get_mins()
    assert 3697418565283905118 in e1.get_mins()


def test_protein_mh(track_abundance):
    e1 = MinHash(n=5, ksize=6, is_protein=True,
                    track_abundance=track_abundance)
    e2 = MinHash(n=5, ksize=6, is_protein=True,
                    track_abundance=track_abundance)

    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    for i in range(len(seq) - 5):
        kmer = seq[i:i + 6]
        e2.add(kmer)

    assert e1.get_mins() == e2.get_mins()
    assert 901193879228338100 in e1.get_mins()


def test_pickle(track_abundance):
    import pickle
    from io import BytesIO

    e1 = MinHash(n=5, ksize=6, is_protein=False,
                 track_abundance=track_abundance)

    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)
    e1.add_sequence(seq)

    fp = BytesIO()
    pickle.dump(e1, fp)

    fp2 = BytesIO(fp.getvalue())
    e2 = pickle.load(fp2)

    assert e1.get_mins(with_abundance=track_abundance) == \
           e2.get_mins(with_abundance=track_abundance)
    assert e1.num == e2.num
    assert e1.ksize == e2.ksize
    assert e1.is_protein == e2.is_protein
    assert e1.max_hash == e2.max_hash
    assert e1.seed == e2.seed


def test_bad_construct_1(track_abundance):
    try:
        e1 = MinHash(ksize=6, is_protein=False,
                        track_abundance=track_abundance)
        assert 0, "require n in constructor"
    except TypeError:
        pass


def test_bad_construct_2(track_abundance):
    try:
        e1 = MinHash(n=100, is_protein=False,
                        track_abundance=track_abundance)
        assert 0, "require ksize in constructor"
    except TypeError:
        pass


def test_abund_similarity():
    E1 = MinHash(n=5, ksize=20, track_abundance=True)
    E2 = MinHash(n=5, ksize=20, track_abundance=True)

    for i in [1]:
        E1.add_hash(i)
    for i in [1, 2]:
        E2.add_hash(i)

    assert round(E1.similarity(E1)) == 1.0
    assert round(E1.similarity(E2), 2) == 0.5

    assert round(E1.similarity(E1, ignore_abundance=True)) == 1.0
    assert round(E1.similarity(E2, ignore_abundance=True), 2) == 1.0


def test_abund_similarity_zero():
    E1 = MinHash(n=5, ksize=20, track_abundance=True)
    E2 = MinHash(n=5, ksize=20, track_abundance=True)

    for i in [1]:
        E1.add_hash(i)

    assert E1.similarity(E2) == 0.0
