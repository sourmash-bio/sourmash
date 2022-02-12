"""
Legacy tests from when there were Estimator objects and not just MinHash
objects.
"""

import pytest
from sourmash import MinHash

import sourmash_tst_utils as utils

# below, 'track_abundance' is toggled to both True and False by py.test --
# see conftest.py.


def test_jaccard_1(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.add_hash(i)

    # here the union is [1, 2, 3, 4, 5]
    # and the intersection is [1, 2, 3, 4] => 4/5.

    assert round(E1.jaccard(E2), 2) == round(4 / 5.0, 2)
    assert round(E2.jaccard(E1), 2) == round(4 / 5.0, 2)


def test_jaccard_2_difflen(track_abundance):
    E1 = MinHash(n=5, ksize=20, track_abundance=track_abundance)
    E2 = MinHash(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.add_hash(i)
    for i in [1, 2, 3, 4]:
        E2.add_hash(i)

    print(E1.jaccard(E2))
    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 5.0


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
        e2.add_kmer(seq[i:i + 4])

    assert e1.hashes.keys() == e2.hashes.keys()
    print(e1.hashes.keys())
    assert 726311917625663847 in e1.hashes.keys()
    assert 3697418565283905118 in e1.hashes.keys()


def test_protein_mh(track_abundance):
    e1 = MinHash(n=5, ksize=2, is_protein=True,
                    track_abundance=track_abundance)
    e2 = MinHash(n=5, ksize=2, is_protein=True,
                    track_abundance=track_abundance)

    # ok, so this is confusing, but: we are adding _DNA_ kmers here,
    # and translating. so, add_sequence and add_kmer actually both add
    # 6-mers.
    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    for i in range(len(seq) - 5):
        kmer = seq[i:i + 6]
        e2.add_kmer(kmer)

    assert e1.hashes.keys() == e2.hashes.keys()
    assert 901193879228338100 in e1.hashes.keys()


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

    assert e1.hashes == e2.hashes
    assert e1.num == e2.num
    assert e1.ksize == e2.ksize
    assert e1.is_protein == e2.is_protein
    assert e1.scaled == e2.scaled
    assert e1.scaled == 0
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
    assert round(E1.similarity(E2, ignore_abundance=True), 2) == 0.5


def test_abund_similarity_zero():
    E1 = MinHash(n=5, ksize=20, track_abundance=True)
    E2 = MinHash(n=5, ksize=20, track_abundance=True)

    for i in [1]:
        E1.add_hash(i)

    assert E1.similarity(E2) == 0.0


####

def test_jaccard_on_real_data():
    from sourmash.signature import load_signatures

    afile = 'n10000/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz'
    a = utils.get_test_data(afile)
    sig1 = list(load_signatures(a))[0]
    mh1 = sig1.minhash

    bfile = 'n10000/GCF_000006945.1_ASM694v1_genomic.fna.gz.sig.gz'
    b = utils.get_test_data(bfile)
    sig2 = list(load_signatures(b))[0]
    mh2 = sig2.minhash

    assert mh1.similarity(mh2) == 0.0183
    assert mh2.similarity(mh1) == 0.0183

    mh1 = mh1.downsample(num=1000)
    mh2 = mh2.downsample(num=1000)
    assert mh1.similarity(mh2) == 0.011
    assert mh2.similarity(mh1) == 0.011

    mh1 = mh1.downsample(num=100)
    mh2 = mh2.downsample(num=100)
    assert mh1.similarity(mh2) == 0.01
    assert mh2.similarity(mh1) == 0.01

    mh1 = mh1.downsample(num=10)
    mh2 = mh2.downsample(num=10)
    assert mh1.similarity(mh2) == 0.0
    assert mh2.similarity(mh1) == 0.0


def test_scaled_on_real_data():
    from sourmash.signature import load_signatures

    afile = 'scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz'
    a = utils.get_test_data(afile)
    sig1 = list(load_signatures(a))[0]
    mh1 = sig1.minhash

    bfile = 'scaled100/GCF_000006945.1_ASM694v1_genomic.fna.gz.sig.gz'
    b = utils.get_test_data(bfile)
    sig2 = list(load_signatures(b))[0]
    mh2 = sig2.minhash

    assert round(mh1.similarity(mh2), 5) == 0.01644
    assert round(mh2.similarity(mh1), 5) == 0.01644

    mh1 = mh1.downsample(scaled=100)
    mh2 = mh2.downsample(scaled=100)
    assert round(mh1.similarity(mh2), 5) == 0.01644
    assert round(mh2.similarity(mh1), 5) == 0.01644

    mh1 = mh1.downsample(scaled=1000)
    mh2 = mh2.downsample(scaled=1000)
    assert round(mh1.similarity(mh2), 5) == 0.01874
    assert round(mh2.similarity(mh1), 5) == 0.01874

    mh1 = mh1.downsample(scaled=10000)
    mh2 = mh2.downsample(scaled=10000)

    assert mh1.similarity(mh2) == 0.01
    assert mh2.similarity(mh1) == 0.01


def test_scaled_on_real_data_2():
    from sourmash.signature import load_signatures

    afile = 'scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz'
    a = utils.get_test_data(afile)
    sig1 = list(load_signatures(a))[0]
    mh1 = sig1.minhash

    bfile = 'scaled100/GCF_000006945.1_ASM694v1_genomic.fna.gz.sig.gz'
    b = utils.get_test_data(bfile)
    sig2 = list(load_signatures(b))[0]
    mh2 = sig2.minhash

    assert round(mh1.similarity(mh2), 5) == 0.01644
    assert round(mh2.similarity(mh1), 5) == 0.01644

    mh1 = mh1.downsample(scaled=1000)
    mh2 = mh2.downsample(scaled=1000)

    assert round(mh1.similarity(mh2), 4) == 0.0187
    assert round(mh2.similarity(mh1), 4) == 0.0187

    mh1 = mh1.downsample(scaled=10000)
    mh2 = mh2.downsample(scaled=10000)
    assert round(mh1.similarity(mh2), 3) == 0.01
    assert round(mh2.similarity(mh1), 3) == 0.01

    mh1 = mh1.downsample(scaled=100000)
    mh2 = mh2.downsample(scaled=100000)
    assert round(mh1.similarity(mh2), 2) == 0.01
    assert round(mh2.similarity(mh1), 2) == 0.01


def test_downsample_scaled_with_num():
    from sourmash.signature import load_signatures

    afile = 'scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz'
    a = utils.get_test_data(afile)
    sig1 = list(load_signatures(a))[0]
    mh1 = sig1.minhash

    with pytest.raises(ValueError) as exc:
        mh = mh1.downsample(num=500)

    assert 'cannot downsample a scaled MinHash using num' in str(exc.value)
