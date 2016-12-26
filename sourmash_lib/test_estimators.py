from . import Estimators

# below, 'track_abundance' is toggled to both True and False by py.test --
# see conftest.py.


def test_jaccard_1(track_abundance):
    E1 = Estimators(n=5, ksize=20, track_abundance=track_abundance)
    E2 = Estimators(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.mh.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 5.0


def test_jaccard_2_difflen(track_abundance):
    E1 = Estimators(n=5, ksize=20, track_abundance=track_abundance)
    E2 = Estimators(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4]:
        E2.mh.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 4.0


def test_common_1(track_abundance):
    E1 = Estimators(n=5, ksize=20, track_abundance=track_abundance)
    E2 = Estimators(n=5, ksize=20, track_abundance=track_abundance)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.mh.add_hash(i)

    assert E1.count_common(E2) == 4
    assert E2.count_common(E1) == 4


def test_dna_mh(track_abundance):
    e1 = Estimators(n=5, ksize=4, track_abundance=track_abundance)
    e2 = Estimators(n=5, ksize=4, track_abundance=track_abundance)

    seq = 'ATGGCAGTGACGATGCCAG'
    e1.add_sequence(seq)
    for i in range(len(seq) - 3):
        e2.add(seq[i:i + 4])

    assert e1.mh.get_mins() == e2.mh.get_mins()
    print(e1.mh.get_mins())
    assert 726311917625663847 in e1.mh.get_mins()
    assert 3697418565283905118 in e1.mh.get_mins()


def test_protein_mh(track_abundance):
    e1 = Estimators(n=5, ksize=6, protein=True, track_abundance=track_abundance)
    e2 = Estimators(n=5, ksize=6, protein=True, track_abundance=track_abundance)

    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    for i in range(len(seq) - 5):
        kmer = seq[i:i + 6]
        e2.add(kmer)

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert 901193879228338100 in e1.mh.get_mins()


def test_pickle(track_abundance):
    import pickle
    from io import BytesIO

    e1 = Estimators(n=5, ksize=6, protein=False, track_abundance=track_abundance)
    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    fp = BytesIO()
    pickle.dump(e1, fp)

    fp2 = BytesIO(fp.getvalue())
    e2 = pickle.load(fp2)

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert e1.num == e2.num
    assert e1.ksize == e2.ksize
    assert e1.is_protein == e2.is_protein


def test_bad_construct_1(track_abundance):
    try:
        e1 = Estimators(ksize=6, protein=False, track_abundance=track_abundance)
        assert 0, "require n in constructor"
    except ValueError:
        pass


def test_bad_construct_2(track_abundance):
    try:
        e1 = Estimators(n=100, protein=False, track_abundance=track_abundance)
        assert 0, "require ksize in constructor"
    except ValueError:
        pass
