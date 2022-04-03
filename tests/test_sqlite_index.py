"Tests for SqliteIndex"

import pytest

import sourmash
from sourmash.index.sqlite_index import SqliteIndex
from sourmash import load_one_signature, SourmashSignature
from sourmash.picklist import SignaturePicklist, PickStyle

import sourmash_tst_utils as utils


def test_sqlite_index_search():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    sqlidx = SqliteIndex(":memory:")
    sqlidx.insert(ss2)
    sqlidx.insert(ss47)
    sqlidx.insert(ss63)

    # now, search for sig2
    sr = sqlidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = sqlidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = sqlidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47

    # search for sig63 with high threshold => 1 match
    sr = sqlidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_sqlite_index_prefetch():
    # prefetch does basic things right:
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    sqlidx = SqliteIndex(":memory:")
    sqlidx.insert(ss2)
    sqlidx.insert(ss47)
    sqlidx.insert(ss63)

    # search for ss2
    results = []
    for result in sqlidx.prefetch(ss2, threshold_bp=0):
        results.append(result)

    assert len(results) == 1
    assert results[0].signature == ss2

    # search for ss47 - expect two results
    results = []
    for result in sqlidx.prefetch(ss47, threshold_bp=0):
        results.append(result)

    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_sqlite_index_prefetch_empty():
    # check that an exception is raised upon for an empty database
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    sqlidx = SqliteIndex(":memory:")

    # since this is a generator, we need to actually ask for a value to
    # get exception raised.
    g = sqlidx.prefetch(ss2, threshold_bp=0)
    with pytest.raises(ValueError) as e:
        next(g)

    assert "no signatures to search" in str(e.value)


def test_sqlite_index_gather():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    sqlidx = SqliteIndex(":memory:")
    sqlidx.insert(ss2)
    sqlidx.insert(ss47)
    sqlidx.insert(ss63)

    matches = sqlidx.gather(ss2)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss2

    matches = sqlidx.gather(ss47)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47


def test_index_search_subj_scaled_is_lower():
    # check that subject sketches are appropriately downsampled
    sigfile = utils.get_test_data('scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz')
    ss = sourmash.load_one_signature(sigfile)

    # double check :)
    assert ss.minhash.scaled == 100

    # build a new query that has a scaled of 1000
    qs = SourmashSignature(ss.minhash.downsample(scaled=1000))

    # create Index to search
    sqlidx = SqliteIndex(":memory:")
    sqlidx.insert(ss)

    # search!
    results = list(sqlidx.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_sqlite_index_save_load(runtmp):
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    filename = runtmp.output('foo')
    sqlidx = SqliteIndex(filename)
    sqlidx.insert(ss2)
    sqlidx.insert(ss47)
    sqlidx.insert(ss63)

    sqlidx.close()

    sqlidx2 = SqliteIndex.load(filename)

    # now, search for sig2
    sr = sqlidx2.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2


def test_sqlite_gather_threshold_1():
    # test gather() method, in some detail
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    sqlidx = SqliteIndex(":memory:")

    sqlidx.insert(sig47)
    sqlidx.insert(sig63)
    sqlidx.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        sqlidx.gather(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    results = sqlidx.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == ":memory:"

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        sqlidx.gather(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    results = sqlidx.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == ":memory:"

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        sqlidx.gather(SourmashSignature(new_mh), threshold_bp=5000)


def test_sqlite_gather_threshold_5():
    # test gather() method above threshold
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    sqlidx = SqliteIndex(":memory:")

    sqlidx.insert(sig47)
    sqlidx.insert(sig63)
    sqlidx.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # add five hashes
    for i in range(5):
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())

    # should get a result with no threshold (any match at all is returned)
    results = sqlidx.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == ':memory:'

    # now, check with a threshold_bp that should be meet-able.
    results = sqlidx.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == ':memory:'


def test_sqlite_index_multik_select():
    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # select most specifically
    sqlidx2 = sqlidx.select(ksize=31, moltype='DNA')
    assert len(sqlidx2) == 1

    # all are DNA:
    sqlidx2 = sqlidx.select(moltype='DNA')
    assert len(sqlidx2) == 3


def test_sqlite_index_num_select():
    # this will fail on 'num' select, which is not allowed
    sqlidx = SqliteIndex(":memory:")
    with pytest.raises(ValueError):
        sqlidx.select(num=100)


def test_sqlite_index_abund_select():
    # this will fail on 'track_abundance' select, which is not allowed
    sqlidx = SqliteIndex(":memory:")
    with pytest.raises(ValueError):
        sqlidx.select(track_abundance=True)


def test_sqlite_index_moltype_select():
    # @CTB
    return

    # this loads multiple ksizes (19, 31) and moltypes (DNA, protein, hp, etc)
    filename = utils.get_test_data('prot/all.zip')
    siglist = sourmash.load_file_as_signatures(filename)

    sqlidx = SqliteIndex(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # select most specific DNA
    sqlidx2 = sqlidx.select(ksize=31, moltype='DNA')
    assert len(sqlidx2) == 2

    # select most specific protein
    sqlidx2 = sqlidx.select(ksize=19, moltype='protein')
    assert len(sqlidx2) == 2

    # can leave off ksize, selects all ksizes
    sqlidx2 = sqlidx.select(moltype='DNA')
    assert len(sqlidx2) == 2

    # can leave off ksize, selects all ksizes
    sqlidx2 = sqlidx.select(moltype='protein')
    assert len(sqlidx2) == 2

    # try hp
    sqlidx2 = sqlidx.select(moltype='hp')
    assert len(sqlidx2) == 2

    # try dayhoff
    sqlidx2 = sqlidx.select(moltype='dayhoff')
    assert len(sqlidx2) == 2

    # select something impossible
    sqlidx2 = sqlidx.select(ksize=4)
    assert len(sqlidx2) == 0


def test_sqlite_index_picklist_select():
    # test select with a picklist

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['f3a90d4e'])

    # select on picklist
    sqlidx2 = sqlidx.select(picklist=picklist)
    assert len(sqlidx2) == 1
    ss = list(sqlidx2.signatures())[0]
    assert ss.minhash.ksize == 31
    assert ss.md5sum().startswith('f3a90d4e55')


def test_sqlite_index_picklist_select_exclude():
    # test select with a picklist, but exclude

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8', pickstyle=PickStyle.EXCLUDE)
    picklist.init(['f3a90d4e'])

    # select on picklist
    sqlidx2 = sqlidx.select(picklist=picklist)
    assert len(sqlidx2) == 2
    md5s = set()
    ksizes = set()
    for ss in list(sqlidx2.signatures()):
        md5s.add(ss.md5sum())
        ksizes.add(ss.minhash.ksize)
    assert md5s == set(['f372e47893edd349e5956f8b0d8dcbf7','43f3b48e59443092850964d355a20ac0'])
    assert ksizes == set([21,51])


def test_sqlite_jaccard_ordering():
    # this tests a tricky situation where for three sketches A, B, C,
    # |A intersect B| is greater than |A intersect C|
    # _but_
    # |A jaccard B| is less than |A intersect B|
    a = sourmash.MinHash(ksize=31, n=0, scaled=2)
    b = a.copy_and_clear()
    c = a.copy_and_clear()

    a.add_many([1, 2, 3, 4])
    b.add_many([1, 2, 3] + list(range(10, 30)))
    c.add_many([1, 5])

    def _intersect(x, y):
        return x.intersection_and_union_size(y)[0]

    print('a intersect b:', _intersect(a, b))
    print('a intersect c:', _intersect(a, c))
    print('a jaccard b:', a.jaccard(b))
    print('a jaccard c:', a.jaccard(c))
    assert _intersect(a, b) > _intersect(a, c)
    assert a.jaccard(b) < a.jaccard(c)

    # thresholds to use:
    assert a.jaccard(b) < 0.15
    assert a.jaccard(c) > 0.15

    # now - make signatures, try out :)
    ss_a = sourmash.SourmashSignature(a, name='A')
    ss_b = sourmash.SourmashSignature(b, name='B')
    ss_c = sourmash.SourmashSignature(c, name='C')

    sqlidx = SqliteIndex(":memory:")
    sqlidx.insert(ss_a)
    sqlidx.insert(ss_b)
    sqlidx.insert(ss_c)

    sr = sqlidx.search(ss_a, threshold=0.15)
    print(sr)
    assert len(sr) == 2
    assert sr[0].signature == ss_a
    assert sr[1].signature == ss_c
