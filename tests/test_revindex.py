"""
Tests of the RevIndex and DiskRevIndex classes.
"""
import pytest
import sourmash_tst_utils as utils
import shutil

from sourmash.index import revindex
from sourmash.index.revindex import RevIndex, DiskRevIndex
from sourmash.signature import load_one_signature_from_json
from sourmash.search import JaccardSearch, SearchType

##
## test a slightly outre version of JaccardSearch - this is a test of the
## JaccardSearch 'collect' protocol, in particular...
##


class JaccardSearchBestOnly_ButIgnore(JaccardSearch):
    "A class that ignores certain results, but still does all the pruning."

    def __init__(self, ignore_list):
        super().__init__(SearchType.JACCARD, threshold=0.1)
        self.ignore_list = ignore_list

    # a collect function that _ignores_ things in the ignore_list
    def collect(self, score, match):
        print("in collect; current threshold:", self.threshold)
        for q in self.ignore_list:
            print("ZZZ", match, match.similarity(q))
            if match.similarity(q) == 1.0:
                print("yes, found.")
                return False

        # update threshold if not perfect match, which could help prune.
        self.threshold = score
        return True


def test_revindex_index_search():
    # confirm that RevIndex works
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_revindex_best_containment():
    # check that RevIndex.best_containment works.
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    match = lidx.best_containment(ss2)
    assert match
    assert match.score == 1.0
    assert match.signature == ss2

    match = lidx.best_containment(ss47)
    assert match
    assert match.score == 1.0
    assert match.signature == ss47


def test_revindex_gather_ignore():
    # check that RevIndex gather ignores things properly.
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47, ksize=31)
    ss63 = load_one_signature_from_json(sig63, ksize=31)

    # construct an index...
    lidx = RevIndex(template=ss2.minhash, signatures=[ss2, ss47, ss63])

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(lidx.find(search_fn, ss47))
    results = [ss.signature for ss in results]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_rocksdb_load(runtmp):
    # check loading from non .rocksdb directories
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")

    newpath = runtmp.output('foo.dir')
    shutil.copytree(rocksdb_path, newpath)
    
    db = DiskRevIndex(newpath)
    print(db)
    assert len(db) == 3, len(db)


def test_rocksdb_len():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)


def test_rocksdb_signatures():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    xx = list(db.signatures())
    assert len(xx) == 3
    for ss in xx:
        print(ss.name)
    # victory!


def test_rocksdb_signatures_with_internal():
    # check that 'internal' matches enumeration order.
    # CTB note: this is, for now, an important implementation detail,
    # implemented by the Python layer but corresponding to the Rust
    # behavior. Be careful about changing it :). The better thing
    # to do would be to export the manifest directly from Rust...
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    xx = list(db._signatures_with_internal())
    assert len(xx) == 3
    for n, (ss, internal) in enumerate(xx):
        assert n == int(internal)
        print(ss.name)
    # victory!


def test_rocksdb_best_containment():
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    result = db.best_containment(ss47)
    print(result)
    assert round(result.score, 5) == 0.48851, result
    assert (
        result.signature.name == "NC_011663.1 Shewanella baltica OS223, complete genome"
    ), result.signature.name


def test_rocksdb_prefetch():
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    matches = list(db.prefetch(ss47, threshold_bp=0))
    print(matches)
    match = matches[0]
    assert match.signature.name.startswith("NC_009665.1 Shewanella baltica OS185")
    assert round(match.score, 5) == 1.0

    match = matches[1]
    assert match.signature.name.startswith("NC_011663.1 Shewanella baltica OS223")
    assert round(match.score, 5) == 0.48851

    assert len(matches) == 2


def test_rocksdb_ksize_wrong():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    with pytest.raises(ValueError):
        db.select(ksize=21)


def test_rocksdb_ksize():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print("xxx", db, db.select(ksize=31))
    assert db == db.select(ksize=31)


def test_create_dataset_picklist_1():
    dataset_picks = revindex.DiskRevIndex_DatasetPicklist([0, 1])

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)

    # no picklist
    xx = list(db.search(ss2, threshold=0))
    assert len(xx) == 1

    # forcibly set picklist for now
    db.idx_picklist = dataset_picks

    # picklist including match:
    xx = list(db.search(ss2, threshold=0))
    assert len(xx) == 1
    assert xx[0].score == 1.0


def test_create_dataset_picklist_2():
    dataset_picks = revindex.DiskRevIndex_DatasetPicklist([0, 1])

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    # no picklist, 2 matches
    xx = list(db.search(ss47, threshold=0))
    assert len(xx) == 2

    # forcibly set picklist for now
    db.idx_picklist = dataset_picks

    # picklist, 1 exact match
    xx = list(db.search(ss47, threshold=0))
    assert len(xx) == 1
    assert xx[0].score == 1.0


def test_create_dataset_picklist_3():
    dataset_picks = revindex.DiskRevIndex_DatasetPicklist([0, 1])

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    sig63 = utils.get_test_data("63.fa.sig")
    ss63 = load_one_signature_from_json(sig63, ksize=31)

    # no picklist
    xx = list(db.search(ss63, threshold=0))
    assert len(xx) == 2

    # forcibly set picklist for now
    db.idx_picklist = dataset_picks

    # picklist, 1 inexact match
    xx = list(db.search(ss63, threshold=0))
    assert len(xx) == 1
    assert round(xx[0].score, 3) == 0.321
