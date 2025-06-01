"""
Tests of the MemRevIndex and DiskRevIndex classes.
"""

import pytest
import sourmash_tst_utils as utils
import shutil

from sourmash.index import revindex
from sourmash.index.revindex import MemRevIndex, DiskRevIndex
from sourmash.signature import load_one_signature_from_json
from sourmash.search import JaccardSearch, SearchType
from sourmash import SourmashSignature, minhash
from sourmash.picklist import SignaturePicklist

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


def test_mem_revindex_empty():
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)
    lidx = MemRevIndex(template=ss2.minhash)

    with pytest.raises(ValueError):
        list(lidx.signatures())


def test_mem_revindex_num():
    mh = minhash.MinHash(500, 4)

    with pytest.raises(ValueError):
        MemRevIndex(template=mh)


def test_mem_revindex_basic():
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)

    db = MemRevIndex(template=ss2.minhash)
    db.insert(ss2)

    db = db.select(ksize=31, scaled=1000, moltype="DNA")
    assert len(db) == 1
    assert db.location is None


def test_mem_revindex_manifest():
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)

    db = MemRevIndex(template=ss2.minhash)
    db.insert(ss2)
    db = db.select(ksize=31, scaled=1000, moltype="DNA")

    manifest = db.manifest

    assert len(db) == len(manifest)


def test_mem_revindex_index_search():
    # confirm that RevIndex works
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2
    assert sr[0][2] is None

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63
    assert sr[1][2] is None  # location

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47
    assert sr[1][2] is None  # location

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_mem_revindex_index_search_picklist(runtmp):
    # confirm that MemRevIndex works w/picklists
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    pl = SignaturePicklist("ident")
    pl.init(values=["CP001071.1"])
    lidx = lidx.select(picklist=pl)

    assert len(lidx) == 1
    assert len(list(lidx.signatures())) == 1

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2
    assert sr[0][2] is None

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    assert len(sr) == 0

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    assert len(sr) == 0


def test_mem_revindex_index_search_retrieve_orig():
    # confirm that RevIndex returns the original, not the downsampled sketches
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    # store downsampled:
    ds_mh = ss2.minhash.downsample(scaled=10_000)
    lidx = MemRevIndex(template=ds_mh)
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


def test_mem_revindex_best_containment():
    # check that RevIndex.best_containment works.
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
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


def test_mem_revindex_gather_ignore():
    # check that RevIndex gather ignores things properly.
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47, ksize=31)
    ss63 = load_one_signature_from_json(sig63, ksize=31)

    # construct an index...
    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

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


def test_mem_revindex_insert_after_init():
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx._init_inner()

    # should not work!
    with pytest.raises(Exception):
        lidx.insert(ss47)


def test_mem_revindex_union_found():
    # confirm that RevIndex works
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    # construct with 2 & 63
    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss63)

    # gather with 47
    counter = lidx.counter_gather(ss47, 0)

    # check found hashes
    ident_mh = counter.union_found
    print("XXX", ident_mh, ident_mh.scaled)
    print("YYY", ss47.minhash, ss47.minhash.scaled)
    print(ident_mh.contained_by(ss47.minhash))
    print(ss47.minhash.contained_by(ident_mh))
    assert ident_mh.contained_by(ss47.minhash) == 1.0
    assert round(ss47.minhash.contained_by(ident_mh), 5) == 0.48851


def test_mem_revindex_index_check_errors():
    # check various errors
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    empty_mh = ss2.minhash.copy_and_clear()
    empty_ss = SourmashSignature(empty_mh)

    with pytest.raises(ValueError):
        lidx.search(empty_ss)

    with pytest.raises(TypeError):
        lidx.search(ss2, threshold=None)

    with pytest.raises(TypeError):
        lidx.search(ss2, threshold=None)

    with pytest.raises(ValueError):
        lidx.counter_gather(empty_ss)

    with pytest.raises(ValueError):
        lidx.prefetch(empty_ss)


def test_mem_revindex_index_check_nomatches():
    # check what happens if no matches
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    lidx = MemRevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    pl = SignaturePicklist("ident")
    pl.init(values=["CP001071.1"])
    lidx = lidx.select(picklist=pl)

    assert len(lidx) == 1
    assert len(list(lidx.signatures())) == 1

    assert lidx.peek(ss47.minhash) == []

    with pytest.raises(ValueError):
        lidx.counter_gather(ss47)


def test_disk_revindex_basic():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    db = db.select(ksize=31, scaled=1000, moltype="DNA")
    assert len(db) == 3
    assert db.location == rocksdb_path


def test_disk_revindex_manifest():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    db = db.select(ksize=31, scaled=1000, moltype="DNA")
    assert len(db) == 3
    assert db.location == rocksdb_path

    mf = db.manifest

    assert len(db) == len(mf)
    print(list(mf.rows))


def test_disk_revindex_prefetch_to_revindex():
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    ri = db.counter_gather(ss47, threshold_bp=0)
    assert len(list(ri.signatures())) == 2


def test_disk_revindex_ksize_wrong():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    with pytest.raises(ValueError):
        db.select(ksize=21)


def test_disk_revindex_load(runtmp):
    # check loading from non .rocksdb directories
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")

    newpath = runtmp.output("foo.dir")
    shutil.copytree(rocksdb_path, newpath)

    db = DiskRevIndex(newpath)
    print(db)
    assert len(db) == 3, len(db)


def test_disk_revindex_len():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)


def test_disk_revindex_signatures():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    xx = list(db.signatures())
    assert len(xx) == 3
    for ss in xx:
        print(ss.name)
    # victory!


def test_disk_revindex_signatures_with_internal():
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


def test_disk_revindex_signatures_with_picklist():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    dataset_picks = revindex.RevIndex_DatasetPicklist([0, 1])
    db._idx_picklist = dataset_picks
    assert len(db) == 2, len(db)

    xx = list(db.signatures())
    assert len(xx) == 2
    for ss in xx:
        print(ss.name)
    # victory!


def test_disk_revindex_signatures_with_internal_with_picklist():
    # check that 'internal' matches enumeration order, and ignores
    # picklists.

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    dataset_picks = revindex.RevIndex_DatasetPicklist([0, 1])
    db._idx_picklist = dataset_picks
    assert len(db) == 2, len(db)  # len pays attention to picklist...

    # BUT: picklist is ignored by signatures_with_internal.
    xx = list(db._signatures_with_internal())
    assert len(xx) == 3
    for n, (ss, internal) in enumerate(xx):
        assert n == int(internal)
        print(ss.name)
    # victory!


def test_disk_revindex_best_containment():
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
    assert result.location == rocksdb_path


def test_disk_revindex_prefetch():
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    matches = list(db.prefetch(ss47, threshold_bp=0))
    print(matches)
    match = matches[0]
    assert match.signature.name.startswith("NC_009665.1 Shewanella baltica OS185")
    assert round(match.score, 5) == 1.0
    assert match.location == rocksdb_path

    match = matches[1]
    assert match.signature.name.startswith("NC_011663.1 Shewanella baltica OS223")
    assert round(match.score, 5) == 0.48851
    assert match.location == rocksdb_path

    assert len(matches) == 2


def test_disk_revindex_ksize_wrong():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    with pytest.raises(ValueError):
        db.select(ksize=21)


def test_disk_revindex_ksize():
    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print("xxx", db, db.select(ksize=31))
    assert db == db.select(ksize=31)


def test_create_dataset_picklist_1():
    dataset_picks = revindex.RevIndex_DatasetPicklist([0, 1])

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
    db._idx_picklist = dataset_picks

    # picklist including match:
    xx = list(db.search(ss2, threshold=0))
    assert len(xx) == 1
    assert xx[0].score == 1.0


def test_create_dataset_picklist_2():
    dataset_picks = revindex.RevIndex_DatasetPicklist([0, 1])

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
    db._idx_picklist = dataset_picks

    # picklist, 1 exact match
    xx = list(db.search(ss47, threshold=0))
    assert len(xx) == 1
    assert xx[0].score == 1.0


def test_create_dataset_picklist_3():
    dataset_picks = revindex.RevIndex_DatasetPicklist([0, 1])

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
    db._idx_picklist = dataset_picks

    # picklist, 1 inexact match
    xx = list(db.search(ss63, threshold=0))
    assert len(xx) == 1
    assert round(xx[0].score, 3) == 0.321


def test_create_dataset_picklist_4():
    # what if picklist -> empty? + prefetch?
    dataset_picks = revindex.RevIndex_DatasetPicklist([])

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
    db._idx_picklist = dataset_picks

    # picklist, 0 matches
    xx = list(db.prefetch(ss63, threshold=0))
    assert len(xx) == 0


def test_create_dataset_picklist_5():
    # make sure _signatures_with_internal reports everything.
    dataset_picks = revindex.RevIndex_DatasetPicklist([0])

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)
    print(db)
    assert len(db) == 3, len(db)

    db._idx_picklist = dataset_picks

    siglist = list(db._signatures_with_internal())
    assert len(siglist) == 3


def test_disk_revindex_against_bad_abund_rocksdb(runtmp):
    # check against a RocksDB that contains sketches w/abund,
    # created by branchwater. An alternative would be to write
    # a test that uses the internal Python API for gather.
    db = utils.get_test_data("track_abund/47+63.abund.rocksdb")
    metag = utils.get_test_data("SRR606249.sig.gz")

    runtmp.sourmash("gather", metag, db)


def test_disk_revindex_prefetch_to_cg_colors_1():
    sig63 = utils.get_test_data("63.fa.sig")
    ss63 = load_one_signature_from_json(sig63, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(ss63, threshold_bp=0)
    sr, isect_mh = cg.peek(ss63.minhash)
    assert sr.score == 1.0


def test_disk_revindex_prefetch_to_cg_colors_2():
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(ss47, threshold_bp=0)
    sr, isect_mh = cg.peek(ss47.minhash)
    assert round(sr.score, 5) == 0.48851

    cg.consume(isect_mh)
    assert cg.peek(ss47.minhash) == []


def test_disk_revindex_prefetch_to_cg_colors_3():
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(ss2, threshold_bp=0)
    sr, isect_mh = cg.peek(ss2.minhash)
    print(sr)
    assert sr.score == 1.0

    cg.consume(isect_mh)


def test_disk_revindex_prefetch_to_cg_colors_4():
    # test peek/consume on db that contains 63 and 2, but not 47
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    sig63 = utils.get_test_data("63.fa.sig")
    ss63 = load_one_signature_from_json(sig63, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    combined_mh = ss47.minhash.to_mutable().copy()
    combined_mh += ss63.minhash
    combined_ss = SourmashSignature(combined_mh)

    cg = db.counter_gather(combined_ss, threshold_bp=0)
    sr, isect_mh = cg.peek(ss47.minhash)
    assert round(sr.score, 5) == 0.48851

    cg.consume(isect_mh)

    sr, intersect_mh = cg.peek(ss63.minhash)
    print(sr)
    assert sr
    assert len(intersect_mh) == 5238


def test_disk_revindex_prefetch_to_cg_colors_5():
    # test peek/consume on db that contains 47, 63 and 2
    metag_path = utils.get_test_data("SRR606249.sig.gz")
    metag = load_one_signature_from_json(metag_path, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(metag, threshold_bp=0)

    # round 1
    sr, isect_mh = cg.peek(metag.minhash)
    assert sr.signature.name.startswith("NC_011663.1")
    print(sr.signature.name)
    assert round(sr.score, 5) == 0.01048

    assert len(list(cg.signatures())) == 3
    cg.consume(isect_mh)
    assert len(list(cg.signatures())) == 2

    # round 2
    mh = metag.minhash.to_mutable()
    mh.remove_many(isect_mh)

    sr, isect_mh = cg.peek(mh)
    assert sr.signature.name.startswith("CP001071.1")
    print(sr.signature.name)
    assert round(sr.score, 5) == 0.00529

    assert len(list(cg.signatures())) == 2
    cg.consume(isect_mh)
    assert len(list(cg.signatures())) == 1

    # round 3
    mh.remove_many(isect_mh)

    sr, isect_mh = cg.peek(mh)
    assert sr.signature.name.startswith("NC_009665.1")
    print(sr.signature.name)
    assert round(sr.score, 5) == 0.00435

    assert len(list(cg.signatures())) == 1
    cg.consume(isect_mh)
    assert len(list(cg.signatures())) == 0


def test_disk_revindex_prefetch_to_cg_colors_6():
    # test signatures on db that contains 47, 63 and 2
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(ss47, threshold_bp=0)

    siglist = list(cg.signatures())
    assert len(siglist) == 2


def test_disk_revindex_prefetch_to_cg_colors_7():
    # test signatures on db that contains 47, 63 and 2
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = load_one_signature_from_json(sig2, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(ss2, threshold_bp=0)

    siglist = list(cg.signatures())
    assert len(siglist) == 1


def test_disk_revindex_prefetch_to_cg_colors_8():
    # test peek/consume on db that contains 47, 63 and 2
    metag_path = utils.get_test_data("SRR606249.sig.gz")
    metag = load_one_signature_from_json(metag_path, ksize=31)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    cg = db.counter_gather(metag, threshold_bp=0)

    siglist = list(cg.signatures())
    assert len(siglist) == 3


def test_disk_revindex_union_found():
    # test union_found on db that contains 63 and 2, but not 47
    sig47 = utils.get_test_data("47.fa.sig")
    ss47 = load_one_signature_from_json(sig47, ksize=31)

    sig63 = utils.get_test_data("63.fa.sig")
    load_one_signature_from_json(sig63, ksize=31)

    rocksdb_path = utils.get_test_data("2sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    # gather with 47
    counter = db.counter_gather(ss47, 0)

    # check found hashes
    ident_mh = counter.union_found
    assert ident_mh.contained_by(ss47.minhash) == 1.0
    assert round(ss47.minhash.contained_by(ident_mh), 5) == 0.48851


def test_disk_revindex_index_search_picklist(runtmp):
    # confirm that disk-based RevIndex search works w/picklists
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = load_one_signature_from_json(sig2, ksize=31)
    ss47 = load_one_signature_from_json(sig47)
    ss63 = load_one_signature_from_json(sig63)

    rocksdb_path = utils.get_test_data("3sigs.branch_0913.rocksdb")
    db = DiskRevIndex(rocksdb_path)

    pl = SignaturePicklist("ident")
    pl.init(values=["CP001071.1"])
    db = db.select(picklist=pl)

    assert len(db) == 1
    assert len(list(db.signatures())) == 1

    # now, search for sig2
    sr = db.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; should be no result.
    sr = db.search(ss47, threshold=0.1)
    assert len(sr) == 0

    # search for sig63 with lower threshold; should be no result.
    sr = db.search(ss63, threshold=0.1)
    assert len(sr) == 0
