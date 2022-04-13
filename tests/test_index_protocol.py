"""
Tests for the 'Index' class and protocol. All Index classes should support
this functionality.
"""

import pytest

import sourmash
from sourmash import SourmashSignature
from sourmash.index import (LinearIndex, ZipFileLinearIndex,
                            LazyLinearIndex, MultiIndex,
                            StandaloneManifestIndex, LazyLoadedIndex)
from sourmash.index.sqlite_index import SqliteIndex
from sourmash.index.revindex import RevIndex
from sourmash.sbt import SBT, GraphFactory
from sourmash.manifest import CollectionManifest
from sourmash.lca.lca_db import LCA_Database, load_single_database

import sourmash_tst_utils as utils


def _load_three_sigs():
    # utility function - load & return these three sigs.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    return [ss2, ss47, ss63]


def build_linear_index(runtmp):
    ss2, ss47, ss63 = _load_three_sigs()

    lidx = LinearIndex()
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    return lidx


def build_lazy_linear_index(runtmp):
    lidx = build_linear_index(runtmp)
    return LazyLinearIndex(lidx)


def build_sbt_index(runtmp):
    ss2, ss47, ss63 = _load_three_sigs()
    
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=2)

    root.insert(ss2)
    root.insert(ss47)
    root.insert(ss63)

    return root


def build_sbt_index_save_load(runtmp):
    root = build_sbt_index(runtmp)
    out = runtmp.output('xyz.sbt.zip')
    root.save(out)

    return sourmash.load_file_as_index(out)


def build_zipfile_index(runtmp):
    from sourmash.sourmash_args import SaveSignatures_ZipFile

    location = runtmp.output('index.zip')
    with SaveSignatures_ZipFile(location) as save_sigs:
        for ss in _load_three_sigs():
            save_sigs.add(ss)

    idx = ZipFileLinearIndex.load(location)
    return idx


def build_multi_index(runtmp):
    siglist = _load_three_sigs()
    lidx = LinearIndex(siglist)

    mi = MultiIndex.load([lidx], [None], None)
    return mi


def build_standalone_manifest_index(runtmp):
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    siglist = [(ss2, sig2), (ss47, sig47), (ss63, sig63)]

    rows = []
    rows.extend((CollectionManifest.make_manifest_row(ss, loc) for ss, loc in siglist ))
    mf = CollectionManifest(rows)
    mf_filename = runtmp.output("mf.csv")
    
    mf.write_to_filename(mf_filename)

    idx = StandaloneManifestIndex.load(mf_filename)
    return idx


def build_lca_index(runtmp):
    siglist = _load_three_sigs()
    db = LCA_Database(31, 1000, 'DNA')
    for ss in siglist:
        db.insert(ss)

    return db


def build_lca_index_save_load(runtmp):
    db = build_lca_index(runtmp)
    outfile = runtmp.output('db.lca.json')
    db.save(outfile)

    return sourmash.load_file_as_index(outfile)


def build_lca_index_save_load(runtmp):
    db = build_lca_index(runtmp)
    outfile = runtmp.output('db.lca.json')
    db.save(outfile)

    return sourmash.load_file_as_index(outfile)


def build_sqlite_index(runtmp):
    filename = runtmp.output('idx.sqldb')
    db = SqliteIndex.create(filename)

    siglist = _load_three_sigs()
    for ss in siglist:
        db.insert(ss)

    return db


def build_lazy_loaded_index(runtmp):
    db = build_lca_index(runtmp)
    outfile = runtmp.output('db.lca.json')
    db.save(outfile)

    mf = CollectionManifest.create_manifest(db._signatures_with_internal())
    return LazyLoadedIndex(outfile, mf)


def build_revindex(runtmp):
    ss2, ss47, ss63 = _load_three_sigs()

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    return lidx


def build_lca_index_save_load_sql(runtmp):
    db = build_lca_index(runtmp)
    outfile = runtmp.output('db.lca.json')
    db.save(outfile, format='sql')

    x = load_single_database(outfile)
    db_load = x[0]

    return db_load


#
# create a fixture 'index_obj' that is parameterized by all of these
# building functions.
#

@pytest.fixture(params=[build_linear_index,
                        build_lazy_linear_index,
                        build_sbt_index,
                        build_zipfile_index,
                        build_multi_index,
                        build_standalone_manifest_index,
                        build_lca_index,
                        build_sbt_index_save_load,
                        build_lca_index_save_load,
                        build_sqlite_index,
                        build_lazy_loaded_index,
                        build_lca_index_save_load_sql,
#                        build_revindex,
                        ]
)
def index_obj(request, runtmp):
    build_fn = request.param

    # build on demand
    return build_fn(runtmp)


###
### generic Index tests go here
###


def test_index_search_exact_match(index_obj):
    # search for an exact match
    ss2, ss47, ss63 = _load_three_sigs()

    sr = index_obj.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0].signature.minhash == ss2.minhash
    assert sr[0].score == 1.0


def test_index_search_lower_threshold(index_obj):
    # search at a lower threshold/multiple results with ss47
    ss2, ss47, ss63 = _load_three_sigs()

    sr = index_obj.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0].signature.minhash == ss47.minhash
    assert sr[0].score == 1.0
    assert sr[1].signature.minhash == ss63.minhash
    assert round(sr[1].score, 2) == 0.32


def test_index_search_lower_threshold_2(index_obj):
    # search at a lower threshold/multiple results with ss63
    ss2, ss47, ss63 = _load_three_sigs()

    sr = index_obj.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0].signature.minhash == ss63.minhash
    assert sr[0].score == 1.0
    assert sr[1].signature.minhash == ss47.minhash
    assert round(sr[1].score, 2) == 0.32


def test_index_search_higher_threshold_2(index_obj):
    # search at a higher threshold/one match
    ss2, ss47, ss63 = _load_three_sigs()

    # search for sig63 with high threshold => 1 match
    sr = index_obj.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0].signature.minhash == ss63.minhash
    assert sr[0].score == 1.0


def test_index_search_containment(index_obj):
    # search for containment at a low threshold/multiple results with ss63
    ss2, ss47, ss63 = _load_three_sigs()

    sr = index_obj.search(ss63, do_containment=True, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0].signature.minhash == ss63.minhash
    assert sr[0].score == 1.0
    assert sr[1].signature.minhash == ss47.minhash
    assert round(sr[1].score, 2) == 0.48


def test_index_signatures(index_obj):
    # signatures works?
    siglist = list(index_obj.signatures())

    ss2, ss47, ss63 = _load_three_sigs()
    assert len(siglist) == 3

    # check md5sums, since 'in' doesn't always work
    md5s = set(( ss.md5sum() for ss in siglist ))
    assert ss2.md5sum() in md5s
    assert ss47.md5sum() in md5s
    assert ss63.md5sum() in md5s


def test_index_len(index_obj):
    # len works?
    assert len(index_obj) == 3


def test_index_bool(index_obj):
    # bool works?
    assert bool(index_obj)


def test_index_select_basic(index_obj):
    # select does the basic thing ok
    idx = index_obj.select(ksize=31, moltype='DNA', abund=False,
                           containment=True, scaled=1000, num=0, picklist=None)

    assert len(idx) == 3
    siglist = list(idx.signatures())
    assert len(siglist) == 3

    # check md5sums, since 'in' doesn't always work
    md5s = set(( ss.md5sum() for ss in siglist ))
    ss2, ss47, ss63 = _load_three_sigs()
    assert ss2.md5sum() in md5s
    assert ss47.md5sum() in md5s
    assert ss63.md5sum() in md5s


def test_index_select_nada(index_obj):
    # select works ok when nothing matches!

    # CTB: currently this EITHER raises a ValueError OR returns an empty
    # Index object, depending on implementation. :think:
    # See: https://github.com/sourmash-bio/sourmash/issues/1940
    try:
        idx = index_obj.select(ksize=21)
    except ValueError:
        idx = LinearIndex([])

    assert len(idx) == 0
    siglist = list(idx.signatures())
    assert len(siglist) == 0


def test_index_prefetch(index_obj):
    # test basic prefetch
    ss2, ss47, ss63 = _load_three_sigs()

    # search for ss2
    results = []
    for result in index_obj.prefetch(ss2, threshold_bp=0):
        results.append(result)

    assert len(results) == 1
    assert results[0].signature.minhash == ss2.minhash

    # search for ss47 - expect two results
    results = []
    for result in index_obj.prefetch(ss47, threshold_bp=0):
        results.append(result)

    assert len(results) == 2
    assert results[0].signature.minhash == ss47.minhash
    assert results[1].signature.minhash == ss63.minhash


def test_index_gather(index_obj):
    # test basic gather
    ss2, ss47, ss63 = _load_three_sigs()

    matches = index_obj.gather(ss2)
    assert len(matches) == 1
    assert matches[0].score == 1.0
    assert matches[0].signature.minhash == ss2.minhash

    matches = index_obj.gather(ss47)
    assert len(matches) == 1
    assert matches[0].score == 1.0
    assert matches[0].signature.minhash == ss47.minhash


def test_linear_gather_threshold_1(index_obj):
    # test gather() method, in some detail
    ss2, ss47, ss63 = _load_three_sigs()

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(ss2.minhash.hashes))
    new_mh = ss2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        index_obj.gather(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    results = index_obj.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        index_obj.gather(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    results = index_obj.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        index_obj.gather(SourmashSignature(new_mh), threshold_bp=5000)


def test_gather_threshold_5(index_obj):
    # test gather() method, in some detail
    ss2, ss47, ss63 = _load_three_sigs()

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(ss2.minhash.hashes.keys()))
    new_mh = ss2.minhash.copy_and_clear()

    # add five hashes
    for i in range(5):
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())

    # should get a result with no threshold (any match at all is returned)
    results = index_obj.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # now, check with a threshold_bp that should be meet-able.
    results = index_obj.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash
