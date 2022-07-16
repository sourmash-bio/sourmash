"""
Tests for the 'Index' class and protocol. All Index classes should support
this functionality.
"""

import pytest
import glob

import sourmash
from sourmash import SourmashSignature
from sourmash.index import (LinearIndex, ZipFileLinearIndex,
                            LazyLinearIndex, MultiIndex,
                            StandaloneManifestIndex,
                            IndexSearchResult)
from sourmash.index import CounterGather
from sourmash.index.sqlite_index import SqliteIndex
from sourmash.index.revindex import RevIndex
from sourmash.sbt import SBT, GraphFactory
from sourmash.manifest import CollectionManifest, BaseCollectionManifest
from sourmash.lca.lca_db import LCA_Database, load_single_database
from sourmash.minhash import (flatten_and_intersect_scaled,
                              flatten_and_downsample_scaled)

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


def build_sqlite_index(runtmp):
    filename = runtmp.output('idx.sqldb')
    db = SqliteIndex.create(filename)

    siglist = _load_three_sigs()
    for ss in siglist:
        db.insert(ss)

    return db


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


def test_index_signatures_with_location(index_obj):
    # signatures_with_location works?
    siglist = list(index_obj.signatures_with_location())

    ss2, ss47, ss63 = _load_three_sigs()
    assert len(siglist) == 3

    # check md5sums, since 'in' doesn't always work
    md5s = set(( ss.md5sum() for ss, loc in siglist ))
    assert ss2.md5sum() in md5s
    assert ss47.md5sum() in md5s
    assert ss63.md5sum() in md5s


def test_index_len(index_obj):
    # len works?
    assert len(index_obj) == 3


def test_index_bool(index_obj):
    # bool works?
    assert bool(index_obj)


def test_index_location(index_obj):
    # location works?
    assert str(index_obj.location)


def test_index_manifest(index_obj):
    # manifest is either None or a BaseCollectionManifest
    manifest = index_obj.manifest
    if manifest is not None:
        assert isinstance(manifest, BaseCollectionManifest)


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


def test_index_best_containment(index_obj):
    # test basic containment search
    ss2, ss47, ss63 = _load_three_sigs()

    match = index_obj.best_containment(ss2)
    assert match
    assert match.score == 1.0
    assert match.signature.minhash == ss2.minhash

    match = index_obj.best_containment(ss47)
    assert match
    assert match.score == 1.0
    assert match.signature.minhash == ss47.minhash


def test_index_best_containment_threshold_1(index_obj):
    # test best_containment() method, in some detail
    ss2, ss47, ss63 = _load_three_sigs()

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(ss2.minhash.hashes))
    new_mh = ss2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        index_obj.best_containment(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    result = index_obj.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        index_obj.best_containment(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    result = index_obj.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        index_obj.best_containment(SourmashSignature(new_mh), threshold_bp=5000)


def test_best_containment_threshold_5(index_obj):
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
    result = index_obj.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash

    # now, check with a threshold_bp that should be meet-able.
    result = index_obj.best_containment(SourmashSignature(new_mh), threshold_bp=5000)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == ss2.minhash


###
### CounterGather tests
###


def create_basic_counter_gather(runtmp):
    "Construct a CounterGather class."
    return CounterGather


class CounterGather_LinearIndex:
    """
    Provides an (inefficient) CounterGather-style class, for
    protocol testing purposes.
    """
    def __init__(self, orig_query):
        "Constructor - take a SourmashSignature that is the original query."
        orig_query_mh = orig_query.minhash
        if orig_query_mh.scaled == 0:
            raise ValueError

        # Index object used to actually track matches.
        self.idx = LinearIndex()
        self.orig_query_mh = orig_query_mh.copy().flatten()
        self.query_started = 0
        self.scaled = orig_query_mh.scaled
        self.locations = {}

    def add(self, ss, *, location=None, require_overlap=True):
        "Insert potential match."
        if self.query_started:
            raise ValueError("cannot add more signatures to counter after peek/consume")

        # skip duplicates
        md5 = ss.md5sum()
        if md5 in self.locations:
            return

        # confirm that this match has an overlap...
        add_mh = ss.minhash.flatten()
        overlap = self.orig_query_mh.count_common(add_mh, downsample=True)

        # ...figure out what scaled we are operating at now...
        if overlap:
            self.downsample(add_mh.scaled)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

        # ...and add to the Index, while also tracking location!
        self.idx.insert(ss)
        self.locations[md5] = location

    def signatures(self):
        "Yield all signatures"
        return self.idx.signatures()

    def downsample(self, scaled):
        "Track highest scaled across all possible matches."
        if scaled > self.scaled:
            self.scaled = scaled
        return self.scaled

    def peek(self, cur_query_mh, *, threshold_bp=0):
        """
        Find best match to current query within this CounterGather object.
        """
        self.query_started = 1

        scaled = self.downsample(cur_query_mh.scaled)
        cur_query_mh = flatten_and_downsample_scaled(cur_query_mh, scaled)

        # no hashes remaining? exit.
        if not self.orig_query_mh or not cur_query_mh:
            return []

        # verify current query is a subset of the original.
        if cur_query_mh.contained_by(self.orig_query_mh, downsample=True) < 1:
            raise ValueError("current query not a subset of original query")

        # did we get a match?
        res = self.idx.peek(cur_query_mh, threshold_bp=threshold_bp)
        if not res:
            return []
        sr, intersect_mh = res

        # got match - replace location & return.
        match = sr.signature
        md5 = match.md5sum()
        location = self.locations[md5]
        return IndexSearchResult(sr.score, match, location), intersect_mh

    def consume(self, *args, **kwargs):
        self.query_started = 1
        return self.idx.consume(*args, **kwargs)


class CounterGather_LCA:
    """
    Provides an alternative implementation of a CounterGather-style class,
    based on LCA_Database. This is currently just for protocol
    and API testing purposes.
    """
    def __init__(self, query):
        from sourmash.lca.lca_db import LCA_Database

        query_mh = query.minhash
        if query_mh.scaled == 0:
            raise ValueError("must use scaled MinHash")

        self.orig_query_mh = query_mh
        lca_db = LCA_Database(query_mh.ksize, query_mh.scaled,
                              query_mh.moltype)
        self.db = lca_db
        self.siglist = {}
        self.locations = {}
        self.query_started = 0

    def add(self, ss, *, location=None, require_overlap=True):
        "Add this signature into the counter."
        if self.query_started:
            raise ValueError("cannot add more signatures to counter after peek/consume")

        overlap = self.orig_query_mh.count_common(ss.minhash, True)
        if overlap:
            self.downsample(ss.minhash.scaled)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

        self.db.insert(ss)

        md5 = ss.md5sum()
        self.siglist[md5] = ss
        self.locations[md5] = location

    def signatures(self):
        "Yield all signatures."
        for ss in self.siglist.values():
            yield ss

    def downsample(self, scaled):
        "Track highest scaled across all possible matches."
        if scaled > self.db.scaled:
            self.db.downsample_scaled(scaled)
        return self.db.scaled

    def peek(self, query_mh, *, threshold_bp=0):
        "Return next possible match."
        from sourmash import SourmashSignature

        self.query_started = 1
        scaled = self.downsample(query_mh.scaled)
        query_mh = query_mh.downsample(scaled=scaled)

        if not self.orig_query_mh or not query_mh:
            return []

        if query_mh.contained_by(self.orig_query_mh, downsample=True) < 1:
            raise ValueError("current query not a subset of original query")

        query_ss = SourmashSignature(query_mh)

        # returns search_result, intersect_mh
        try:
            result = self.db.best_containment(query_ss, threshold_bp=threshold_bp)
        except ValueError:
            result = None

        if not result:
            return []

        cont = result.score
        match = result.signature

        intersect_mh = flatten_and_intersect_scaled(result.signature.minhash,
                                                    query_mh)

        md5 = result.signature.md5sum()
        location = self.locations[md5]

        new_sr = IndexSearchResult(cont, match, location)
        return [new_sr, intersect_mh]

    def consume(self, intersect_mh):
        self.query_started = 1


@pytest.fixture(params=[CounterGather,
                        CounterGather_LinearIndex,
                        CounterGather_LCA,
                        ]
)
def counter_gather_constructor(request):
    build_fn = request.param

    # build on demand
    return build_fn


def test_counter_get_signatures(counter_gather_constructor):
    # test .signatures() method
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(10, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(15, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    siglist = list(counter.signatures())
    assert len(siglist) == 3
    assert match_ss_1 in siglist
    assert match_ss_2 in siglist
    assert match_ss_3 in siglist


def _consume_all(query_mh, counter, threshold_bp=0):
    results = []
    query_mh = query_mh.to_mutable()

    last_intersect_size = None
    while 1:
        result = counter.peek(query_mh, threshold_bp=threshold_bp)
        if not result:
            break

        sr, intersect_mh = result
        print(sr.signature.name, len(intersect_mh))
        if last_intersect_size:
            assert len(intersect_mh) <= last_intersect_size

        last_intersect_size = len(intersect_mh)

        counter.consume(intersect_mh)
        query_mh.remove_many(intersect_mh.hashes)

        results.append((sr, len(intersect_mh)))

    return results


def test_counter_gather_1(counter_gather_constructor):
    # check a contrived set of non-overlapping gather results,
    # generated via CounterGather
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(10, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(15, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_b(counter_gather_constructor):
    # check a contrived set of somewhat-overlapping gather results,
    # generated via CounterGather. Here the overlaps are structured
    # so that the gather results are the same as those in
    # test_counter_gather_1(), even though the overlaps themselves are
    # larger.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_c_with_threshold(counter_gather_constructor):
    # check a contrived set of somewhat-overlapping gather results,
    # generated via CounterGather. Here the overlaps are structured
    # so that the gather results are the same as those in
    # test_counter_gather_1(), even though the overlaps themselves are
    # larger.
    # use a threshold, here.

    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter,
                           threshold_bp=3)

    expected = (['match1', 10],
                ['match2', 5])
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_d_diff_scaled(counter_gather_constructor):
    # test as above, but with different scaled.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear().downsample(scaled=10)
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().downsample(scaled=20)
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_d_diff_scaled_query(counter_gather_constructor):
    # test as above, but with different scaled for QUERY.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))

    match_mh_1 = query_mh.copy_and_clear().downsample(scaled=10)
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().downsample(scaled=20)
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # downsample query now -
    query_ss = SourmashSignature(query_mh.downsample(scaled=100), name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_e_abund_query(counter_gather_constructor):
    # test as above, but abund query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1, track_abundance=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear().flatten()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().flatten()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().flatten()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    # must flatten before peek!
    results = _consume_all(query_ss.minhash.flatten(), counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_f_abund_match(counter_gather_constructor):
    # test as above, but abund query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1, track_abundance=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh.flatten(), name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    # must flatten before peek!
    results = _consume_all(query_ss.minhash.flatten(), counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_2(counter_gather_constructor):
    # check basic set of gather results on semi-real data,
    # generated via CounterGather
    testdata_combined = utils.get_test_data('gather/combined.sig')
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_ss = sourmash.load_one_signature(testdata_combined, ksize=21)
    subject_sigs = [ (sourmash.load_one_signature(t, ksize=21), t)
                     for t in testdata_sigs ]

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    for ss, loc in subject_sigs:
        counter.add(ss, location=loc)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['NC_003198.1', 487],
                ['NC_000853.1', 192],
                ['NC_011978.1', 169],
                ['NC_002163.1', 157],
                ['NC_003197.2', 152],
                ['NC_009486.1', 92],
                ['NC_006905.1', 76],
                ['NC_011080.1', 59],
                ['NC_011274.1', 42],
                ['NC_006511.1', 31],
                ['NC_011294.1', 7],
                ['NC_004631.1', 2])
    assert len(results) == len(expected)

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]
        print(sr_name, size)

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_exact_match(counter_gather_constructor):
    # query == match
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter; provide a location override, too.
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    results = _consume_all(query_ss.minhash, counter)
    assert len(results) == 1
    (sr, intersect_mh) = results[0]

    assert sr.score == 1.0
    assert sr.signature == query_ss
    assert sr.location == 'somewhere over the rainbow'


def test_counter_gather_multiple_identical_matches(counter_gather_constructor):
    # test multiple identical matches being inserted, with only one return
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # create counter...
    counter = counter_gather_constructor(query_ss)

    # now add multiple identical matches.
    match_mh = query_mh.copy_and_clear()
    match_mh.add_many(range(5, 15))

    for name in 'match1', 'match2', 'match3':
        match_ss = SourmashSignature(match_mh, name=name)
        counter.add(match_ss, location=name)

    results = _consume_all(query_ss.minhash, counter)
    assert len(results) == 1

    sr, overlap_count = results[0]
    assert sr.score == 0.5
    assert overlap_count == 10

    # any one of the three is valid
    assert sr.location in ('match1', 'match2', 'match3')


def test_counter_gather_add_after_peek(counter_gather_constructor):
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    counter.peek(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, location="try again")


def test_counter_gather_add_after_consume(counter_gather_constructor):
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    counter.consume(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, location="try again")


def test_counter_gather_consume_empty_intersect(counter_gather_constructor):
    # check that consume works fine when there is an empty signature.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    # nothing really happens here :laugh:, just making sure there's no error
    counter.consume(query_ss.minhash.copy_and_clear())


def test_counter_gather_empty_initial_query(counter_gather_constructor):
    # check empty initial query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1, require_overlap=False)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_num_query(counter_gather_constructor):
    # check num query
    query_mh = sourmash.MinHash(n=500, ksize=31)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    with pytest.raises(ValueError):
        counter_gather_constructor(query_ss)


def test_counter_gather_empty_cur_query(counter_gather_constructor):
    # test empty cur query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    results = _consume_all(cur_query_mh, counter)
    assert results == []


def test_counter_gather_add_num_matchy(counter_gather_constructor):
    # test add num query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh = sourmash.MinHash(n=500, ksize=31)
    match_mh.add_many(range(0, 20))
    match_ss = SourmashSignature(match_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    with pytest.raises(ValueError):
        counter.add(match_ss, location='somewhere over the rainbow')


def test_counter_gather_bad_cur_query(counter_gather_constructor):
    # test cur query that is not subset of original query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(query_ss, location='somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    cur_query_mh.add_many(range(20, 30))
    with pytest.raises(ValueError):
        counter.peek(cur_query_mh)


def test_counter_gather_add_no_overlap(counter_gather_constructor):
    # check adding match with no overlap w/query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(10, 20))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    with pytest.raises(ValueError):
        counter.add(match_ss_1)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_big_threshold(counter_gather_constructor):
    # check 'peek' with a huge threshold
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = counter_gather_constructor(query_ss)
    counter.add(match_ss_1)

    # impossible threshold:
    threshold_bp=30*query_ss.minhash.scaled
    results = counter.peek(query_ss.minhash, threshold_bp=threshold_bp)
    assert results == []


def test_counter_gather_empty_counter(counter_gather_constructor):
    # check empty counter
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    # empty counter!
    counter = counter_gather_constructor(query_ss)

    assert counter.peek(query_ss.minhash) == []
