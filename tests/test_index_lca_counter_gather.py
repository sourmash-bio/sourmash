import pytest
import glob

import sourmash
from sourmash import SourmashSignature
from sourmash.index import CounterGather
from sourmash.lca import LCA_Database
import sourmash_tst_utils as utils


class CounterGather_LCA:
    def __init__(self, mh):
        if mh.scaled == 0:
            raise ValueError("must use scaled MinHash")

        self.orig_query_mh = mh
        lca_db = LCA_Database(mh.ksize, mh.scaled, mh.moltype)
        self.db = lca_db
        self.siglist = []
        self.locations = []
        self.query_started = 0

    def add(self, ss, location=None, require_overlap=True):
        if self.query_started:
            raise ValueError("cannot add more signatures to counter after peek/consume")

        overlap = self.orig_query_mh.count_common(ss.minhash, True)
        if not overlap and require_overlap:
            raise ValueError("no overlap between query and signature!?")

        self.db.insert(ss)
        self.siglist.append(ss)
        self.locations.append(location)

    def peek(self, query_mh, threshold_bp=0):
        self.query_started = 1
        if not self.orig_query_mh or not query_mh:
            return []

        if query_mh.contained_by(self.orig_query_mh, downsample=True) < 1:
            raise ValueError("current query not a subset of original query")

        query_ss = SourmashSignature(query_mh)

        # returns search_result, intersect_mh
        try:
            result = self.db.gather(query_ss, threshold_bp=threshold_bp)
        except ValueError:
            result = None

        if not result:
            return []

        sr = result[0]
        match_mh = sr.signature.minhash
        scaled = max(query_mh.scaled, match_mh.scaled)
        match_mh = match_mh.downsample(scaled=scaled).flatten()
        query_mh = query_mh.downsample(scaled=scaled)
        intersect_mh = match_mh & query_mh

        return [sr, intersect_mh]

    def consume(self, intersect_mh):
        self.query_started = 1


def _consume_all(query_mh, counter, threshold_bp=0):
    results = []
    query_mh = query_mh.to_mutable()

    last_intersect_size = None
    while 1:
        result = counter.peek(query_mh, threshold_bp)
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


def test_counter_gather_1():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_b():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_c_with_threshold():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_d_diff_scaled():
    # test as above, but with different scaled.
    # @CTB munged the scaleds... won't work with LCA database.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=30)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_d_diff_scaled_query():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_e_abund_query():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_1_f_abund_match():
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
    counter = CounterGather_LCA(query_ss.minhash)
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


def test_counter_gather_2():
    # check basic set of gather results on semi-real data,
    # generated via CounterGather
    testdata_combined = utils.get_test_data('gather/combined.sig')
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_ss = sourmash.load_one_signature(testdata_combined, ksize=21)
    subject_sigs = [ (sourmash.load_one_signature(t, ksize=21), t)
                     for t in testdata_sigs ]

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    for ss, loc in subject_sigs:
        counter.add(ss, loc)

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


def test_counter_gather_exact_match():
    # query == match
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    results = _consume_all(query_ss.minhash, counter)
    assert len(results) == 1
    (sr, intersect_mh) = results[0]

    assert sr.score == 1.0
    assert sr.signature == query_ss
    # assert sr.location == 'somewhere over the rainbow' @CTB


def test_counter_gather_add_after_peek():
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    counter.peek(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, "try again")


def test_counter_gather_add_after_consume():
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    counter.consume(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, "try again")


def test_counter_gather_consume_empty_intersect():
    # check that consume works fine when there is an empty signature.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    # nothing really happens here :laugh:, just making sure there's no error
    counter.consume(query_ss.minhash.copy_and_clear())


def test_counter_gather_empty_initial_query():
    # check empty initial query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(match_ss_1, require_overlap=False)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_num_query():
    # check num query
    query_mh = sourmash.MinHash(n=500, ksize=31)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    with pytest.raises(ValueError):
        counter = CounterGather_LCA(query_ss.minhash)


def test_counter_gather_empty_cur_query():
    # test empty cur query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    results = _consume_all(cur_query_mh, counter)
    assert results == []


def test_counter_gather_add_num_matchy():
    # test add num query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh = sourmash.MinHash(n=500, ksize=31)
    match_mh.add_many(range(0, 20))
    match_ss = SourmashSignature(match_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    with pytest.raises(ValueError):
        counter.add(match_ss, 'somewhere over the rainbow')


def test_counter_gather_bad_cur_query():
    # test cur query that is not subset of original query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    cur_query_mh.add_many(range(20, 30))
    with pytest.raises(ValueError):
        counter.peek(cur_query_mh)


def test_counter_gather_add_no_overlap():
    # check adding match with no overlap w/query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(10, 20))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    with pytest.raises(ValueError):
        counter.add(match_ss_1)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_big_threshold():
    # check 'peek' with a huge threshold
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather_LCA(query_ss.minhash)
    counter.add(match_ss_1)

    # impossible threshold:
    threshold_bp=30*query_ss.minhash.scaled
    results = counter.peek(query_ss.minhash, threshold_bp=threshold_bp)
    assert results == []


def test_counter_gather_empty_counter():
    # check empty counter
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    # empty counter!
    counter = CounterGather_LCA(query_ss.minhash)

    assert counter.peek(query_ss.minhash) == []



