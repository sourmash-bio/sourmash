"Tests for search.py code."

import pytest
import numpy as np
import sourmash_tst_utils as utils

from sourmash import search, SourmashSignature, MinHash, load_one_signature
from sourmash.search import (make_jaccard_search_query,
                             make_containment_query,
                             SearchResult, PrefetchResult, GatherResult)
from sourmash.index import LinearIndex


def test_make_jaccard_search_query():
    search_obj = make_jaccard_search_query(threshold=0)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_cont():
    search_obj = make_jaccard_search_query(do_containment=True,
                                           threshold=0)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_max_cont():
    search_obj = make_jaccard_search_query(do_max_containment=True,
                                           threshold=0)

    assert search_obj.score_fn == search_obj.score_max_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_best_only():
    search_obj = make_jaccard_search_query(best_only=True)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert type(search_obj) == search.JaccardSearchBestOnly


def test_make_jaccard_search_query_no_threshold_none():
    search_obj = make_jaccard_search_query(threshold=None)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_cont_and_max_cont():
    with pytest.raises(TypeError) as exc:
        search_obj = make_jaccard_search_query(do_containment=True,
                                               do_max_containment=True)

    assert str(exc.value) == "'do_containment' and 'do_max_containment' cannot both be True"


def test_cont_requires_scaled():
    search_obj = make_jaccard_search_query(do_containment=True)
    assert search_obj.require_scaled
    
    mh = MinHash(n=500, ksize=31)
    with pytest.raises(TypeError) as exc:
        search_obj.check_is_compatible(SourmashSignature(mh))
    assert str(exc.value) == "this search requires a scaled signature"


def test_search_requires_flat():
    search_obj = make_jaccard_search_query()
    
    mh = MinHash(n=500, ksize=31, track_abundance=True)
    with pytest.raises(TypeError) as exc:
        search_obj.check_is_compatible(SourmashSignature(mh))
    assert str(exc.value) == "this search cannot be done with an abund signature"


def test_score_jaccard_similarity():
    search_obj = make_jaccard_search_query()

    assert search_obj.score_fn(None, 100, None, 200) == 0.5


def test_score_jaccard_containment():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(100, 50, None, 0) == 0.5


def test_score_jaccard_containment_zero_query_size():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(0, 100, None, None) == 0


def test_score_jaccard_max_containment_1():
    search_obj = make_jaccard_search_query(do_max_containment=True)

    assert search_obj.score_fn(150, 75, 100, None) == 0.75


def test_score_jaccard_max_containment_2():
    search_obj = make_jaccard_search_query(do_max_containment=True)

    assert search_obj.score_fn(100, 75, 150, None) == 0.75


def test_score_jaccard_max_containment_zero_query_size():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(0, 100, None, None) == 0


def test_collect():
    search_obj = make_jaccard_search_query(threshold=0)
    search_obj.collect(1.0, None)
    assert search_obj.threshold == 0


def test_collect_best_only():
    search_obj = make_jaccard_search_query(threshold=0, best_only=True)
    search_obj.collect(1.0, None)
    assert search_obj.threshold == 1.0


def test_make_containment_query():
    # test basic make_containment_query call
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    search_obj = make_containment_query(mh, 5e4)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0.5


def test_make_containment_query_no_threshold():
    # test basic make_containment_query call
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    search_obj = make_containment_query(mh, None)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_containment_query_num_minhash():
    # will fail on non-scaled minhash
    mh = MinHash(n=500, ksize=31)

    for i in range(100):
        mh.add_hash(i)

    with pytest.raises(TypeError) as exc:
        search_obj = make_containment_query(mh, 5e4)

    assert str(exc.value) == "query signature must be calculated with scaled"


def test_make_containment_query_empty_minhash():
    # will fail on non-scaled minhash
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    with pytest.raises(TypeError) as exc:
        search_obj = make_containment_query(mh, -1)

    assert str(exc.value) == "threshold_bp must be non-negative"


def test_make_containment_query_high_threshold():
    # will fail on non-scaled minhash
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    # effective threshold > 1; raise ValueError
    with pytest.raises(ValueError):
        search_obj = make_containment_query(mh, 200000)


class FakeIndex(LinearIndex):
    _signatures = []
    filename = "something_or_other"

    def __init__(self, validator_fn):
        self.validator = validator_fn

    def find(self, search_fn, query, *args, **kwargs):
        if self.validator:
            self.validator(search_fn, query, args, kwargs)
        else:
            assert 0, "what are we even doing here?"
        return []


def test_index_search_passthru():
    # check that kwargs are passed through from 'search' to 'find'
    query = None

    def validate_kwarg_passthru(search_fn, query, args, kwargs):
        assert "this_kw_arg" in kwargs
        assert kwargs["this_kw_arg"] == 5

    idx = FakeIndex(validate_kwarg_passthru)

    idx.search(query, threshold=0.0, this_kw_arg=5)


def test_index_containment_passthru():
    # check that kwargs are passed through from 'search' to 'find'
    query = None

    def validate_kwarg_passthru(search_fn, query, args, kwargs):
        assert "this_kw_arg" in kwargs
        assert kwargs["this_kw_arg"] == 5

    idx = FakeIndex(validate_kwarg_passthru)

    idx.search(query, threshold=0.0, this_kw_arg=5)


def test_search_with_abund_query():
    mh = MinHash(n=0, ksize=31, scaled=1, track_abundance=True)
    query = SourmashSignature(mh)

    with pytest.raises(TypeError):
        search.search_databases_with_abund_query(query, [],
                                                 threshold=0,
                                                 do_containment=True)

    with pytest.raises(TypeError):
        search.search_databases_with_abund_query(query, [],
                                                 threshold=0,
                                                 do_max_containment=True)


def test_scaledSearchResult():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    ss4763 = ss4763.to_mutable()
    ss4763.filename = ss4763_file

    scaled = ss47.minhash.scaled

    res = SearchResult(ss47, ss4763, cmp_scaled=scaled, similarity= ss47.contained_by(ss4763))

    assert res.query_name == ss47.name
    assert res.match_name == ss4763.name
    assert res.query_scaled == ss47.minhash.scaled == 1000
    assert res.match_scaled == ss4763.minhash.scaled == 1000
    assert res.cmp_scaled == 1000
    assert res.query_abundance == ss47.minhash.track_abundance
    assert res.match_abundance == ss4763.minhash.track_abundance
#    assert res.query_bp == len(ss47.minhash) * scaled
#    assert res.match_bp == len(ss4763.minhash) * scaled
    assert res.ksize == 31
    assert res.moltype == 'DNA'
    assert res.query_filename == '47.fa'
    assert res.match_filename == ss4763_file
    assert res.query_md5 == ss47.md5sum()
    assert res.match_md5 == ss4763.md5sum()
 #   assert res.query_n_hashes == len(ss47.minhash)
 #   assert res.match_n_hashes == len(ss4763.minhash)
    assert res.md5 == ss4763.md5sum()
    assert res.name == ss4763.name
    assert res.filename == ss4763.filename
    queryc_ani = ss47.containment_ani(ss4763)
    matchc_ani = ss4763.containment_ani(ss47)
    # check that we _can_ get avg_containment_ani
    assert res.cmp.avg_containment_ani == np.mean([queryc_ani.ani, matchc_ani.ani])

def test_numSearchResult():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('num/47.fa.sig')
    ss63_file = utils.get_test_data('num/63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss63 = load_one_signature(ss63_file, ksize=31, select_moltype='dna')
    ss63 = ss63.to_mutable()
    ss63.filename = ss63_file

    assert ss47.minhash.num and ss63.minhash.num

    res = SearchResult(ss47, ss63, similarity= ss47.jaccard(ss63))
    print(res.cmp_num)
    assert res.mh1.num
    assert res.cmp.cmp_num == 500
    assert res.query_name == ss47.name
    assert res.match_name == ss63.name
    assert res.query_num == ss47.minhash.num == 500
    assert res.match_num == ss63.minhash.num == 500
    assert res.query_abundance == ss47.minhash.track_abundance
    assert res.match_abundance == ss63.minhash.track_abundance
    assert res.ksize == 31
    assert res.moltype == 'DNA'
    assert res.query_filename == '47.fa'
    assert res.match_filename == ss63_file
    assert res.query_md5 == ss47.md5sum()
    assert res.match_md5 == ss63.md5sum()
    assert res.md5 == ss63.md5sum()
    assert res.name == ss63.name
    assert res.filename == ss63.filename

    # check that we can't get ani
    with pytest.raises(TypeError) as exc:
        res.estimate_search_ani()
    assert("ANI can only be estimated from scaled signatures.") in str(exc)

    # get result as dictionary (of just items to write)
    resD = res.resultdict
    assert resD["filename"] == res.filename
    assert resD["name"] == res.name
    assert resD["similarity"] == res.similarity


def test_SearchResult_incompatible_sigs():
    ss47_file = utils.get_test_data('num/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    with pytest.raises(TypeError) as exc:
        SearchResult(ss47, ss4763, similarity=10)
    print(str(exc))
    assert "Error: Both sketches must be 'num' or 'scaled'." in str(exc)


def test_SearchResult_notsigs():
    ss47_file = utils.get_test_data('num/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')

    with pytest.raises(AttributeError) as exc:
        SearchResult(ss47_file, ss4763_file, similarity=10)
    print(str(exc))
    assert "'str' object has no attribute 'minhash'" in str(exc)


def test_SearchResult_no_similarity():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    with pytest.raises(ValueError) as exc:
        SearchResult(ss47, ss4763)
    print(str(exc))
    assert "Error: Must provide 'similarity' for SearchResult." in str(exc)


def test_PrefetchResult():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    ss4763 = ss4763.to_mutable()
    ss4763.filename = ss4763_file

    scaled = ss47.minhash.scaled

    intersect_mh = ss47.minhash.intersection(ss4763.minhash)
    intersect_bp = len(intersect_mh) * scaled
    jaccard=ss4763.jaccard(ss47)
    max_containment=ss4763.max_containment(ss47)
    f_match_query=ss47.contained_by(ss4763)
    f_query_match=ss4763.contained_by(ss47)
    queryc_ani = ss47.containment_ani(ss4763)
    matchc_ani = ss4763.containment_ani(ss47)

    res = PrefetchResult(ss47, ss4763, cmp_scaled = scaled)

    assert res.query_name == ss47.name
    assert res.match_name == ss4763.name
    assert res.query_scaled == ss47.minhash.scaled == 1000
    assert res.match_scaled == ss4763.minhash.scaled == 1000
    assert res.cmp_scaled == 1000
    assert res.query_abundance == ss47.minhash.track_abundance
    assert res.match_abundance == ss4763.minhash.track_abundance
    assert res.query_bp == len(ss47.minhash) * scaled
    assert res.match_bp == len(ss4763.minhash) * scaled
    assert res.ksize == 31
    assert res.moltype == 'DNA'
    assert res.query_filename == '47.fa'
    assert res.match_filename == ss4763_file
    assert res.query_md5 == ss47.md5sum()
    assert res.match_md5 == ss4763.md5sum()
    assert res.query_n_hashes == len(ss47.minhash)
    assert res.match_n_hashes == len(ss4763.minhash)
    assert res.md5 == ss4763.md5sum()
    assert res.name == ss4763.name
    assert res.intersect_bp == intersect_bp
    assert res.jaccard == jaccard
    assert res.max_containment == max_containment
    assert res.f_query_match == f_query_match
    assert res.f_match_query == f_match_query

    # check ani
    assert res.query_containment_ani == queryc_ani.ani
    print(queryc_ani.ani)
    print(matchc_ani.ani)
    assert res.match_containment_ani == matchc_ani.ani
    assert res.max_containment_ani == max(queryc_ani.ani, matchc_ani.ani)
    assert res.average_containment_ani == np.mean([queryc_ani.ani, matchc_ani.ani])
    assert res.potential_false_negative == False


def test_PrefetchResult_incompatible_sigs():
    ss47_file = utils.get_test_data('num/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    with pytest.raises(TypeError) as exc:
        PrefetchResult(ss47, ss4763)
    print(str(exc))
    assert "Error: prefetch and gather results must be between scaled signatures." in str(exc)


def test_GatherResult():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('track_abund/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    ss4763 = ss4763.to_mutable()
    ss4763.filename = ss4763_file

    scaled = ss47.minhash.scaled

    intersect_mh = ss47.minhash.flatten().intersection(ss4763.minhash)
    remaining_mh = ss4763.minhash.to_mutable()
    remaining_mh.remove_many(intersect_mh)

    intersect_bp = len(intersect_mh) * scaled
    max_containment=ss4763.max_containment(ss47)
    f_match_query = ss47.contained_by(ss4763)
    orig_query_abunds = ss47.minhash.hashes
    queryc_ani = ss47.containment_ani(ss4763)
    matchc_ani = ss4763.containment_ani(ss47)

    # make some fake vals to check
    gather_result_rank = 1
    sum_abunds = 1000

    res = GatherResult(ss47, ss4763, cmp_scaled=scaled,
                        gather_querymh=remaining_mh,
                        gather_result_rank=gather_result_rank,
                        total_weighted_hashes = sum_abunds,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)

    assert res.query_name == ss47.name
    assert res.match_name == ss4763.name
    assert res.query_scaled == ss47.minhash.scaled == 1000
    assert res.match_scaled == ss4763.minhash.scaled == 1000
    assert res.cmp_scaled == 1000
    assert res.query_abundance == ss47.minhash.track_abundance
    assert res.match_abundance == ss4763.minhash.track_abundance
    assert res.query_bp == len(ss47.minhash) * scaled
    assert res.match_bp == len(ss4763.minhash) * scaled
    assert res.ksize == 31
    assert res.moltype == 'DNA'
    assert res.query_filename == 'podar-ref/47.fa'
    assert res.match_filename == ss4763_file
    assert res.query_md5 == ss47.md5sum()
    assert res.match_md5 == ss4763.md5sum()
    assert res.query_n_hashes == len(ss47.minhash)
    assert res.match_n_hashes == len(ss4763.minhash)
    assert res.query_bp == ss47.minhash.unique_dataset_hashes
    assert res.match_bp == ss4763.minhash.unique_dataset_hashes
    assert res.md5 == ss4763.md5sum()
    assert res.name == ss4763.name
    assert res.match_filename == ss4763.filename
    assert res.intersect_bp == intersect_bp
    assert res.max_containment == max_containment

    # check that we can write prefetch result directly from gather
    pf = PrefetchResult(ss47, ss4763, cmp_scaled=scaled)
    assert pf.prefetchresultdict == res.prefetchresultdict

    # check ani
    assert res.query_containment_ani == queryc_ani.ani
    print(queryc_ani.ani)
    print(matchc_ani.ani)
    assert res.match_containment_ani == matchc_ani.ani
    assert res.max_containment_ani == max(queryc_ani.ani, matchc_ani.ani)
    assert res.average_containment_ani == np.mean([queryc_ani.ani, matchc_ani.ani])
    assert res.potential_false_negative == False

    # get write dict version of GatherResult
    resD = res.gatherresultdict
    assert resD["intersect_bp"] == res.intersect_bp


def test_GatherResult_ci():
    # check that values get stored/calculated correctly
    ss47_file = utils.get_test_data('track_abund/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')

    ss4763 = ss4763.to_mutable()
    ss4763.filename = ss4763_file

    scaled = ss47.minhash.scaled

    intersect_mh = ss47.minhash.flatten().intersection(ss4763.minhash)
    remaining_mh = ss4763.minhash.to_mutable()
    remaining_mh.remove_many(intersect_mh)

    orig_query_abunds = ss47.minhash.hashes
    queryc_ani = ss47.containment_ani(ss4763,estimate_ci=True)
    matchc_ani = ss4763.containment_ani(ss47, estimate_ci=True)

    # make some fake vals to check
    gather_result_rank = 1
    sum_abunds = 1000

    res = GatherResult(ss47, ss4763, cmp_scaled=scaled,
                        gather_querymh=remaining_mh,
                        gather_result_rank=gather_result_rank,
                        total_weighted_hashes = sum_abunds,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds,
                        estimate_ani_ci=True)

    # check that we can write prefetch result directly from gather
    pf = PrefetchResult(ss47, ss4763, cmp_scaled=scaled, estimate_ani_ci=True)
    assert pf.prefetchresultdict == res.prefetchresultdict

    # check ani
    assert res.query_containment_ani == queryc_ani.ani
    print(queryc_ani.ani)
    print(matchc_ani.ani)
    assert res.match_containment_ani == matchc_ani.ani
    assert res.match_containment_ani_low == matchc_ani.ani_low
    assert res.match_containment_ani_high == matchc_ani.ani_high
    assert res.max_containment_ani == max(queryc_ani.ani, matchc_ani.ani)
    assert res.average_containment_ani == np.mean([queryc_ani.ani, matchc_ani.ani])
    assert res.potential_false_negative == False

    # get write dict version of GatherResult
    resD = res.gatherresultdict
    assert resD["intersect_bp"] == res.intersect_bp
    assert resD["match_containment_ani_low"] == res.match_containment_ani_low


def test_GatherResult_incompatible_sigs():
    ss47_file = utils.get_test_data('num/47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = ss47.minhash.hashes

    with pytest.raises(TypeError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: prefetch and gather results must be between scaled signatures." in str(exc)


def test_GatherResult_incomplete_input_cmpscaled():
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = ss47.minhash.hashes

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=None,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide comparison scaled value ('cmp_scaled') for GatherResult" in str(exc)


def test_GatherResult_incomplete_input_gathermh():
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = ss47.minhash.hashes

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=None,
                        gather_result_rank=1,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide current gather sketch (remaining hashes) for GatherResult" in str(exc)


def test_GatherResult_incomplete_input_gather_result_rank():
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = ss47.minhash.hashes

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=None,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide 'gather_result_rank' to GatherResult" in str(exc)


def test_GatherResult_incomplete_input_total_weighted_hashes():
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = ss47.minhash.hashes

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = None,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide sum of all abundances ('total_weighted_hashes') to GatherResult" in str(exc)

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = 0,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide sum of all abundances ('total_weighted_hashes') to GatherResult" in str(exc)


def test_GatherResult_incomplete_input_orig_query_abunds():
    ss47_file = utils.get_test_data('47.fa.sig')
    ss4763_file = utils.get_test_data('47+63.fa.sig')
    ss47 = load_one_signature(ss47_file, ksize=31, select_moltype='dna')
    ss4763 = load_one_signature(ss4763_file, ksize=31, select_moltype='dna')
    orig_query_abunds = None

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide original query abundances ('orig_query_abunds') to GatherResult" in str(exc)

    orig_query_abunds = {}

    with pytest.raises(ValueError) as exc:
        GatherResult(ss47, ss4763, cmp_scaled=1000,
                        gather_querymh=ss47.minhash,
                        gather_result_rank=1,
                        total_weighted_hashes = 1,
                        orig_query_len=len(ss47.minhash),
                        orig_query_abunds=orig_query_abunds)
    print(str(exc))
    assert "Error: must provide original query abundances ('orig_query_abunds') to GatherResult" in str(exc)
