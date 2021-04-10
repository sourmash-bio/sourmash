"Tests for search.py code."
import pytest

from sourmash import search, SourmashSignature, MinHash
from sourmash.search import make_jaccard_search_query

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
    search_obj.collect(1.0)
    assert search_obj.threshold == 0


def test_collect_best_only():
    search_obj = make_jaccard_search_query(threshold=0, best_only=True)
    search_obj.collect(1.0)
    assert search_obj.threshold == 1.0
