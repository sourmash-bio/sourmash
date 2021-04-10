"Tests for search.py code."
from sourmash import search
from sourmash.search import make_jaccard_search_query

def test_make_jaccard_search_query():
    search_obj = make_jaccard_search_query(threshold=0)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_no_threshold_none():
    search_obj = make_jaccard_search_query(threshold=None)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0
