"""
Tests for functions in taxonomy submodule.
"""
import pytest

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, build_tree, find_lca,
                                    taxlist, count_lca_for_assignments,
                                    zip_lineage, display_lineage,
                                    make_lineage, is_lineage_match,
                                    pop_to_rank)

# some utility functions for testing
def make_mh(hashvals, ksize=3, scaled=1):
    mh = sourmash.MinHash(n=0, scaled=1, ksize=3)
    mh.add_many(hashvals)
    return mh


def make_sig_and_lin(hashvals, ident, lin, ksize=3, scaled=1):
    mh = make_mh(hashvals)
    sig = sourmash.SourmashSignature(mh, name=ident)
    lineage = lca_utils.make_lineage(lin)
    return mh, sig, lineage

def test_gen_mh():
    mh = make_mh([12345678])
    return mh.copy_and_clear()

def test_gather_at_rank_1():
    # one minhash, one set of ranks
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    lin_db = LineageDB()
    lin_db.insert(ident, lin1)

    gather_results=list(gather_at_rank(mh1, lca_db, lin_db, "class"))
    assert len(gather_results) == 1
    assert gather_results[0][0] == lin1
    assert gather_results[0][1] == 1


def test_gather_at_rank_2():
   #two minhashes, fully shared ranks

    # first sig
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')

    # second sig
    hashval2 = 87654321
    ident2 = 'second'
    mh2, sig2, lin2 = make_sig_and_lin([hashval2], ident2, 'a;b;c')

    # create lca_db w sigs
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident1)
    lca_db.insert(sig2, ident=ident2)

    # make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    lin_db.insert(ident2, lin2)

    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))
    assert len(gather_results) == 1
    assert gather_results[0][0] == lin1
    assert gather_results[0][1] == 2


def test_gather_at_rank_3():
    # two minhashes, totally distinct ranks
    # first sig
    hashval1  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident1, 'a;b;c')

    # second sig
    hashval2 = 87654321
    ident2 = 'second'
    mh2, sig2, lin2 = make_sig_and_lin([hashval2], ident2, 'd;e;f')

    # create lca_db w sig1
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident1)
    lca_db.insert(sig2, ident=ident2)

    # next, make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    lin_db.insert(ident2, lin2)

    # search with combined hashvals
    search_mh = make_mh([hashval1, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))

    assert len(gather_results) == 2
    assert set([gather_results[0][0],gather_results[1][0]]) == set([lin1, lin2])
    assert set([gather_results[0][1],gather_results[1][1]]) == set([1])
