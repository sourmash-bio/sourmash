"""
Tests for functions in taxonomy submodule.
"""
import pytest

import sourmash
import sourmash_tst_utils as utils

from sourmash.tax import tax_utils
from sourmash.tax.tax_utils import (ascending_taxlist, get_ident, load_gather_results,
                                    summarize_gather_at, find_missing_identities)#,
                                    #gather_at_rank)

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import LineagePair#, build_tree, find_lca,
#                                    taxlist, count_lca_for_assignments,
#                                    zip_lineage, display_lineage,
#                                    make_lineage, is_lineage_match,
#                                    pop_to_rank)

# utility functions for testing
def make_mini_gather_results(g_infolist):
    # make mini gather_results
    min_header = ["name","match_ident","f_unique_weighted"]
    gather_results = []
    for g_info in g_infolist:
        inf = dict(zip(min_header, g_info))
        gather_results.append(inf)
    return gather_results

def make_mini_taxonomy(tax_info):
    #pass in list of tuples: (name, lineage)
    taxD = {}
    for (name,lin) in tax_info:
        taxD[name] = lca_utils.make_lineage(lin)
    return taxD

## tests
def test_ascending_taxlist_1():
    assert list(ascending_taxlist()) ==  ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

def test_ascending_taxlist_2():
    assert list(ascending_taxlist(include_strain=False)) ==  ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def test_get_ident():
    ident = "GCF_001881345.1"
    n_id = tax_utils.get_ident(ident)
    assert n_id == "GCF_001881345"

def test_load_gather_results():
    gather_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results = tax_utils.load_gather_results([gather_csv])
    assert len(gather_results) == 4

def test_find_missing_identities():
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    n, ids = find_missing_identities(g_res, taxD)
    print("n_missing: ", n)
    print("ids_missing: ", ids)
    assert n == 1
    assert ids == ["gB"]


def test_summarize_gather_at_0():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 1.0)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.0)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.5)]

def test_summarize_gather_at_1():
    """test two matches, diff f_unique_weighted"""
    # make mini gather_results
    gA = ["gA","0.5","0.6"]
    gB = ["gB","0.3","0.1"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 0.7)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),0.7)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.6),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.1)]

def test_summarize_gather_at_over100percent_f_unique_weighted():
    """gather matches that add up to >100% f_unique_weighted"""
    ## @NTP:  currently passes, we should probably make this fail
    # make mini gather_results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.6"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 1.1)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.1)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.6),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5)]

def test_summarize_gather_at_best_only_0():
    """test two matches, diff f_unique_weighted"""
    # make mini gather_results
    gA = ["gA","0.5","0.6"]
    gB = ["gB","0.3","0.1"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res, best_only=True)
    assert sk_sum == ((LineagePair(rank='superkingdom', name='a'),), 0.7)
    phy_sum = summarize_gather_at("phylum", taxD, g_res, best_only=True)
    assert phy_sum == ((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),0.7)
    cl_sum = summarize_gather_at("class", taxD, g_res, best_only=True)
    assert cl_sum == ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.6)

def test_summarize_gather_at_best_only_equal_choose_first():
    """test two matches, equal f_unique_weighted. best_only chooses first"""
    # make mini gather_results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res, best_only=True)
    assert sk_sum == ((LineagePair(rank='superkingdom', name='a'),), 1.0)
    phy_sum = summarize_gather_at("phylum", taxD, g_res, best_only=True)
    assert phy_sum == ((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.0)
    cl_sum = summarize_gather_at("class", taxD, g_res, best_only=True)
    assert cl_sum == ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5)
