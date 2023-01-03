"""
Tests for the 'TaxComparison' classes.
"""

import pytest
from pytest import approx
#import sourmash_tst_utils as utils

from sourmash.tax.taxcomparison import (LineagePair, LineageTuple, BaseLineageInfo, RankLineageInfo,
                                        SummarizedGatherResult, build_tree, GatherRow, TaxResult, QueryTaxResult,
                                        make_mini_taxonomy, make_GatherRow, make_TaxResult, make_QueryTaxResults) #LINSLineageInfo

# sigh, can't make build tree work as easily with both orig LineagePair and LineageInfo.
# What if LineagePair just had the extra info?
from sourmash.tax.taxcomparison import LineageTuple as LineagePair
from sourmash.lca.lca_utils import find_lca, make_lineage
from sourmash.tax.tax_utils import write_summary, aggregate_by_lineage_at_rank, format_for_krona, write_krona

taxranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def test_LineagePair():
    lin = LineagePair(rank="rank1", name='name1')
    print(lin)


def test_LineageTuple():
    lin1 = LineageTuple(rank="rank1", name='name1', taxid=1)
    print(lin1)


def test_LineageTuple_taxid_as_str():
    with pytest.raises(TypeError) as exc:
        LineageTuple(rank="rank1", name='name1', taxid="1")
    print(exc)
    assert "taxid not an int" in str(exc)


def test_RankLineageInfo_taxlist():
    taxinf = RankLineageInfo()
    assert taxinf.taxlist == taxranks
    assert taxinf.ascending_taxlist == taxranks[::-1] 


def test_RankLineageInfo_init_lineage_str_1():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c', '', '', '', '', '']


def test_BaseLineageInfo_init_lineage_str_1():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    taxinf = BaseLineageInfo(lineage_str=x, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c']


def test_BaseLineageInfo_init_lineage_str_lineage_dict_test_eq():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    rankD = {"A": "a", "B": "b", "C": "c"}
    lin1 = BaseLineageInfo(lineage_str=x, ranks=ranks)
    lin2 = BaseLineageInfo(lineage_dict=rankD, ranks=ranks)
    assert lin1 == lin2 


def test_RankLineageInfo_init_lineage_str_lineage_dict_test_eq():
    x = "a;b;c"
    rankD = {"superkingdom": "a", "phylum": "b", "class": "c"}
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_dict=rankD)
    assert lin1 == lin2 


def test_BaseLineageInfo_reinit_lineage():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    lin1 = BaseLineageInfo(lineage_str=x, ranks=ranks)
    print(lin1)
    with pytest.raises(ValueError) as exc:
        lin1.init_empty()
    assert "lineage not empty" in str(exc)


def test_RankLineageInfo_reinit_lineage():
    x = "a;b;c"
    lin1 = RankLineageInfo(lineage_str=x)
    print(lin1)
    with pytest.raises(ValueError) as exc:
        lin1.init_empty()
    assert "lineage not empty" in str(exc)


def test_BaseLineageInfo_init_lineage_str_no_ranks():
    x = "a;b;c"
    with pytest.raises(ValueError) as exc:
        BaseLineageInfo(lineage_str=x)
    print(exc)
    assert "Cannot initialize BaseLineageInfo. Please provide lineage or rank info." in str(exc)


def test_RankLineageInfo_init_lineage_str_1_truncate():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', 'c']


def test_RankLineageInfo_init_lineage_str_2():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', 'c' '', '', '', '', '']


def test_RankLineageInfo_init_lineage_str_2_truncate():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', '', 'c']


def test_RankLineageInfo_init_lineage_with_incorrect_rank():
    x = [ LineagePair('superkingdom', 'a'), LineagePair("NotARank", ''), LineagePair('class', 'c') ]
    with pytest.raises(ValueError) as exc:
        RankLineageInfo(lineage=x)
    print(str(exc))
    assert f"Rank 'NotARank' not present in " in str(exc)


def test_BaseLineageInfo_init_lineage_dict_1():
    x = {'rank1': 'name1', 'rank2': 'name2'}
    taxinf = BaseLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', 'name2']


def test_BaseLineageInfo_init_lineage_dict_withtaxid():
    x = {'rank1': {'name': 'name1', 'taxid': 1}, 'rank2': {'name':'name2', 'taxid': 2}}
    taxinf = BaseLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', 'name2']
    assert taxinf.zip_taxid()== ['1', '2']


def test_RankLineageInfo_init_lineage_dict_1():
    x = {'superkingdom': 'name1', 'class': 'name2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_withtaxid():
    x = {'superkingdom': {'name': 'name1', 'taxid': 1}, 'class': {'name':'name2', 'taxid': 2}}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']
    assert taxinf.zip_taxid()== ['1', '', '2', '', '', '', '', '']


def test_zip_lineage_1():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage() == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage(truncate_empty=True))
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage() == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_3_truncate():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_display_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.display_lineage() == "a;b"


def test_display_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.display_lineage() == "a;;c"


def test_display_taxid_1():
    x = [ LineageTuple('superkingdom', 'a', 1), LineageTuple('phylum', 'b', 2) ]
    taxinf = RankLineageInfo(lineage=x)
    print(taxinf)
    assert taxinf.display_taxid() == "1;2"

def test_display_taxid_2():
    x = [ LineageTuple('superkingdom', 'name1', 1), LineageTuple(None, ''), LineageTuple    ('class', 'name2',2) ]
    taxinf = RankLineageInfo(lineage=x)
    print(taxinf)
    assert taxinf.display_taxid() == "1;;2"


def test_is_lineage_match_1():
    # basic behavior: match at order and above, but not at family or below.
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    assert lin1.is_lineage_match(lin2, 'superkingdom')
    assert lin2.is_lineage_match(lin1, 'superkingdom')
    assert lin1.is_lineage_match(lin2, 'phylum')
    assert lin2.is_lineage_match(lin1, 'phylum')
    assert lin1.is_lineage_match(lin2, 'class')
    assert lin2.is_lineage_match(lin1, 'class')
    assert lin1.is_lineage_match(lin2, 'order')
    assert lin2.is_lineage_match(lin1, 'order')
    
    assert not lin1.is_lineage_match(lin2, 'family')
    assert not lin2.is_lineage_match(lin1, 'family')
    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_2():
    # match at family, and above, levels; no genus or species to match
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    assert lin1.is_lineage_match(lin2, 'superkingdom')
    assert lin2.is_lineage_match(lin1, 'superkingdom')
    assert lin1.is_lineage_match(lin2, 'phylum')
    assert lin2.is_lineage_match(lin1, 'phylum')
    assert lin1.is_lineage_match(lin2, 'class')
    assert lin2.is_lineage_match(lin1, 'class')
    assert lin1.is_lineage_match(lin2, 'order')
    assert lin2.is_lineage_match(lin1, 'order')
    assert lin1.is_lineage_match(lin2, 'family')
    assert lin2.is_lineage_match(lin1, 'family')

    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_3():
    # one lineage is empty
    lin1 = RankLineageInfo()
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    
    assert not lin1.is_lineage_match(lin2, 'superkingdom')
    assert not lin2.is_lineage_match(lin1, 'superkingdom')
    assert not lin1.is_lineage_match(lin2, 'phylum')
    assert not lin2.is_lineage_match(lin1, 'phylum')
    assert not lin1.is_lineage_match(lin2, 'class')
    assert not lin2.is_lineage_match(lin1, 'class')
    assert not lin1.is_lineage_match(lin2, 'order')
    assert not lin2.is_lineage_match(lin1, 'order')
    assert not lin1.is_lineage_match(lin2, 'family')
    assert not lin2.is_lineage_match(lin1, 'family')
    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_incorrect_ranks():
    #test comparison with incompatible ranks
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e', ranks=taxranks[::-1])
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    with pytest.raises(ValueError) as exc:
        lin1.is_lineage_match(lin2, 'superkingdom')
    print(str(exc))
    assert 'Cannot compare lineages from taxonomies with different ranks.' in str(exc)


def test_pop_to_rank_1():
    # basic behavior - pop to order?
    lin1 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d')
    lin2 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d;f__f')

    print(lin1)
    popped = lin2.pop_to_rank('order')
    print(popped)
    assert popped == lin1


def test_pop_to_rank_2():
    # what if we're already above rank?
    lin2 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d;f__f')
    print(lin2.pop_to_rank('species'))
    assert lin2.pop_to_rank('species') == lin2


def test_pop_to_rank_rank_not_avail():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(ValueError) as exc:
        lin1.pop_to_rank("NotARank")
    print(str(exc))
    assert "Desired Rank 'NotARank' not available for this lineage" in str(exc)
    


def test_build_tree_base():
    lin1 = BaseLineageInfo(lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2')])
    print(lin1)
    tree = build_tree(lin1)
    assert tree == { LineagePair('rank1', 'name1'):
                         { LineagePair('rank2', 'name2') : {}} }


def test_build_tree_1():
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2')])
    print(lin1)
    tree = build_tree(lin1)
    assert tree == { LineagePair('rank1', 'name1'):
                         { LineagePair('rank2', 'name2') : {}} }


def test_build_tree_2():
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2a')])
    lin2 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2b')])
    tree = build_tree([lin1,lin2])

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }

def test_build_tree_2_LineagePairs():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }


def test_build_tree_3():
    # empty 'rank2' name
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', '')])
    tree = build_tree(lin1)
    assert tree == { LineagePair('rank1', 'name1'): {} }


def test_build_tree_3_LineagePairs():
    # empty 'rank2'
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"], 
                           lineage=[LineagePair('rank1', 'name1')])
    tree = build_tree(lin1)
    assert tree == { LineagePair('rank1', 'name1'): {} }


def test_build_tree_4():
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2a')])
    lin2 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2b')])
    tree = build_tree(lin1)

    tree = build_tree(lin2, tree)

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }

def test_build_tree_4_LineagePairs():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                      ])

    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ], tree)

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }

def test_build_tree_5():
    with pytest.raises(ValueError):
        tree = build_tree([])

def test_build_tree_5b():
    with pytest.raises(ValueError):
        tree = build_tree("")

def test_find_lca():
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2')])
    tree = build_tree(lin1)
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2'),), 0)


def test_find_lca_LineagePairs():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2'),), 0)


def test_find_lca_2():
    lin1 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2a')]) 
    lin2 = RankLineageInfo(ranks = ['rank1', "rank2"],
                           lineage=[LineagePair('rank1', 'name1'),
                                    LineagePair('rank2', 'name2b')]) 

    tree = build_tree([lin1, lin2])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'),), 2)


def test_find_lca_2_LineagePairs():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'),), 2)


def test_find_lca_3():
    lin1 = RankLineageInfo(lineage_str="a;b;c")
    lin2 = RankLineageInfo(lineage_str="a;b")

    tree = build_tree([lin1, lin2])
    lca, reason = find_lca(tree)
    assert lca == tuple(lin1.filled_lineage)           # find most specific leaf node


def test_lineage_at_rank_norank():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(TypeError) as exc:
        lin1.lineage_at_rank()
    print(str(exc))
    assert "lineage_at_rank() missing 1 required positional argument: 'rank'" in str(exc)


def test_lineage_at_rank_rank_not_avail():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(ValueError) as exc:
        lin1.lineage_at_rank("NotARank")
    print(str(exc))
    assert "Desired Rank 'NotARank' not available for this lineage" in str(exc)
    

def test_lineage_at_rank_1():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage_at_rank('superkingdom'))
    
    assert lin1.lineage_at_rank('superkingdom') == [LineageTuple(rank='superkingdom', name='d__a', taxid=None)]
    print(lin1.lineage_at_rank('class'))
    assert lin1.lineage_at_rank('class') == [LineageTuple(rank='superkingdom', name='d__a', taxid=None), 
                                             LineageTuple(rank='phylum', name='p__b', taxid=None), 
                                             LineageTuple(rank='class', name='c__c', taxid=None)]
    

def test_lineage_at_rank_below_rank():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage_at_rank('superkingdom'))
    # if rank is not provided, we only return the filled lineage, to follow original pop_to_rank behavior.

    print(lin1.lineage_at_rank('genus'))
    assert lin1.lineage_at_rank('genus') == [LineageTuple(rank='superkingdom', name='d__a', taxid=None), 
                                             LineageTuple(rank='phylum', name='p__b', taxid=None), 
                                             LineageTuple(rank='class', name='c__c', taxid=None),
                                             LineageTuple(rank='order', name='o__d', taxid=None), 
                                             LineageTuple(rank='family', name='f__f', taxid=None)]



def test_TaxResult_old_gather():
    # gather does not contain query_name column
    gA = {"name": "gA.1 name"}
    gRow = make_GatherRow(gA)
    gRow.query_bp = None  # reset query_bp
    print("query_bp: ", gRow.query_bp)
    with pytest.raises(ValueError) as exc:
        TaxResult(raw=gRow)
    print(str(exc))
    assert "Error: Please run gather with sourmash >= 4.4 for taxonomic summarization." in str(exc)


def test_TaxResult_get_ident_1():
    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    print(taxres.match_ident)
    assert taxres.match_ident == "gA"


def test_TaxResult_get_ident_keep_full():
    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA, keep_full_ident=True)
    print("raw ident: ", taxres.raw.name)
    print("keep_full?: ", taxres.keep_full_identifiers)
    print("keep_version?: ",taxres.keep_identifier_versions)
    print("final ident: ", taxres.match_ident)
    assert taxres.match_ident == "gA.1 name"


def test_TaxResult_get_ident_keep_version():
    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA, keep_ident_version=True)
    print("raw ident: ", taxres.raw.name)
    print("keep_full?: ", taxres.keep_full_identifiers)
    print("keep_version?: ",taxres.keep_identifier_versions)
    print("final ident: ", taxres.match_ident)
    assert taxres.match_ident == "gA.1"


def test_TaxResult_get_match_lineage_1():
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD)
    assert taxres.lineageInfo.display_lineage() == "a;b;c"


def test_TaxResult_get_match_lineage_skip_ident():
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gA'])
    print("skipped_ident?: ", taxres.skipped_ident)
    print("missed_ident?: ", taxres.missed_ident)
    assert taxres.skipped_ident == True
    assert taxres.lineageInfo.display_lineage() == ""


def test_TaxResult_get_match_lineage_missed_ident():
    gA_tax = ("gA.1", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gB'])
    print("skipped_ident?: ", taxres.skipped_ident)
    print("missed_ident?: ", taxres.missed_ident)
    assert taxres.skipped_ident == False
    assert taxres.missed_ident == True
    assert taxres.lineageInfo.display_lineage() == ""


def test_TaxResult_get_match_lineage_missed_ident_fail_on_missing():
    gA_tax = ("gA.1", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    with pytest.raises(ValueError) as exc:
        taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gB'], fail_on_missing_taxonomy=True)
    print(str(exc))
    assert "Error: ident 'gA' is not in the taxonomy database." in str(exc)


def test_QueryTaxResult_1():
    "basic functionality: initialize and add a taxresult"
    tax_info = [("gA", "a;b;c")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    taxres = make_TaxResult(taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    q_res.add_taxresult(taxres)
    # check that new querytaxres is compatible with taxres
    assert q_res.is_compatible(taxres)
    # check that a few thngs were set properly and/or are not yet set.
    assert q_res.query_name == "q1"
    assert q_res.query_info.query_bp == 100
    assert len(q_res.raw_taxresults) == 1
    assert q_res.skipped_idents == set()
    assert q_res.missed_idents == set()
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_add_incompatible():
    "basic functionality: initialize and add a taxresult"
    tax_info = [("gA", "a;b;c")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    taxres = make_TaxResult(taxD=taxD)
    taxres2 = make_TaxResult({'query_name': 'q2'}, taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # check that new querytaxres is compatible with taxres and not taxres2
    assert q_res.is_compatible(taxres)
    assert not q_res.is_compatible(taxres2)
    q_res.add_taxresult(taxres)
    with pytest.raises(ValueError) as exc:
        q_res.add_taxresult(taxres2)
    print(str(exc))
    assert "Error: Cannot add TaxResult: query information does not match." in str(exc)


def test_QueryTaxResult_add_without_tax_info():
    "initialize and add a taxresult with missed ident"
    taxres = make_TaxResult() # do not add taxonomic info
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    print("attempted to add lineage info?: ", taxres.match_lineage_attempted)
    with pytest.raises(ValueError) as exc:
        q_res.add_taxresult(taxres)
    print(str(exc))
    assert "Error: Cannot add TaxResult. Please use get_match_lineage() to add taxonomic lineage information first." in str(exc)
    
    
def test_QueryTaxResult_add_skipped_ident():
    "initialize and add a taxresult with skipped ident"
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])
    taxres = make_TaxResult(taxD=taxD, skip_idents = ['gA'])
#    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gA'])
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    q_res.add_taxresult(taxres)
    assert len(q_res.skipped_idents) == 1
    assert len(q_res.raw_taxresults) == 1
    assert q_res.missed_idents == set()
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_add_missed_ident():
    "initialize and add a taxresult with missed ident"
    gA_tax = ("gB", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])
    taxres = make_TaxResult(taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # add taxonomic info to taxres
    q_res.add_taxresult(taxres)
    assert len(q_res.missed_idents) == 1
    assert len(q_res.raw_taxresults) == 1
    assert q_res.skipped_idents == set()
    assert q_res.summarized_lineage_results == {} 


def test_QueryTaxResult_track_missed_and_skipped():
    "make sure missed and skipped idents are being tracked"
    # make taxonomy
    tax_info = [("gA", "a;b;c"), ("gB", "a;b;d")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    # make results
    taxres = make_TaxResult()
    taxres2 = make_TaxResult({"name": 'gB'}) # skipped
    taxres3 = make_TaxResult({"name": 'gB'}) # skipped
    taxres4 = make_TaxResult({"name": 'gC'}) # skipped
    taxres5 = make_TaxResult({"name": 'gD'}) # missed
    taxres6 = make_TaxResult({"name": 'gE'}) # missed
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # add taxonomic info to taxres, add to q_res
    for n, tr in enumerate([taxres, taxres2, taxres3, taxres4, taxres5, taxres6]):
        tr.get_match_lineage(tax_assignments=taxD, skip_idents=['gB', 'gC'])
        print("num: ", n)
        print("skipped?: ", tr.skipped_ident)
        print("missed?: ", tr.missed_ident)
        q_res.add_taxresult(tr)
    assert len(q_res.raw_taxresults) == 6
    print(q_res.n_skipped)
    print(q_res.n_missed)
    assert q_res.n_missed == 2
    assert q_res.n_skipped == 3
    assert 'gB' in q_res.skipped_idents
    assert len(q_res.skipped_idents) == 2
    assert 'gD' in q_res.missed_idents
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_track_missed_and_skipped_using_fn():
    "make sure missed and skipped idents are being tracked. Same as above but use helper fn."
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}, {"name": 'gB'}, {"name": 'gC'}, {"name": 'gD'}, {"name": 'gE'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, skip_idents=['gB', 'gC'])
    # should have 6 results for default query 'q1' 
    print(gres.keys())
    q_res = next(iter(gres.values()))
    assert len(q_res.raw_taxresults) == 6
    print(q_res.n_skipped)
    print(q_res.n_missed)
    assert q_res.n_missed == 2
    assert q_res.n_skipped == 3
    assert 'gB' in q_res.skipped_idents
    assert len(q_res.skipped_idents) == 2
    assert 'gD' in q_res.missed_idents
    assert q_res.summarized_lineage_results == {}
    

def test_QueryTaxResult_summarize_up_ranks_1():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    #print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): 0.4}
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.4}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.2}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 40}
    print(q_res.sum_uniq_weighted['class'])
    assert q_res.sum_uniq_weighted['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.2, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='d', taxid=None)): 0.2}
    assert q_res.sum_uniq_to_query['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.1, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='d', taxid=None)): 0.1}
    assert q_res.sum_uniq_bp['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 20, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='d', taxid=None)): 20}


def test_QueryTaxResult_summarize_up_ranks_2():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): approx(0.3)}
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.3)}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.15)}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): 30}
    print(q_res.sum_uniq_weighted['class'])
    assert q_res.sum_uniq_weighted['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): approx(0.2), 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): approx(0.1)}
    assert q_res.sum_uniq_to_query['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): 0.1, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): 0.05}
    assert q_res.sum_uniq_bp['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): 20, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): 10}


def test_QueryTaxResult_summarize_up_ranks_missing_lineage():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    #print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): 0.2}
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.2}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.1}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 20}
    print(q_res.sum_uniq_weighted['class'])
    assert q_res.sum_uniq_weighted['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.2}
    assert q_res.sum_uniq_to_query['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.1}
    assert q_res.sum_uniq_bp['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 20}


def test_QueryTaxResult_summarize_up_ranks_skipped_lineage():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, skip_idents=['gB'])
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    #print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): 0.2}
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.2}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 0.1}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                  LineageTuple(rank='phylum', name='b', taxid=None)): 20}
    print(q_res.sum_uniq_weighted['class'])
    assert q_res.sum_uniq_weighted['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.2}
    assert q_res.sum_uniq_to_query['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 0.1}
    assert q_res.sum_uniq_bp['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                 LineageTuple(rank='phylum', name='b', taxid=None), 
                                                 LineageTuple(rank='class', name='c', taxid=None)): 20}


def test_QueryTaxResult_summarize_up_ranks_perfect_match():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{'f_unique_to_query': 1.0}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 1
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_to_query['superkingdom'])
    assert q_res.sum_uniq_to_query['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): approx(1.0)}
    assert 'gA' in q_res.perfect_match


def test_QueryTaxResult_summarize_up_ranks_already_summarized():
    "summarize up ranks: error, already summarized"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{'f_unique_to_query': 1.0}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    with pytest.raises(ValueError) as exc:
        q_res.summarize_up_ranks()
    print(str(exc))
    assert "Error: already summarized" in str(exc)
    

def test_QueryTaxResult_summarize_up_ranks_already_summarized_force():
    "summarize up ranks: already summarized but force"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    q_res.summarize_up_ranks(force_resummarize=True)

    #check that all results are still good 
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == q_res.ranks[::-1]
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {(LineageTuple(rank='superkingdom', name='a', taxid=None),): approx(0.3)}
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.3)}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.15)}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): 30}
    print(q_res.sum_uniq_weighted['class'])
    assert q_res.sum_uniq_weighted['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): approx(0.2), 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): approx(0.1)}
    assert q_res.sum_uniq_to_query['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): 0.1, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): 0.05}
    assert q_res.sum_uniq_bp['class'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='c', taxid=None)): 20, 
                                                (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                                LineageTuple(rank='phylum', name='b', taxid=None), 
                                                LineageTuple(rank='class', name='d', taxid=None)): 10}
                            

def test_QueryTaxResult_summarize_up_ranks_single_rank():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks(single_rank='phylum')
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == ['phylum']
    print(q_res.sum_uniq_weighted.keys())
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.3)}
    assert q_res.sum_uniq_to_query['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): approx(0.15)}
    assert q_res.sum_uniq_bp['phylum'] == {(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                                LineageTuple(rank='phylum', name='b', taxid=None)): 30}
    assert q_res.summarized_ranks == ['phylum']


def test_QueryTaxResult_build_summarized_result_1():
    "basic functionality: build summarized_result"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    sk = [SummarizedGatherResult(query_name='q1',query_md5='md5', query_filename='query_fn', 
                             rank='superkingdom', fraction=0.2, f_weighted_at_rank=0.4, 
                             lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),), 
                             bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2),
                             total_weighted_hashes=0),
          SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                             rank='superkingdom', fraction=0.8, f_weighted_at_rank=0.6,
                             lineage=(), bp_match_at_rank=60, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                  rank='phylum', fraction=0.2, f_weighted_at_rank=0.4, 
                                  lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                           LineageTuple(rank='phylum', name='b', taxid=None)),
                                  bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2),
                                  total_weighted_hashes=0),
           SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                  rank='phylum', fraction=0.8, f_weighted_at_rank=0.6,
                                  lineage=(), bp_match_at_rank=60, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                 rank='class', fraction=0.1, f_weighted_at_rank=0.2, 
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                          LineageTuple(rank='phylum', name='b', taxid=None), 
                                          LineageTuple(rank='class', name='c', taxid=None)), 
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.93, rel=1e-2),
                                  total_weighted_hashes=0),
          SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                 rank='class', fraction=0.1, f_weighted_at_rank=0.2,
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),
                                          LineageTuple(rank='phylum', name='b', taxid=None),
                                          LineageTuple(rank='class', name='d', taxid=None)),
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.93, rel=1e-2),
                                  total_weighted_hashes=0),
          SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                 rank='class', fraction=0.8, f_weighted_at_rank=0.6,
                                 lineage=(), bp_match_at_rank=60, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.4)
    assert q_res.total_f_classified['class'] == approx(0.2)
    assert q_res.total_bp_classified['superkingdom'] == 40
    assert q_res.best_only_result == False


def test_QueryTaxResult_build_summarized_result_2():
    """test two queries, build summarized result for each"""
    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax, gB_tax])
    # make gather results
    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.5,'f_unique_to_query': 0.5,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.4,'f_unique_to_query': 0.3,'unique_intersect_bp': 30},
                      {'query_name': 'queryB', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    
    for query_name, q_res in gres.items():
        q_res.build_summarized_result() # summarize and build result
        sk = q_res.summarized_lineage_results['superkingdom']
        phy = q_res.summarized_lineage_results['phylum']
        assert len(sk) == 2
        assert sk[0].lineage == (LineageTuple(rank='superkingdom', name='a', taxid=None),)
        print(phy)
        if query_name == 'queryA':
            # check superkingdom results
            assert sk[0].fraction == approx(0.8)
            assert sk[0].f_weighted_at_rank == approx(0.9)
            assert sk[0].bp_match_at_rank == 80
            assert sk[1].fraction == approx(0.2)
            assert sk[1].f_weighted_at_rank == approx(0.1)
            assert sk[1].bp_match_at_rank == 20
            assert sk[1].lineage == ()
            # check phylum results
            assert len(phy) == 3
            assert phy[0].fraction == approx(0.5)
            assert phy[0].f_weighted_at_rank == approx(0.5)
            assert phy[0].bp_match_at_rank == 50
            assert phy[0].lineage == (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                      LineageTuple(rank='phylum', name='b', taxid=None))
            assert phy[1].fraction == approx(0.3)
            assert phy[1].f_weighted_at_rank == approx(0.4)
            assert phy[1].bp_match_at_rank == 30
            assert phy[1].lineage == (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                      LineageTuple(rank='phylum', name='c', taxid=None))
            assert phy[2].fraction == approx(0.2)
            assert phy[2].f_weighted_at_rank == approx(0.1)
            assert phy[2].bp_match_at_rank == 20
            assert phy[2].lineage == ()
        if query_name == 'queryB':
            # check superkingdom results
            assert sk[0].fraction == approx(0.3)
            assert sk[0].f_weighted_at_rank == approx(0.3)
            assert sk[0].bp_match_at_rank == 30
            assert sk[1].fraction == approx(0.7)
            assert sk[1].f_weighted_at_rank == approx(0.7)
            assert sk[1].bp_match_at_rank == 70
            assert sk[1].lineage == ()
            # check phylum results
            assert len(phy) == 2
            assert phy[0].fraction == approx(0.3)
            assert phy[0].f_weighted_at_rank == approx(0.3)
            assert phy[0].bp_match_at_rank == 30
            assert phy[0].lineage == (LineageTuple(rank='superkingdom', name='a', taxid=None),
                                      LineageTuple(rank='phylum', name='c', taxid=None))
            assert phy[1].fraction == approx(0.7)
            assert phy[1].f_weighted_at_rank == approx(0.7)
            assert phy[1].bp_match_at_rank == 70
            assert phy[1].lineage == ()


def test_QueryTaxResult_build_summarized_result_missing_lineage():
    "build summarized_result with missing lineage"
    taxD = make_mini_taxonomy([("gA", "a;b;c")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    print(q_res.summarized_lineage_results['superkingdom'])

    sk = [SummarizedGatherResult(query_name='q1', rank='superkingdom', fraction=0.1, 
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),), 
                                 query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2, 
                                 bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), 
                                 total_weighted_hashes=0), 
          SummarizedGatherResult(query_name='q1', rank='superkingdom', fraction=0.9, lineage=(),
                                 query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.8,
                                 bp_match_at_rank=80, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(query_name='q1', rank='phylum', fraction=0.1, 
                                  lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                           LineageTuple(rank='phylum', name='b', taxid=None)),
                                  query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2,
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), total_weighted_hashes=0),
           SummarizedGatherResult(query_name='q1', rank='phylum', fraction=0.9, lineage=(), query_md5='md5',
                                  query_filename='query_fn', f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                  query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(query_name='q1', rank='class', fraction=0.1,
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),
                                          LineageTuple(rank='phylum', name='b', taxid=None),
                                          LineageTuple(rank='class', name='c', taxid=None)),
                                  query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2,
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), total_weighted_hashes=0),
          SummarizedGatherResult(query_name='q1', rank='class', fraction=0.9, lineage=(), query_md5='md5',
                                 query_filename='query_fn', f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                 query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.2)
    assert q_res.total_f_classified['class'] == approx(0.1)
    assert q_res.total_bp_classified['superkingdom'] == 20
    assert q_res.best_only_result == False


def test_QueryTaxResult_build_summarized_result_skipped_lineage():
    "build summarized_result with skipped lineage"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, skip_idents=['gB'])
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    print(q_res.summarized_lineage_results['superkingdom'])

    sk = [SummarizedGatherResult(query_name='q1', rank='superkingdom', fraction=0.1, 
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),), 
                                 query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2, 
                                 bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), 
                                 total_weighted_hashes=0), 
          SummarizedGatherResult(query_name='q1', rank='superkingdom', fraction=0.9, lineage=(),
                                 query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.8,
                                 bp_match_at_rank=80, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(query_name='q1', rank='phylum', fraction=0.1, 
                                  lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                           LineageTuple(rank='phylum', name='b', taxid=None)),
                                  query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2,
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), total_weighted_hashes=0),
           SummarizedGatherResult(query_name='q1', rank='phylum', fraction=0.9, lineage=(), query_md5='md5',
                                  query_filename='query_fn', f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                  query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(query_name='q1', rank='class', fraction=0.1,
                                 lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),
                                          LineageTuple(rank='phylum', name='b', taxid=None),
                                          LineageTuple(rank='class', name='c', taxid=None)),
                                  query_md5='md5', query_filename='query_fn', f_weighted_at_rank=0.2,
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2), total_weighted_hashes=0),
          SummarizedGatherResult(query_name='q1', rank='class', fraction=0.9, lineage=(), query_md5='md5',
                                 query_filename='query_fn', f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                 query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.2)
    assert q_res.total_f_classified['class'] == approx(0.1)
    assert q_res.total_bp_classified['superkingdom'] == 20
    assert q_res.best_only_result == False


def test_QueryTaxResult_build_summarized_result_over100percent():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_to_query': 0.95}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    assert len(q_res.raw_taxresults) == 2
    with pytest.raises(ValueError) as exc:
        q_res.build_summarized_result()
    print(str(exc))
    assert "The tax summary of query 'q1' is 1.05, which is > 100% of the query!! This should not be possible." in str(exc)


def test_QueryTaxResult_build_summarized_result_best_only():
    "basic functionality: build summarized_result"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_summarized_result(best_only=True)
    print(q_res.summarized_lineage_results.keys())
    print(q_res.summarized_lineage_results['superkingdom'])
    sk = [SummarizedGatherResult(query_name='q1',query_md5='md5', query_filename='query_fn', 
                              rank='superkingdom', fraction=0.2, f_weighted_at_rank=0.4, 
                              lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None),), 
                              bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2),
                              total_weighted_hashes=0),
           SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                              rank='superkingdom', fraction=0.8, f_weighted_at_rank=0.6,
                              lineage=(), bp_match_at_rank=60, query_ani_at_rank=None,
                              total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                   rank='phylum', fraction=0.2, f_weighted_at_rank=0.4, 
                                   lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                            LineageTuple(rank='phylum', name='b', taxid=None)),
                                   bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2),
                                   total_weighted_hashes=0),
            SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                   rank='phylum', fraction=0.8, f_weighted_at_rank=0.6,
                                   lineage=(), bp_match_at_rank=60, query_ani_at_rank=None,
                                   total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    # best_only picks first if equally good
    cl = [SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                  rank='class', fraction=0.1, f_weighted_at_rank=0.2, 
                                  lineage=(LineageTuple(rank='superkingdom', name='a', taxid=None), 
                                           LineageTuple(rank='phylum', name='b', taxid=None), 
                                           LineageTuple(rank='class', name='c', taxid=None)), 
                                   bp_match_at_rank=20, query_ani_at_rank=approx(0.93, rel=1e-2),
                                   total_weighted_hashes=0), 
           SummarizedGatherResult(query_name='q1', query_md5='md5', query_filename='query_fn',
                                  rank='class', fraction=0.9, f_weighted_at_rank=0.8,
                                  lineage=(), bp_match_at_rank=80, query_ani_at_rank=None, total_weighted_hashes=0)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.4)
    assert q_res.total_f_classified['phylum'] == approx(0.2)
    assert q_res.total_f_classified['class'] == approx(0.1)
    assert q_res.total_bp_classified['superkingdom'] == 40
    assert q_res.total_bp_classified['class'] == 20
    assert q_res.best_only_result == True
