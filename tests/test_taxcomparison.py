"""
Tests for the 'TaxComparison' classes.
"""

#import numpy as np
import pytest
#import sourmash_tst_utils as utils

from sourmash.tax.taxcomparison import LineagePair, LineageTuple, BaseLineageInfo, RankLineageInfo, LINSLineageInfo, build_tree

# sigh, can't make build tree work as easily with both LineagePair and LineageInfo.
# What if LineagePair just had the extra info?
from sourmash.tax.taxcomparison import LineageTuple as LineagePair
from sourmash.lca.lca_utils import find_lca

standard_taxranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
strain_taxranks = standard_taxranks + ['strain']

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
    assert taxinf.taxlist() == standard_taxranks
    assert taxinf.ascending_taxlist() == standard_taxranks[::-1] 


def test_RankLineageInfo_taxlist_with_strain():
    taxinf = RankLineageInfo(include_strain = True)
    assert taxinf.taxlist() == strain_taxranks
    assert taxinf.ascending_taxlist() == strain_taxranks[::-1]


def test_RankLineageInfo_init_lineage_str_1():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_str_1_truncate():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', 'c']


def test_RankLineageInfo_init_lineage_str_2():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', 'c' '', '', '', '', '']


def test_RankLineageInfo_init_lineage_str_2_truncate():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x, include_strain=True)
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
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_strain():
    x = {'superkingdom': 'name1', 'class': 'name2'}
    taxinf = RankLineageInfo(lineage_dict=x, include_strain=True)
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
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '']
    assert taxinf.zip_taxid()== ['1', '', '2', '', '', '', '']


def test_zip_lineage_1():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x, include_strain=True)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage() == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x, include_strain=True)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage(truncate_empty=True))
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x, include_strain=True)
    assert taxinf.zip_lineage() == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_3_truncate():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x, include_strain=True)
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x, include_strain=True)
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


#def test_LINSLineageInfo_taxlist():
#    taxinf = LINSLineageInfo(num_positions=10)
#    assert taxinf.taxlist() == [0]*10
#    assert taxinf.ascending_taxlist() == standard_taxranks[::-1] 