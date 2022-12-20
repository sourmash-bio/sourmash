"""
Tests for the 'TaxComparison' classes.
"""

#import numpy as np
import pytest
#import sourmash_tst_utils as utils

from sourmash.tax.taxcomparison import LineagePair, LineageTuple, LineageInfoRanks, LineageInfoLINS


standard_taxranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
strain_taxranks = standard_taxranks + ['strain']

def test_LineageInfoRanks_taxlist():
    taxinf = LineageInfoRanks()
    assert taxinf.taxlist() == standard_taxranks
    assert taxinf.ascending_taxlist() == standard_taxranks[::-1] 


def test_LineageInfoRanks_taxlist_with_strain():
    taxinf = LineageInfoRanks(include_strain = True)
    assert taxinf.taxlist() == strain_taxranks
    assert taxinf.ascending_taxlist() == strain_taxranks[::-1]


def test_LineageInfoRanks_init_lineage_str():
    x = "a;b;c"
    taxinf = LineageInfoRanks(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c', '', '', '', '', '']

def test_LineageInfoRanks_init_lineage_str_truncate():
    x = "a;b;c"
    taxinf = LineageInfoRanks(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', 'c']

def test_LineageInfoRanks_init_lineage_str_2():
    x = "a;b;;c"
    taxinf = LineageInfoRanks(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', 'c' '', '', '', '', '']

def test_LineageInfoRanks_init_lineage_str_2_truncate():
    x = "a;b;;c"
    taxinf = LineageInfoRanks(lineage_str=x, include_strain=True)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', '', 'c']


def test_LineageInfoRanks_init_lineage_with_incorrect_rank():
    x = [ LineagePair('superkingdom', 'a'), LineagePair("NotARank", ''), LineagePair('class', 'c') ]
    with pytest.raises(ValueError) as exc:
        LineageInfoRanks(lineage=x)
    print(str(exc))
    assert f"Rank 'NotARank' not present in " in str(exc)


def test_zip_lineage_1():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = LineageInfoRanks(lineage=x, include_strain=True) 
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage() == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineageTuple('superkingdom', 'a'), LineageTuple('phylum', 'b') ]
    taxinf = LineageInfoRanks(lineage=x, include_strain=True) 
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage(truncate_empty=True))
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = LineageInfoRanks(lineage=x, include_strain=True) 
    assert taxinf.zip_lineage() == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_3_truncate():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = LineageInfoRanks(lineage=x, include_strain=True) 
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    taxinf = LineageInfoRanks(lineage=x, include_strain=True) 
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']

#    assert 'incomplete lineage at phylum - is class instead' in str(e.value)


def test_display_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    taxinf = LineageInfoRanks(lineage=x) 
    assert taxinf.display_lineage() == "a;b"


def test_display_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = LineageInfoRanks(lineage=x) 
    assert taxinf.display_lineage() == "a;;c"