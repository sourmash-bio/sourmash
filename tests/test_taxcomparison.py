"""
Tests for the 'TaxComparison' classes.
"""

#import numpy as np
import pytest

from sourmash.tax.taxcomparison import LineagePair, LineageTuple, LineageInfoRanks, LineageInfoLINS

#import sourmash_tst_utils as utils

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

    #with pytest.raises(TypeError) as exc:
    #    cmp.angular_similarity
    #print(str(exc))
    #assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)

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

