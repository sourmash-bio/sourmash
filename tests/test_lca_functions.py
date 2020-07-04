"""
Tests for functions in lca submodule.
"""
import pytest

from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, build_tree, find_lca,
                                    taxlist, count_lca_for_assignments,
                                    zip_lineage, display_lineage,
                                    make_lineage, is_lineage_match,
                                    pop_to_rank)


class FakeLCA_Database(object):
    def __init__(self):
        self._assignments = {}

    def _set_lineage_assignment(self, hashval, assignment):
        self._assignments[hashval] = assignment

    def get_lineage_assignments(self, hashval):
        if hashval in self._assignments:
            return self._assignments[hashval]
        else:
            return None


def test_taxlist_1():
    assert list(taxlist()) == ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']


def test_taxlist_2():
    assert list(taxlist(include_strain=False)) == ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def test_zip_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    assert zip_lineage(x) == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    assert zip_lineage(x, truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    assert zip_lineage(x) == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_3_truncate():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    assert zip_lineage(x, truncate_empty=True) == ['a', '', 'c']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    with pytest.raises(ValueError) as e:
        zip_lineage(x)

    assert 'incomplete lineage at phylum - is class instead' in str(e.value)


def test_display_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    assert display_lineage(x) == "a;b", display_lineage(x)


def test_display_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    assert display_lineage(x) == "a;;c", display_lineage(x)


def test_build_tree():
    tree = build_tree([[LineagePair('rank1', 'name1'),
                        LineagePair('rank2', 'name2')]])
    assert tree == { LineagePair('rank1', 'name1'):
                         { LineagePair('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }


def test_build_tree_3():                  # empty 'rank2' name
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', '')]])
    assert tree == { LineagePair('rank1', 'name1'): {} }


def test_build_tree_4():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                      ])

    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ], tree)

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }

def test_build_tree_5():
    with pytest.raises(ValueError):
        tree = build_tree([])


def test_find_lca():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2'),), 0)


def test_find_lca_2():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'),), 2)


def test_find_lca_3():
    lin1 = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b')

    tree = build_tree([lin1, lin2])
    lca, reason = find_lca(tree)
    assert lca == lin1                    # find most specific leaf node


def test_gather_assignments_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin ])


def test_gather_assignments_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin, lin2 ])


def test_gather_assignments_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval, hashval2], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin, lin2 ])
    assert assignments[hashval2] == set([ lin ])


def test_count_lca_for_assignments_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert len(counts) == 1
    assert counts[lin] == 1


def test_count_lca_for_assignments_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert counts[lin] == 0
    assert counts[lin2] == 0

    assert len(counts) == 1
    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 1


def test_count_lca_for_assignments_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval, hashval2], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 1
    assert counts[lin2] == 0

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 1


def test_count_lca_for_assignments_abund_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    hashval_counts = dict()
    hashval_counts[hashval] = 3

    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts.keys(), [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 1
    assert counts[lin] == 3


def test_count_lca_for_assignments_abund_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    hashval_counts = dict()
    hashval_counts[hashval] = 3

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert counts[lin] == 0
    assert counts[lin2] == 0

    assert len(counts) == 1
    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 3           # yes!


def test_count_lca_for_assignments_abund_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 5               # makes sense
    assert counts[lin2] == 0              # makes sense

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 2           # yes!

def test_count_lca_for_assignments_abund_4():
    # test basic mechanics of gather_assignments function with three lineages
    # and three hashvals
    hashval = 12345678
    hashval2 = 87654321
    hashval3 = 34567891
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5
    hashval_counts[hashval3] = 3

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')
    lin3 = lca_utils.make_lineage('a;b;d;e')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ])) # lca: a;b
    db._set_lineage_assignment(hashval2, set([ lin ])) # lca: a;b;c
    db._set_lineage_assignment(hashval3, set([ lin2, lin3 ])) # a;b;d;e

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 3
    assert counts[lin] == 5               # makes sense b/c hashval2
    assert counts[lin2] == 0              # a;b;d (lin2) + a;b;d;e (lin3) -->a;b;d;e (lin3) only
    assert counts[lin3] == 3              # hashval3

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 2           # yes, b/c hashval

def test_count_lca_for_assignments_abund_5():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals when linages match but one has lower taxo detail
    hashval = 12345678
    hashval2 = 87654321
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5

    lin = lca_utils.make_lineage('a;b;d')
    lin2 = lca_utils.make_lineage('a;b;d;e')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 5               # makes sense
    assert counts[lin2] == 2              # lin+lin2 yield just lin2


def test_is_lineage_match_1():
    # basic behavior: match at order and above, but not at family or below.
    lin1 = make_lineage('d__a;p__b;c__c;o__d;f__e')
    lin2 = make_lineage('d__a;p__b;c__c;o__d;f__f')

    assert is_lineage_match(lin1, lin2, 'superkingdom')
    assert is_lineage_match(lin1, lin2, 'phylum')
    assert is_lineage_match(lin1, lin2, 'class')
    assert is_lineage_match(lin1, lin2, 'order')
    assert not is_lineage_match(lin1, lin2, 'family')
    assert not is_lineage_match(lin1, lin2, 'genus')
    assert not is_lineage_match(lin1, lin2, 'species')


def test_is_lineage_match_2():
    # match at family, and above, levels; no genus or species to match
    lin1 = make_lineage('d__a;p__b;c__c;o__d;f__f')
    lin2 = make_lineage('d__a;p__b;c__c;o__d;f__f')

    assert is_lineage_match(lin1, lin2, 'superkingdom')
    assert is_lineage_match(lin1, lin2, 'phylum')
    assert is_lineage_match(lin1, lin2, 'class')
    assert is_lineage_match(lin1, lin2, 'order')
    assert is_lineage_match(lin1, lin2, 'family')
    assert not is_lineage_match(lin1, lin2, 'genus')
    assert not is_lineage_match(lin1, lin2, 'species')


def test_is_lineage_match_3():
    # one lineage is empty
    lin1 = make_lineage('')
    lin2 = make_lineage('d__a;p__b;c__c;o__d;f__f')

    assert not is_lineage_match(lin1, lin2, 'superkingdom')
    assert not is_lineage_match(lin1, lin2, 'family')
    assert not is_lineage_match(lin1, lin2, 'order')
    assert not is_lineage_match(lin1, lin2, 'class')
    assert not is_lineage_match(lin1, lin2, 'phylum')
    assert not is_lineage_match(lin1, lin2, 'genus')
    assert not is_lineage_match(lin1, lin2, 'species')


def test_pop_to_rank_1():
    # basic behavior - pop to order?
    lin1 = make_lineage('d__a;p__b;c__c;o__d')
    lin2 = make_lineage('d__a;p__b;c__c;o__d;f__f')

    print(lin1)
    print(pop_to_rank(lin2, 'order'))
    assert pop_to_rank(lin2, 'order') == lin1


def test_pop_to_rank_2():
    # what if we're already above rank?
    lin2 = make_lineage('d__a;p__b;c__c;o__d;f__f')

    print(pop_to_rank(lin2, 'species'))
    assert pop_to_rank(lin2, 'species') == lin2
