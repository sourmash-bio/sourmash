"""
Utility functions for lowest-common-ancestor analysis tools.
"""
from __future__ import print_function, division
from os.path import exists
from collections import namedtuple, defaultdict, Counter

from .lca_db import LCA_Database, load_single_database, load_databases


__all__ = ['taxlist', 'zip_lineage', 'build_tree', 'find_lca',
           'load_single_database', 'load_databases', 'gather_assignments',
           'count_lca_for_assignments', 'LineagePair', 'display_lineage',
           'count_lca_for_assignments_abund']

try:                                      # py2/py3 compat
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

from sourmash.logging import notify, error, debug

# type to store an element in a taxonomic lineage
LineagePair = namedtuple('LineagePair', ['rank', 'name'])


def check_files_exist(*files):
    ret = True
    not_found = []
    for f in files:
        if not exists(f):
            not_found.append(f)
            ret = False

    if len(not_found):
        error('Error! Could not find the following files.'
              ' Make sure the file paths are specified correctly.\n{}'.format('\n'.join(not_found)))

    return ret


# ordered list of taxonomic ranks
def taxlist(include_strain=True):
    """
    Provide an ordered list of taxonomic ranks.
    """
    for k in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'species']:
        yield k
    if include_strain:
        yield 'strain'


# produce an ordered list of tax names from lineage
def zip_lineage(lineage, include_strain=True, truncate_empty=False):
    """
    Given an iterable of LineagePair objects, return list of lineage names.

    This utility function handles species/strain and empty lineage entries
    gracefully.

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    >>> zip_lineage(x)
    ['a', 'b', '', '', '', '', '', '']

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    >>> zip_lineage(x)
    ['a', '', 'c', '', '', '', '', '']
    """

    empty = LineagePair(None, '')

    pairs = zip_longest(taxlist(include_strain=include_strain),
                        lineage, fillvalue=empty)
    pairs = list(pairs)

    # eliminate empty if so requested
    if truncate_empty:
        last_lineage_tup = pairs[-1][1]
        while pairs and last_lineage_tup == empty:
            pairs.pop(-1)
            if pairs:
                last_lineage_tup = pairs[-1][1]

    row = []
    for taxrank, lineage_tup in pairs:
        # validate non-empty tax, e.g. superkingdom/phylum/class in order.
        if lineage_tup != empty and lineage_tup.rank != taxrank:
            raise ValueError('incomplete lineage at {} - is {} instead'.format(taxrank, lineage_tup.rank))

        row.append(lineage_tup.name)
    return row


def display_lineage(lineage, include_strain=True, truncate_empty=True):
    return ";".join(zip_lineage(lineage,
                                include_strain=include_strain,
                                truncate_empty=truncate_empty))


# filter function toreplace blank/na/null with 'unassigned'
filter_null = lambda x: 'unassigned' if x.strip() in \
  ('[Blank]', 'na', 'null', '') else x
null_names = set(['[Blank]', 'na', 'null'])


def make_lineage(lineage_str):
    "Turn a ; or ,-separated set of lineages into a tuple of LineagePair objs."
    lin = lineage_str.split(';')
    if len(lin) == 1:
        lin = lineage.split(',')
    lin = [ LineagePair(rank, n) for (rank, n) in zip(taxlist(), lin) ]
    lin = tuple(lin)

    return lin


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of LineagePair objects
    in 'assignments'.  This tree can then be used to find lowest common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    if not assignments:
        raise ValueError("empty assignment passed to build_tree")

    for assignment in assignments:
        node = tree

        for lineage_tup in assignment:
            if lineage_tup.name:
                child = node.get(lineage_tup, {})
                node[lineage_tup] = child

                # shift -> down in tree
                node = child

    return tree


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return (lineage_tup, reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    lineage = []
    while 1:
        if len(node) == 1:                # descend to only child; track path
            lineage_tup = next(iter(node.keys()))
            lineage.append(lineage_tup)
            node = node[lineage_tup]
        elif len(node) == 0:              # at leaf; end
            return tuple(lineage), 0
        else:                             # len(node) > 1 => confusion!!
            return tuple(lineage), len(node)


def gather_assignments(hashvals, dblist):
    """
    Gather assignments from across all the databases for all the hashvals.

    Ignores counts of the hashvals.
    """
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            if lineages:
                assignments[hashval].update(lineages)

    return assignments


def count_lca_for_assignments(assignments):
    """
    For each hashval, count the LCA across its assignments.

    Ignores counts for each hashval
    """
    counts = Counter()
    for hashval in assignments:
        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = find_lca(tree)
        counts[lca] += 1

    return counts


def count_lca_for_assignments_abund(assignments, hashval_counts):
    """
    For each hashval, count the LCA across its assignments.

    Use hashval_counts to weight lineages.
    """
    counts = Counter()
    for hashval in assignments:
        assert hashval in hashval_counts

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = find_lca(tree)
        counts[lca] += hashval_counts[hashval]

    return counts
