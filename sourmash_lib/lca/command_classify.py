#! /usr/bin/env python
"""
...

TODO:
* check if we've already seen this md5sum?
"""
import sys
import argparse
import csv
from collections import defaultdict, Counter
import itertools
import pprint
import json

import sourmash_lib
from ..logging import notify, error

DEFAULT_THRESHOLD=5                  # how many counts of a taxid at min

taxlist = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
           'species']

_print_debug = False
def debug(*args):
    if _print_debug:
        print(*args)


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of (rank, name) tuples
    in 'assignments'.  This tree can then be used to find least common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    for assignment in assignments:
        node = tree

        for rank, name in assignment:
            if name:
                child = node.get((rank, name), {})
                node[(rank, name)] = child

                # shift -> down in tree
                node = child

    return tree


def test_build_tree():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def test_build_tree_3():                  # empty 'rank2' name
    tree = build_tree([[('rank1', 'name1'), ('rank2', '')]])
    assert tree == { ('rank1', 'name1'): {} }


def test_build_tree_4():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                      ])

    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                      ], tree)

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return ((rank, name), reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    cur = ('root', 'root')
    while 1:
        if len(node) == 1:                # descend to only child
            cur = next(iter(node.keys()))
            node = node[cur]
        elif len(node) == 0:              # at leaf; end
            return cur, 0
        else:                             # len(node) > 1 => confusion!!
            return cur, len(node)


def test_find_lca():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == (('rank2', 'name2'), 0)


def test_find_lca_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == (('rank1', 'name1'), 2)


def build_reverse_tree(assignments, initial=None):
    """
    Builds a child -> parent dictionary (a reverse DAG) from lists of
    (rank, name) tuples in 'assignments'.
    """
    if initial is None:
        parents = {}
    else:
        parents = initial

    for assignment in assignments:
        last_node = ('root', 'root')
        for rank, name in assignment:
            if name:
                parents[(rank, name)] = last_node
                last_node = (rank, name)

    return parents


def test_build_reverse_tree():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2')]])

    print(parents)
    assert parents == { ('rank2', 'name2'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_2():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 [('rank1', 'name1'), ('rank2', 'name2b')],
                                 ])

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_3():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 ])
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                                 ], parents)

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_4():          # empty 'rank2' name
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', '')]])

    print(parents)
    assert parents == { ('rank1', 'name1'): ('root', 'root') }


class LCA_Database(object):
    def __init__(self):
        self.lineage_dict = None
        self.hashval_to_lineage_id = None
        self.ksize = None
        self.scaled = None
        self.signatures_to_lineage = None

    def load(self, db_name):
        with open(db_name, 'rt') as fp:
            load_d = json.load(fp)
            version = load_d['version']
            assert version == '1.0'

            type = load_d['type']
            assert type == 'sourmash_lca'

            ksize = load_d['ksize']
            scaled = load_d['scaled']

            lineage_dict_2 = load_d['lineages']
            lineage_dict = {}
            for k, v in lineage_dict_2.items():
                vv = []
                for rank in taxlist:
                    name = v.get(rank, '')
                    vv.append((rank, name))

                lineage_dict[int(k)] = tuple(vv)

            hashval_to_lineage_id_2 = load_d['hashval_assignments']
            hashval_to_lineage_id = {}
            for k, v in hashval_to_lineage_id_2.items():
                hashval_to_lineage_id[int(k)] = v

            signatures_to_lineage = load_d['signatures_to_lineage']

        self.lineage_dict = lineage_dict
        self.hashval_to_lineage_id = hashval_to_lineage_id
        self.ksize = ksize
        self.scaled = scaled
        self.signatures_to_lineage = signatures_to_lineage


def classify_signature(query_sig, dblist, threshold):
    # gather assignments from across all the databases
    these_assignments = defaultdict(list)
    n_custom = 0
    for hashval in query_sig.minhash.get_mins():
        for lca_db in dblist:
            assignments = lca_db.hashval_to_lineage_id.get(hashval, [])
            for lineage_id in assignments:
                assignment = lca_db.lineage_dict[lineage_id]
                these_assignments[hashval].append(assignment)
                n_custom += 1

    # count number of assignments for each most-specific
    check_counts = Counter()
    for tuple_info in these_assignments.values():
        last_tup = tuple(tuple_info[-1])
        check_counts[last_tup] += 1

    debug('n custom hashvals:', n_custom)
    debug(pprint.pformat(check_counts.most_common()))

    # now convert to trees -> do LCA & counts
    counts = Counter()
    parents = {}
    for hashval in these_assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover least-common-ancestor.
        tuple_info = these_assignments[hashval]
        tree = build_tree(tuple_info)

        # also update a tree that we can ascend from leaves -> parents
        # for all assignments for all hashvals
        parents = build_reverse_tree(tuple_info, parents)

        # now find either a leaf or the first node with multiple
        # children; that's our least-common-ancestor node.
        lca, reason = find_lca(tree)
        counts[lca] += 1

    # ok, we now have the LCAs for each hashval, and their number
    # of counts. Now sum across "significant" LCAs - those above
    # threshold.

    tree = {}
    tree_counts = defaultdict(int)

    debug(pprint.pformat(counts.most_common()))

    n = 0
    for lca, count in counts.most_common():
        if count < threshold:
            break

        n += 1

        xx = []
        parent = lca
        while parent:
            xx.insert(0, parent)
            tree_counts[parent] += count
            parent = parents.get(parent)
        debug(n, count, xx[1:])

        # update tree with this set of assignments
        build_tree([xx], tree)

    if n > 1:
        debug('XXX', n)

    # now find LCA? or whatever.
    lca, reason = find_lca(tree)
    if reason == 0:               # leaf node
        debug('END', lca)
    else:                         # internal node
        debug('MULTI', lca)

    # backtrack to full lineage via parents
    lineage = []
    parent = lca
    while parent != ('root', 'root'):
        lineage.insert(0, parent)
        parent = parents.get(parent)

    debug(parents)
    debug('lineage is:', lineage)

    return lineage


def classify(args):
    p = argparse.ArgumentParser()
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV to this file instead of stdout')
    #p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    if not args.query:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    if args.debug:
        global _print_debug
        _print_debug = True

    ksize_vals = set()
    scaled_vals = set()
    dblist = []

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]
    
    for db_name in args.db:
        notify(u'\r\033[K', end=u'', file=sys.stderr)
        notify('... loading database {}'.format(db_name), end='\r',
              file=sys.stderr)

        lca_db = LCA_Database()
        lca_db.load(db_name)

        ksize_vals.add(lca_db.ksize)
        if len(ksize_vals) > 1:
            raise Exception('multiple ksizes, quitting')
        scaled_vals.add(lca_db.scaled)
        if len(scaled_vals) > 1:
            raise Exception('multiple scaled vals, quitting')

        dblist.append(lca_db)

    notify(u'\r\033[K', end=u'')
    notify('loaded {} databases for LCA use.', len(dblist))

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()
    notify('ksize={} scaled={}', ksize, scaled)
        
    # for each query, gather all the matches across databases, then
    csvfp = csv.writer(sys.stdout)
    if args.output:
        notify("outputting classifications to '{}'", args.output.name)
        csvfp = csv.writer(args.output)
    else:
        notify("outputting classifications to stdout")
    csvfp.writerow(['ID'] + taxlist)

    total_count = 0
    n = 0
    total_n = len(args.query)
    for query_filename in args.query:
        n += 1
        for query_sig in sourmash_lib.load_signatures(query_filename,
                                                      ksize=ksize):
            notify(u'\r\033[K', end=u'')
            notify('... classifying {} (file {} of {})', query_sig.name(), n, total_n, end='\r')
            debug('classifying', query_sig.name())
            total_count += 1

            # make sure we're looking at the same scaled value as database
            query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

            lineage = classify_signature(query_sig, dblist, args.threshold)

            # output!
            row = [query_sig.name()]
            for taxrank, (rank, name) in itertools.zip_longest(taxlist, lineage, fillvalue=('', '')):
                if rank:
                    assert taxrank == rank
                row.append(name)

            csvfp.writerow(row)

    notify(u'\r\033[K', end=u'')
    notify('classified {} signatures total', total_count)


if __name__ == '__main__':
    sys.exit(classify(sys.argv[1:]))
