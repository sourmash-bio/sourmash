#! /usr/bin/env python
"""
...

TODO:
* check if we've already seen this md5sum?
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict, Counter
import itertools
import pprint
try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest
import gzip

import sourmash_lib
from . import lca_utils
from ..logging import notify, error
from .. import sourmash_args
from .lca_utils import set_debug, debug

DEFAULT_THRESHOLD=5                  # how many counts of a taxid at min


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
        tree = lca_utils.build_tree(tuple_info)

        # also update a tree that we can ascend from leaves -> parents
        # for all assignments for all hashvals
        parents = lca_utils.build_reverse_tree(tuple_info, parents)

        # now find either a leaf or the first node with multiple
        # children; that's our least-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)
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
        lca_utils.build_tree([xx], tree)

    if n > 1:
        debug('XXX', n)

    status = 'nomatch'
    if not tree:
        return [('', '')], status

    # now find LCA? or whatever.
    lca, reason = lca_utils.find_lca(tree)
    if reason == 0:               # leaf node
        status = 'found'
        debug('END', lca)
    else:                         # internal node
        status = 'disagree'
        debug('MULTI', lca)

    # backtrack to full lineage via parents
    lineage = []
    parent = lca
    while parent != ('root', 'root'):
        lineage.insert(0, parent)
        parent = parents.get(parent)

    debug(parents)
    debug('lineage is:', lineage)

    return lineage, status


def classify(args):
    p = argparse.ArgumentParser()
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV to this file instead of stdout')
    p.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath directories.')
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
        set_debug(args.debug)

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

        lca_db = lca_utils.LCA_Database()
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

    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.query))
    else:
        inp_files = list(args.query)

    # for each query, gather all the matches across databases
    csvfp = csv.writer(sys.stdout)
    if args.output:
        notify("outputting classifications to '{}'", args.output.name)
        csvfp = csv.writer(args.output)
    else:
        notify("outputting classifications to stdout")
    csvfp.writerow(['ID','status'] + lca_utils.taxlist)

    total_count = 0
    n = 0
    total_n = len(inp_files)
    for query_filename in inp_files:
        n += 1
        for query_sig in sourmash_lib.load_signatures(query_filename,
                                                      ksize=ksize):
            notify(u'\r\033[K', end=u'')
            notify('... classifying {} (file {} of {})', query_sig.name(), n, total_n, end='\r')
            debug('classifying', query_sig.name())
            total_count += 1

            # make sure we're looking at the same scaled value as database
            query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

            lineage, status = classify_signature(query_sig, dblist, args.threshold)

            # output!
            row = [query_sig.name(), status]
            for taxrank, (rank, name) in zip_longest(lca_utils.taxlist,
                                                     lineage,
                                                     fillvalue=('', '')):
                if rank:
                    assert taxrank == rank
                row.append(name)

            if not args.output:           # a bit of a hack to clear line
                notify(u'\r\033[K', end=u'')
            csvfp.writerow(row)

    notify(u'\r\033[K', end=u'')
    notify('classified {} signatures total', total_count)


if __name__ == '__main__':
    sys.exit(classify(sys.argv[1:]))
