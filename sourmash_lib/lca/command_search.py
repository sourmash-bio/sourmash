#! /usr/bin/env python
"""
Summarize the taxonomic content of the given signatures, combined.
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict, Counter

import sourmash_lib
from sourmash_lib import sourmash_args
from sourmash_lib.logging import notify, error, print_results
from sourmash_lib.lca import lca_utils
from sourmash_lib.lca.lca_utils import debug, set_debug


DEFAULT_THRESHOLD=5


def search(hashvals, dblist):
    """
    """

    # gather assignments from across all the databases
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            assignments[hashval].update(lineages)

    # now convert to trees -> do LCA & counts
    counts = Counter()
    for hashval in assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        tuple_info = assignments[hashval]
        tree = lca_utils.build_tree(tuple_info)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)
        counts[lca] += 1

    # ok, we now have the LCAs for each hashval, and their number
    # of counts. Now aggregate counts across the tree, going up from
    # the leaves.
    tree = {}

    debug(counts.most_common())

    aggregated_counts = defaultdict(int)
    for lca, count in counts.most_common():
        if count < threshold:
            break

        # climb from the lca to the root.
        while lca:
            aggregated_counts[lca] += count
            lca = lca[:-1]

    debug(aggregated_counts)

    return aggregated_counts


def lca_search_main(args):
    """
    main lca search function.
    """
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('dblist', nargs='+')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='CSV output')
    p.add_argument('--scaled', type=float)
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-k', '--ksize', type=int, default=31)
    args = p.parse_args(args)

    if args.debug:
        set_debug(args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load query
    query_sig = sourmash_lib.load_one_signature(args.query, ksize=args.ksize)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.dblist, args.scaled)
    notify('ksize={} scaled={}', ksize, scaled)
    assert ksize == args.ksize

    # get the correct query downsampling for the database
    query_mh = query_sig.minhash
    query_mh = query_mh.downsample_scaled(scaled)

    # for each hash in query, gather all matches
    counts = Counter()
    for hashval in query_mh.get_mins():
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            if lineages:
                lineages = tuple(sorted(lineages))
                counts[lineages] += 1

    total = len(query_mh.get_mins())

    debug(counts)

    for lineages, count in counts.most_common():
        if count < args.threshold:
            break

        debug('ZZZ', lineages)

        if not lineages:
            display = ['(root)']
        else:
            display = []
            for lineage in lineages:
                debug('ZZZZ', lineage)
                if lineage:
                    lineage = lca_utils.zip_lineage(lineage,
                                                    truncate_empty=True)
                    lineage = ';'.join(lineage)
                else:
                    lineage = '(root)'
                    assert 0
                display.append(lineage)

        debug(display)

        p = count / total * 100.
        p = '{:.1f}%'.format(p)

        print_results('{:5} {:>5}   {}', p, count, display[0])
        for d in display[1:]:
            print_results('              {}', d)


if __name__ == '__main__':
    sys.exit(lca_search_main(sys.argv[1:]))
