#! /usr/bin/env python
"""
Summarize rank-specific information from LCAs in one or more databases.
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict, OrderedDict, Counter
import json
import pprint

import sourmash_lib
from sourmash_lib import sourmash_args
from sourmash_lib.logging import notify, error
from sourmash_lib.lca import lca_utils
from sourmash_lib.lca.lca_utils import debug, set_debug


def make_lca_counts(dblist):
    """
    """

    # gather all hashvalue assignments from across all the databases
    assignments = defaultdict(set)
    for lca_db in dblist:
        for hashval, lid_list in lca_db.hashval_to_lineage_id.items():
            lineages = [lca_db.lineage_dict[lid] for lid in lid_list]
            assignments[hashval].update(lineages)

    # now convert to trees -> do LCA & counts
    counts = defaultdict(int)
    for hashval, lineages in assignments.items():

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        debug(pprint.pformat(lineages))
        tree = lca_utils.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)
        counts[lca] += 1

    return counts


def rankinfo_main(args):
    """
    main summarization function:
    """
    p = argparse.ArgumentParser()
    p.add_argument('db', nargs='+')
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    if args.debug:
        set_debug(args.debug)

    ksize_vals = set()
    scaled_vals = set()
    dblist = []

    # load all the databases
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
    notify('loaded {} databases.', len(dblist))

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()
    notify('ksize={} scaled={}', ksize, scaled)

    counts = make_lca_counts(dblist)

    counts_by_rank = defaultdict(int)
    for lineage, count in counts.items():
        if lineage:
            lineage_tup = lineage[-1]
            counts_by_rank[lineage_tup.rank] += count

#        p = '{:.1f}%'.format(p)

#        print('{:5} {:>5}   {}'.format(p, count, lineage))

    for rank in lca_utils.taxlist():
        count = counts_by_rank.get(rank, 0)
        print('{},{}'.format(rank, count))


if __name__ == '__main__':
    sys.exit(rankinfo_main(sys.argv[1:]))
