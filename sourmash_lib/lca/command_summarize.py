#! /usr/bin/env python
"""
Summarize the taxonomic content of the given signatures.
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


DEFAULT_THRESHOLD=5


def summarize(hashvals, dblist, threshold):
    """
    Classify 'hashvals' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, counts) where 'lineage' is a lineage tuple
    [(rank, name), ...].
    """

    # gather assignments from across all the databases
    these_assignments = defaultdict(list)
    n_custom = 0
    for hashval in hashvals:
        for lca_db in dblist:
            assignments = lca_db.hashval_to_lineage_id.get(hashval, [])
            for lineage_id in assignments:
                assignment = lca_db.lineage_dict[lineage_id]
                these_assignments[hashval].append(assignment)
                n_custom += 1

    # now convert to trees -> do LCA & counts
    counts = Counter()
    parents = {}
    for hashval in these_assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        tuple_info = these_assignments[hashval]
        tree = lca_utils.build_tree(tuple_info)

        # also update a tree that we can ascend from leaves -> parents
        # for all assignments for all hashvals
        parents = lca_utils.build_reverse_tree(tuple_info, parents)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)
        counts[lca] += 1

    # ok, we now have the LCAs for each hashval, and their number
    # of counts. Now sum across "significant" LCAs - those above
    # threshold - and read out results.

    tree = {}
    tree_counts = defaultdict(int)

    debug(pprint.pformat(counts.most_common()))

    n = 0
    for lca, count in counts.most_common():
        if count < threshold:
            break

        n += 1

        # construct [(rank, name), ...] lineage
        lineage = []
        parent = lca
        debug('lca is', lca)
        while parent != ('root', 'root'):
            debug('inserting', parent)
            lineage.insert(0, parent)
            tree_counts[parent] += count
            parent = parents.get(parent)

        debug(lineage)

        # update tree with this set of assignments
        lca_utils.build_tree([lineage], tree)

    def traverse_and_count(tree, count_d, lineage_so_far=()):
        """
        Convert tree into list of lineages, recursively.
        """
        if not lineage_so_far or lineage_so_far == (('root', 'root')):
            for (rank, name), v in tree.items():
                count_d[((rank, name),)] += tree_counts[(rank, name)]
                lineage_tup = ( (rank, name,), )
                traverse_and_count(v, count_d, lineage_tup)
            return

        for (rank, name), v in tree.items():
            sub_lineage = list(lineage_so_far)
            sub_lineage.append((rank, name))
            sub_lineage = tuple(sub_lineage)
            count_d[sub_lineage] += tree_counts[(rank, name)]
            traverse_and_count(v, count_d, sub_lineage)

    lineage_counts = defaultdict(int)
    traverse_and_count(tree, lineage_counts)
    debug(pprint.pformat(lineage_counts))
    return lineage_counts


def summarize_main(args):
    """
    main summarization function: 
    """
    p = argparse.ArgumentParser()
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath directories.')
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
    notify('loaded {} databases for LCA use.', len(dblist))

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()
    notify('ksize={} scaled={}', ksize, scaled)

    # find all the queries
    notify('finding query signatures...')
    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.query))
    else:
        inp_files = list(args.query)

    # for each query, gather all the hashvals across databases
    total_count = 0
    n = 0
    total_n = len(inp_files)
    hashvals = defaultdict(int)
    for query_filename in inp_files:
        n += 1
        for query_sig in sourmash_lib.load_signatures(query_filename,
                                                      ksize=ksize):
            notify(u'\r\033[K', end=u'')
            notify('... loading {} (file {} of {})', query_sig.name(), n, total_n, end='\r')

            mh = query_sig.minhash.downsample_scaled(scaled)
            for hashval in mh.get_mins():
                hashvals[hashval] += 1

    notify(u'\r\033[K', end=u'')
    notify('loaded signatures from {} files total.', n)

    lineage_counts = summarize(hashvals, dblist, args.threshold)

    for (lineage_tup, count) in lineage_counts.items():
        lineage = ';'.join([ name for (rank, name) in lineage_tup ])
        print(count, lineage)
