#! /usr/bin/env python
"""
Classify individual signature files down to deepest possible node.
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
from sourmash_lib import sourmash_args
from sourmash_lib.logging import notify, error
from sourmash_lib.lca import lca_utils
from sourmash_lib.lca.lca_utils import debug, set_debug, LineagePair

DEFAULT_THRESHOLD=5                  # how many counts of a taxid at min


def classify_signature(query_sig, dblist, threshold):
    """
    Classify 'query_sig' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, status) where 'lineage' is a lineage tuple
    [(rank, name), ...] and 'status' is either 'nomatch', 'found',
    or 'disagree'.

    This function proceeds in two stages:

       * first, build a list of assignments for all the lineages for each
         hashval.  (For e.g. kraken, this is done in the database preparation
         step; here, we do it dynamically each time.
       * then, across all the hashvals, count the number of times each linage
         shows up, and filter out low-abundance ones (under threshold).
         Then, determine the LCA of all of those.

      """
    # gather assignments from across all the databases
    assignments = defaultdict(list)
    for hashval in query_sig.minhash.get_mins():
        for lca_db in dblist:
            lineage = lca_db.get_lineage_assignments(hashval)
            assignments[hashval].extend(lineage)

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
    # of counts. Now sum across "significant" LCAs - those above
    # threshold - and build a tree.

    tree = {}
    tree_counts = defaultdict(int)

    debug(pprint.pformat(counts.most_common()))

    n = 0
    for lca, count in counts.most_common():
        if count < threshold:
            break

        n += 1

        # update tree with this set of assignments
        lca_utils.build_tree([lca], tree)

    if n > 1:
        debug('XXX', n)

    status = 'nomatch'
    if not tree:
        return [('', '')], status

    # now find lowest-common-ancestor of the resulting tree.
    lca, reason = lca_utils.find_lca(tree)
    if reason == 0:               # leaf node
        status = 'found'
        debug('END', lca)
    else:                         # internal node => disagreement
        status = 'disagree'
        debug('MULTI', lca)

    debug('lineage is:', lca)

    return lca, status


def classify(args):
    """
    main single-genome classification function: 
    """
    p = argparse.ArgumentParser()
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV to this file instead of stdout')
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

    # set up output
    csvfp = csv.writer(sys.stdout)
    if args.output:
        notify("outputting classifications to '{}'", args.output.name)
        csvfp = csv.writer(args.output)
    else:
        notify("outputting classifications to stdout")
    csvfp.writerow(['ID','status'] + list(lca_utils.taxlist()))

    # for each query, gather all the matches across databases
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

            # do the classification
            lineage, status = classify_signature(query_sig, dblist,
                                                 args.threshold)

            # output each classification to the spreadsheet
            row = [query_sig.name(), status]
            debug(lineage)
            for taxrank, (rank, name) in zip_longest(lca_utils.taxlist(),
                                                     lineage,
                                                     fillvalue=(None, '')):
                if rank is not None and name:
                    assert taxrank == rank, (taxrank, rank)
                row.append(name)

            # when outputting to stdout, make output intelligible
            if not args.output:
                notify(u'\r\033[K', end=u'')
            csvfp.writerow(row)

    notify(u'\r\033[K', end=u'')
    notify('classified {} signatures total', total_count)


if __name__ == '__main__':
    sys.exit(classify(sys.argv[1:]))
