#! /usr/bin/env python
"""
Classify individual signature files down to deepest possible node.
"""
from __future__ import print_function
import sys
import argparse
import csv

from .. import sourmash_args, load_signatures
from ..logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import check_files_exist
from ..sourmash_args import SourmashArgumentParser

DEFAULT_THRESHOLD=5                  # how many counts of a taxid at min


def classify_signature(query_sig, dblist, threshold):
    """
    Classify 'query_sig' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, status) where 'lineage' is a tuple of LineagePairs
    and 'status' is either 'nomatch', 'found', or 'disagree'.

    This function proceeds in two stages:

       * first, build a list of assignments for all the lineages for each
         hashval.  (For e.g. kraken, this is done in the database preparation
         step; here, we do it dynamically each time.
       * then, across all the hashvals, count the number of times each linage
         shows up, and filter out low-abundance ones (under threshold).
         Then, determine the LCA of all of those.

      """
    # gather assignments from across all the databases
    assignments = lca_utils.gather_assignments(query_sig.minhash.get_mins(),
                                               dblist)

    # now convert to trees -> do LCA & counts
    counts = lca_utils.count_lca_for_assignments(assignments)
    debug(counts.most_common())

    # ok, we now have the LCAs for each hashval, and their number of
    # counts. Now build a tree across "significant" LCAs - those above
    # threshold.

    tree = {}

    for lca, count in counts.most_common():
        if count < threshold:
            break

        # update tree with this set of assignments
        lca_utils.build_tree([lca], tree)

    status = 'nomatch'
    if not tree:
        return [], status

    # now find lowest-common-ancestor of the resulting tree.
    lca, reason = lca_utils.find_lca(tree)
    if reason == 0:               # leaf node
        debug('END', lca)
        status = 'found'
    else:                         # internal node => disagreement
        debug('MULTI', lca)
        status = 'disagree'

    debug('lineage is:', lca)

    return lca, status


def classify(args):
    """
    main single-genome classification function.
    """
    p = SourmashArgumentParser(prog="sourmash lca classify")
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV to this file instead of stdout')
    p.add_argument('--scaled', type=float)
    p.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath directories.')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    if not args.query:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    # have to have two calls as python < 3.5 can only have one expanded list
    if not check_files_exist(*args.query):
        sys.exit(-1)

    if not check_files_exist(*args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)

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
        for query_sig in load_signatures(query_filename, ksize=ksize):
            notify(u'\r\033[K', end=u'')
            notify('... classifying {} (file {} of {})', query_sig.name(),
                   n, total_n, end='\r')
            debug('classifying', query_sig.name())
            total_count += 1

            # make sure we're looking at the same scaled value as database
            query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

            # do the classification
            lineage, status = classify_signature(query_sig, dblist,
                                                 args.threshold)
            debug(lineage)

            # output each classification to the spreadsheet
            row = [query_sig.name(), status]
            row += lca_utils.zip_lineage(lineage)

            # when outputting to stdout, make output intelligible
            if not args.output:
                notify(u'\r\033[K', end=u'')
            csvfp.writerow(row)

    notify(u'\r\033[K', end=u'')
    notify('classified {} signatures total', total_count)


if __name__ == '__main__':
    sys.exit(classify(sys.argv[1:]))
