#! /usr/bin/env python
"""
Classify individual signature files down to deepest possible node.
"""
import sys
import csv

from .. import sourmash_args
from ..sourmash_args import load_file_as_signatures
from ..logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import check_files_exist

DEFAULT_THRESHOLD=5                  # how many counts of a taxid at min


def classify_signature(query_sig, dblist, threshold, majority):
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
    assignments = lca_utils.gather_assignments(query_sig.minhash.hashes,
                                               dblist)

    # now convert to trees -> do LCA & counts
    counts = lca_utils.count_lca_for_assignments(assignments)
    debug(counts.most_common())

    # ok, we now have the LCAs for each hashval, and their number of
    # counts. Now build a tree across "significant" LCAs - those above
    # threshold.

    tree = {}

    if counts and majority:
        majority_vote, count = counts.most_common()[0]
        if count > threshold:
            lca_utils.build_tree([majority_vote], tree)
    else:
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
    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    if not check_files_exist(*args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)

    # find all the queries
    notify('finding query signatures...')
    inp_files = list(args.query)
    if args.query_from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.query_from_file)
        inp_files.extend(more_files)

    if not check_files_exist(*inp_files):
        sys.exit(-1)

    if not inp_files:
        error('Error! must specify at least one query signature with --query or --query-from-file')
        sys.exit(-1)

    # set up output
    csvfp = csv.writer(sys.stdout)
    notify(f"outputting classifications to {args.output}")
    with sourmash_args.FileOutputCSV(args.output) as outfp:
        csvfp = csv.writer(outfp)

        csvfp.writerow(['ID','status'] + list(lca_utils.taxlist()))

        # for each query, gather all the matches across databases
        total_count = 0
        n = 0
        total_n = len(inp_files)
        for query_filename in inp_files:
            n += 1
            for query_sig in load_file_as_signatures(query_filename,
                                                     ksize=ksize):
                notify(u'\r\033[K', end=u'')
                notify(f'... classifying {query_sig} (file {n} of {total_n})', end='\r')
                debug('classifying', query_sig)
                total_count += 1

                # make sure we're looking at the same scaled value as database
                query_sig.minhash = query_sig.minhash.downsample(scaled=scaled)

                # do the classification
                lineage, status = classify_signature(query_sig, dblist,
                                                     args.threshold, args.majority)
                debug(lineage)

                # output each classification to the spreadsheet
                row = [str(query_sig), status]
                row += lca_utils.zip_lineage(lineage)

                # when outputting to stdout, make output intelligible
                if not args.output:
                    notify(u'\r\033[K', end=u'')
                csvfp.writerow(row)

        notify(u'\r\033[K', end=u'')
        notify(f'classified {total_count} signatures total')


if __name__ == '__main__':
    sys.exit(classify(sys.argv[1:]))
