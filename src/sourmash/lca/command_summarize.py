#! /usr/bin/env python
"""
Summarize the taxonomic content of the given signatures, combined.
"""
import sys
import csv
from collections import defaultdict

from .. import sourmash_args
from ..logging import notify, error, print_results, set_quiet, debug
from . import lca_utils
from .lca_utils import check_files_exist
from sourmash.index import MultiIndex


DEFAULT_THRESHOLD=5


def summarize(hashvals, dblist, threshold, ignore_abundance):
    """
    Classify 'hashvals' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, counts) where 'lineage' is a tuple of LineagePairs.
    """

    # gather assignments from across all the databases
    assignments = lca_utils.gather_assignments(hashvals, dblist)

    # now convert to trees -> do LCA & counts
    if not ignore_abundance:
        counts = lca_utils.count_lca_for_assignments(assignments, hashvals)
    else: # flatten
        counts = lca_utils.count_lca_for_assignments(assignments, None)
    debug(counts.most_common())

    # ok, we now have the LCAs for each hashval, and their number
    # of counts. Now aggregate counts across the tree, going up from
    # the leaves.
    aggregated_counts = defaultdict(int)
    for lca, count in counts.most_common():
        if count < threshold:
            break

        if not lca:
            aggregated_counts[lca] += count

        # climb from the lca to the root.
        while lca:
            aggregated_counts[lca] += count
            lca = lca[:-1]

    debug(aggregated_counts)

    return aggregated_counts


def load_singletons_and_count(filenames, ksize, scaled, ignore_abundance):
    "Load individual signatures and count them individually."
    total_count = 0
    n = 0

    total_n = len(filenames)
    for filename in filenames:
        n += 1
        mi = MultiIndex.load_from_path(filename)
        mi = mi.select(ksize=ksize)

        for query_sig, query_filename in mi.signatures_with_location():
            notify(u'\r\033[K', end=u'')
            notify(f'... loading {query_sig} (file {n} of {total_n})',
                   total_n, end='\r')
            total_count += 1

            if ignore_abundance and query_sig.minhash.track_abundance:
                notify("NOTE: discarding abundances in query, since --ignore-abundance")

            # rebuild hashvals individually
            hashvals = defaultdict(int)
            count_signature(query_sig, scaled, hashvals)
            yield query_filename, query_sig, hashvals

    notify(u'\r\033[K', end=u'')
    notify(f'loaded {total_count} signatures from {n} files total.')


def count_signature(sig, scaled, hashvals):
    "Downsample sig to given scaled, count hashvalues."
    mh = sig.minhash.downsample(scaled=scaled)

    if mh.track_abundance:
        abunds = mh.hashes
        for hashval, count in abunds.items():
            hashvals[hashval] += count
    else:
        for hashval in mh.hashes:
            hashvals[hashval] += 1


def output_results(lineage_counts, total_counts, filename=None, sig=None):
    """\
    Output results in ~human-readable format.
    """

    for (lineage, count) in lineage_counts.items():
        if lineage:
            lineage = lca_utils.zip_lineage(lineage, truncate_empty=True)
            lineage = ';'.join(lineage)
        else:
            lineage = '(root)'

        p = count / total_counts * 100.
        p = '{:.1f}%'.format(p)

        print_results('{:5} {:>5}   {}   {}:{} {}'.format(p, count, lineage, filename, sig.md5sum()[:8], sig))

def output_csv(lineage_counts, csv_fp, filename, sig, write_header=True):
    """\
    Output results in CSV.
    """

    w = csv.writer(csv_fp)
    if write_header:
        headers = ['count'] + list(lca_utils.taxlist())
        headers += ['filename', 'sig_name', 'sig_md5']
        w.writerow(headers)

    for (lineage, count) in lineage_counts.items():
        debug('lineage:', lineage)
        row = [count] + lca_utils.zip_lineage(lineage, truncate_empty=False)
        row += [filename, sig.name, sig.md5sum()]
        w.writerow(row)


def summarize_main(args):
    """
    main summarization function.
    """
    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    ignore_abundance = args.ignore_abundance

    # flatten --db and --query lists
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    if not check_files_exist(*args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    if ignore_abundance:
        notify("Ignoring any k-mer abundances in query, since --ignore-abundance given.")

    # find all the queries
    notify('finding query signatures...')
    inp_files = args.query

    if args.query_from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.query_from_file)
        inp_files.extend(more_files)

    if not inp_files:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    if not check_files_exist(*inp_files):
        sys.exit(-1)

    # summarize each signature individually
    csv_fp = None
    write_header = True
    if args.output:
        csv_fp = open(args.output, 'w', newline='')

    try:
        for filename, sig, hashvals in \
          load_singletons_and_count(inp_files, ksize, scaled, ignore_abundance):

            # get the full counted list of lineage counts in this signature
            lineage_counts = summarize(hashvals, dblist, args.threshold,
                                       ignore_abundance)
            if not ignore_abundance:
                total = float(sum(hashvals.values()))
            else:
                total = float(len(hashvals))

            output_results(lineage_counts, total,
                           filename=filename, sig=sig)

            if csv_fp:
                output_csv(lineage_counts, csv_fp, filename, sig,
                           write_header=write_header)
                write_header = False
    finally:
        if csv_fp:
            csv_fp.close()


if __name__ == '__main__':
    sys.exit(summarize_main(sys.argv[1:]))
