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


DEFAULT_THRESHOLD=5


def summarize(hashvals, dblist, threshold, with_abundance):
    """
    Classify 'hashvals' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, counts) where 'lineage' is a tuple of LineagePairs.
    """

    # gather assignments from across all the databases
    assignments = lca_utils.gather_assignments(hashvals, dblist)

    # now convert to trees -> do LCA & counts
    if with_abundance:
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


def load_and_combine(filenames, ksize, scaled, with_abundance, traverse):
    "Load individual signatures and combine them all for classification."
    total_count = 0
    n = 0
    total_n = len(filenames)
    hashvals = defaultdict(int)
    for query_filename in filenames:
        n += 1
        for query_sig in sourmash_args.load_file_as_signatures(query_filename, ksize=ksize, traverse=traverse):
            notify(u'\r\033[K', end=u'')
            notify('... loading {} (file {} of {})', query_sig.name(), n,
                   total_n, end='\r')
            total_count += 1

            if with_abundance and not query_sig.minhash.track_abundance:
                notify("** error: minhash has no abundances, yet --with-abundance specified")
                sys.exit(-1)

            if not with_abundance and query_sig.minhash.track_abundance:
                notify("NOTE: discarding abundances in query, since --with-abundance not given")

            count_signature(query_sig, scaled, hashvals)

    notify(u'\r\033[K', end=u'')
    notify('loaded {} signatures from {} files total.', total_count, n)

    return hashvals


def load_singletons_and_count(filenames, ksize, scaled, with_abundance, traverse):
    "Load individual signatures and count them individually."
    total_count = 0
    n = 0

    # in order to get the right reporting out of this function, we need
    # to do our own traversal to expand the list of filenames, as opposed
    # to using load_file_as_signatures(..., traverse=True)
    if traverse:
        filenames = sourmash_args.traverse_find_sigs(filenames)
        filenames = list(filenames)

    total_n = len(filenames)

    for query_filename in filenames:
        n += 1
        for query_sig in sourmash_args.load_file_as_signatures(query_filename,
                                                               ksize=ksize):
            notify(u'\r\033[K', end=u'')
            notify('... loading {} (file {} of {})', query_sig.name(), n,
                   total_n, end='\r')
            total_count += 1

            if with_abundance and not query_sig.minhash.track_abundance:
                notify("** error: minhash has no abundances, yet --with-abundance specified")
                sys.exit(-1)

            if not with_abundance and query_sig.minhash.track_abundance:
                notify("NOTE: discarding abundances in query, since --with-abundance not given")

            # rebuild hashvals individually
            hashvals = defaultdict(int)
            count_signature(query_sig, scaled, hashvals)
            yield query_filename, query_sig, hashvals

    notify(u'\r\033[K', end=u'')
    notify('loaded {} signatures from {} files total.', total_count, n)


def count_signature(sig, scaled, hashvals):
    "Downsample sig to given scaled, count hashvalues."
    mh = sig.minhash.downsample_scaled(scaled)

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
    if filename or sig:                   # require both
        if not filename and sig:
            raise ValueError("must include both filename and sig arguments")

    for (lineage, count) in lineage_counts.items():
        if lineage:
            lineage = lca_utils.zip_lineage(lineage, truncate_empty=True)
            lineage = ';'.join(lineage)
        else:
            lineage = '(root)'

        p = count / total_counts * 100.
        p = '{:.1f}%'.format(p)

        if filename and sig:
            print_results('{:5} {:>5}   {}   {}:{} {}'.format(p, count, lineage, filename, sig.md5sum()[:8], sig.name()))
        else:
            print_results('{:5} {:>5}   {}'.format(p, count, lineage))


def output_csv(lineage_counts, csv_fp, filename, sig, write_header=True):
    """\
    Output results in CSV.
    """
    if filename or sig:                   # require both
        assert filename and sig

    w = csv.writer(csv_fp)
    if write_header:
        headers = ['count'] + list(lca_utils.taxlist())
        if filename:
            headers += ['filename', 'sig_name', 'sig_md5']
        w.writerow(headers)

    for (lineage, count) in lineage_counts.items():
        debug('lineage:', lineage)
        row = [count] + lca_utils.zip_lineage(lineage, truncate_empty=False)
        if filename:
            row += [filename, sig.name(), sig.md5sum()]
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

    with_abundance = args.with_abundance

    # flatten --db and --query lists
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    if not check_files_exist(*args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    if with_abundance:
        notify("Weighting output by k-mer abundances in query, since --with-abundance given.")

    # find all the queries
    notify('finding query signatures...')
    inp_files = args.query

    if args.query_from_file:
        more_files = sourmash_args.load_file_list_of_signatures(args.query_from_file)
        inp_files.extend(more_files)

    if not inp_files:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    if not check_files_exist(*inp_files):
        sys.exit(-1)

    if args.singleton:
        # summarize each signature individually
        csv_fp = None
        write_header = True
        if args.output:
            csv_fp = open(args.output, 'wt')

        try:
            for filename, sig, hashvals in \
              load_singletons_and_count(inp_files, ksize, scaled, with_abundance, args.traverse_directory):

                # get the full counted list of lineage counts in this signature
                lineage_counts = summarize(hashvals, dblist, args.threshold,
                                           with_abundance)
                if with_abundance:
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

    else:
        # load and merge all the signatures in all the files
        # DEPRECATE for 4.0.
        hashvals = load_and_combine(inp_files, ksize, scaled, with_abundance,
                                    args.traverse_directory)

        # get the full counted list of lineage counts across signatures
        lineage_counts = summarize(hashvals, dblist, args.threshold,
                                   with_abundance)

        # output!
        if with_abundance:
            total = float(sum(hashvals.values()))
        else:
            total = float(len(hashvals))

        output_results(lineage_counts, total)

        # CSV:
        if args.output:
            with sourmash_args.FileOutput(args.output, 'wt') as csv_fp:
                output_csv(lineage_counts, csv_fp, None, None)


if __name__ == '__main__':
    sys.exit(summarize_main(sys.argv[1:]))
