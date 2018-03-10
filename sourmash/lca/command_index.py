#! /usr/bin/env python
"""
Build a lowest-common-ancestor database with given taxonomy and genome sigs.
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict

from .. import sourmash_args, load_signatures
from ..logging import notify, error
from . import lca_utils
from .lca_utils import debug, set_debug, LineagePair


def load_taxonomy_assignments(filename, delimiter=',', start_column=2,
                              use_headers=True, force=False):
    # parse spreadsheet!
    fp = open(filename, 'rtU')
    r = csv.reader(fp, delimiter=delimiter)
    row_headers = ['identifiers']
    row_headers += ['_skip_']*(start_column - 2)
    row_headers += list(lca_utils.taxlist())

    # first check that headers are interpretable.
    if use_headers:
        notify('examining spreadsheet headers...')
        first_row = next(iter(r))

        n_disagree = 0
        for (column, value) in zip(row_headers, first_row):
            if column == '_skip_':
                continue

            if column.lower() != value.lower():
                notify("** assuming column '{}' is {} in spreadsheet",
                       value, column)
                n_disagree += 1
                if n_disagree > 2:
                    error('whoa, too many assumptions. are the headers right?')
                    error('expecting {}', ",".join(row_headers))
                    if not force:
                        sys.exit(-1)
                    notify('...continue, because --force was specified.')

    # convert into a lineage pair
    assignments = {}
    num_rows = 0
    for row in r:
        if row and row[0].strip():        # want non-empty row
            num_rows += 1
            lineage = list(zip(row_headers, row))
            lineage = [ x for x in lineage if x[0] != '_skip_' ]

            ident = lineage[0][1]
            lineage = lineage[1:]

            # clean lineage of null names, replace with 'unassigned'
            lineage = [ (a, lca_utils.filter_null(b)) for (a,b) in lineage ]
            lineage = [ LineagePair(a, b) for (a, b) in lineage ]

            # remove end nulls
            while lineage and lineage[-1].name == 'unassigned':
                lineage = lineage[:-1]

            # store lineage tuple
            if lineage:
                assignments[ident] = tuple(lineage)

    fp.close()

    return assignments, num_rows


def generate_report(record_duplicates, record_no_lineage, record_remnants,
                    unused_lineages, filename):
    """
    Output a report of anomalies from building the index.
    """
    with open(filename, 'wt') as fp:
        print('Duplicate signatures:', file=fp)
        fp.write("\n".join(record_duplicates))
        fp.write("\n")
        print('----\nNo lineage provided for:', file=fp)
        fp.write("\n".join(record_no_lineage))
        fp.write("\n")
        print('----\nNo signatures found for these lineage assignments:', file=fp)
        fp.write('\n'.join(record_remnants))
        fp.write("\n")
        print('----\nUnused lineages:', file=fp)
        for lineage in unused_lineages:
            fp.write(";".join(lca_utils.zip_lineage(lineage)))
            fp.write("\n")


def index(args):
    """
    main function for building an LCA database.
    """
    p = argparse.ArgumentParser(prog="sourmash lca index")
    p.add_argument('csv', help='taxonomy spreadsheet')
    p.add_argument('lca_db_out', help='name to save database to')
    p.add_argument('signatures', nargs='+',
                   help='one or more sourmash signatures')
    p.add_argument('--scaled', default=10000, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-C', '--start-column', default=2, type=int,
                   help='column at which taxonomic assignments start')
    p.add_argument('--tabs', action='store_true',
                   help='input spreadsheet is tab-delimited (default: commas)')
    p.add_argument('--no-headers', action='store_true',
                   help='no headers present in taxonomy spreadsheet')
    p.add_argument('--split-identifiers', action='store_true',
                   help='split names in signatures on whitspace and period')
    p.add_argument('-f', '--force', action='store_true')
    p.add_argument('--traverse-directory', action='store_true',
                   help='load all signatures underneath directories.')
    p.add_argument('--report', help='output a report on anomalies, if any.')
    args = p.parse_args(args)

    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    if args.debug:
        set_debug(args.debug)

    args.scaled = int(args.scaled)

    # first, load taxonomy spreadsheet
    delimiter = ','
    if args.tabs:
        delimiter = '\t'
    assignments, num_rows = load_taxonomy_assignments(args.csv,
                                               delimiter=delimiter,
                                               start_column=args.start_column,
                                               use_headers=not args.no_headers,
                                               force=args.force)

    # convert lineages to numbers.
    next_lineage_index = 0
    lineage_dict = {}

    assignments_idx = {}
    lineage_to_idx = {}
    for (ident, lineage) in assignments.items():
        idx = lineage_to_idx.get(lineage)
        if idx is None:
            idx = next_lineage_index
            next_lineage_index += 1

            lineage_dict[idx] = lineage
            lineage_to_idx[lineage] = idx

        assignments_idx[ident] = idx

    notify('{} distinct lineages in spreadsheet out of {} rows',
           len(lineage_dict), num_rows)

    # load signatures, construct index of hashvals to lineages
    hashval_to_lineage = defaultdict(set)
    md5_to_lineage = {}
    md5_to_name = {}

    notify('finding signatures...')
    if args.traverse_directory:
        yield_all_files = False           # only pick up *.sig files?
        if args.force:
            yield_all_files = True
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures,
                                                          yield_all_files=yield_all_files))
    else:
        inp_files = list(args.signatures)

    n = 0
    total_n = len(inp_files)
    record_duplicates = set()
    record_no_lineage = set()
    record_remnants = set(assignments_idx.keys())
    for filename in inp_files:
        n += 1
        for sig in load_signatures(filename, ksize=args.ksize):
            notify(u'\r\033[K', end=u'')
            notify('... loading signature {} (file {} of {})', sig.name()[:30], n, total_n, end='\r')
            debug(filename, sig.name())

            if sig.md5sum() in md5_to_lineage:
                notify('\nWARNING: in file {}, duplicate md5sum: {}; skipping', filename, sig.md5sum())
                record_duplicates.add(filename)
                continue

            name = sig.name()
            if args.split_identifiers: # hack for NCBI-style names, etc.
                name = name.split(' ')[0].split('.')[0]

            # is this one for which we have a lineage assigned?
            lineage_idx = assignments_idx.get(name)
            if lineage_idx is None:
               notify('\nWARNING: no lineage assignment for {}.', name)
               record_no_lineage.add(name)
            else:
                # remove from our list of remnant lineages
                record_remnants.remove(name)

                # downsample to specified scaled; this has the side effect of
                # making sure they're all at the same scaled value!
                minhash = sig.minhash.downsample_scaled(args.scaled)

                # connect hashvals to lineage
                for hashval in minhash.get_mins():
                    hashval_to_lineage[hashval].add(lineage_idx)

                # store md5 -> lineage too
                md5_to_lineage[sig.md5sum()] = lineage_idx
                md5_to_name[sig.md5sum()] = sig.name()

    notify(u'\r\033[K', end=u'')
    notify('...found {} genomes with lineage assignments!!',
           len(md5_to_lineage))

    # remove those lineages with no genomes associated
    assigned_lineages = set(md5_to_lineage.values())
    lineage_dict_2 = {}
    for idx in assigned_lineages:
        lineage_dict_2[idx] = lineage_dict[idx]

    unused_lineages = set(lineage_dict.values()) - set(lineage_dict_2.values())

    notify('{} assigned lineages out of {} distinct lineages in spreadsheet',
           len(lineage_dict_2), len(lineage_dict))
    lineage_dict = lineage_dict_2

    # now, save!
    db_outfile = args.lca_db_out
    if not (db_outfile.endswith('.lca.json') or db_outfile.endswith('.lca.json.gz')):
        db_outfile += '.lca.json'
    notify('saving to LCA DB: {}'.format(db_outfile))

    db = lca_utils.LCA_Database()
    db.lineage_dict = lineage_dict
    db.hashval_to_lineage_id = hashval_to_lineage
    db.ksize = int(args.ksize)
    db.scaled = int(args.scaled)
    db.signatures_to_lineage_id = md5_to_lineage
    db.signatures_to_name = md5_to_name

    db.save(db_outfile)

    if record_duplicates or record_no_lineage or record_remnants or unused_lineages:
        if record_duplicates:
            notify('WARNING: {} duplicate signatures.', len(record_duplicates))
        if record_no_lineage:
            notify('WARNING: no lineage provided for {} signatures.',
                   len(record_no_lineage))
        if record_remnants:
            notify('WARNING: no signatures for {} lineage assignments.',
                   len(record_remnants))
        if unused_lineages:
            notify('WARNING: {} unused lineages.', len(unused_lineages))

        if args.report:
            notify("generating a report and saving in '{}'", args.report)
            generate_report(record_duplicates, record_no_lineage,
                            record_remnants, unused_lineages, args.report)
        else:
            notify('(You can use --report to generate a detailed report.)')


if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
