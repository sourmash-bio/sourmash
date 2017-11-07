#! /usr/bin/env python
"""
Build a least-common-ancestor database with given taxonomy and genome sigs.
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict, OrderedDict
import json

import sourmash_lib
from . import lca_utils
from ..logging import notify, error
from .. import sourmash_args
from .lca_utils import debug, set_debug


def index(args):
    """
    main function for building an LCA database.
    """
    p = argparse.ArgumentParser()
    p.add_argument('csv')
    p.add_argument('lca_db_out')
    p.add_argument('signatures', nargs='+')
    p.add_argument('--scaled', default=10000, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-C', '--start-column', default=2, type=int,
                   help='column at which taxonomic assignments start')
    p.add_argument('--split-identifiers', action='store_true',
                   help='split names in signatures on whitspace and period')
    p.add_argument('-f', '--force', action='store_true')
    p.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath directories.')
    args = p.parse_args(args)

    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    if args.debug:
        set_debug(args.debug)

    scaled = int(args.scaled)
    ksize = int(args.ksize)

    # parse spreadsheet!
    fp = open(args.csv, 'rt')
    r = csv.reader(fp)
    row_headers = ['identifiers']
    row_headers += ['_skip_']*(args.start_column - 2)
    row_headers += lca_utils.taxlist

    # first check that headers are interpretable.
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
                if not args.force:
                    sys.exit(-1)
                notify('...continue, because --force was specified.')

    # convert 
    assignments = {}
    num_rows = 0
    for row in r:
        if row and row[0].strip():        # want non-empty row
            num_rows += 1
            lineage = list(zip(row_headers, row))
            lineage = [ x for x in lineage if x[0] != '_skip_' ]

            ident = lineage[0][1]
            lineage = lineage[1:]

            # clean lineage of null names
            lineage = [(a,b) for (a,b) in lineage if b not in lca_utils.null_names]

            # store lineage tuple
            assignments[ident] = tuple(lineage)

    fp.close()

    # convert lineages to numbers.
    next_lineage_index = 0
    lineage_dict = {}

    assignments_idx = {}
    lineage_to_idx = {}
    for (ident, lineage_tuple) in assignments.items():
        idx = lineage_to_idx.get(lineage_tuple)
        if idx is None:
            idx = next_lineage_index
            next_lineage_index += 1

            lineage_dict[idx] = lineage_tuple
            lineage_to_idx[lineage_tuple] = idx

        assignments_idx[ident] = idx

    notify('{} distinct lineages in spreadsheet out of {} rows',
           len(lineage_dict), num_rows)

    # load signatures, construct index of hashvals to lineages
    hashval_to_lineage = defaultdict(set)
    md5_to_lineage = {}

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
    for filename in inp_files:
        n += 1
        for sig in sourmash_lib.load_signatures(filename, ksize=args.ksize):
            notify(u'\r\033[K', end=u'')
            notify('... loading signature {} (file {} of {})', sig.name()[:30], n, total_n, end='\r')
            debug(filename, sig.name())

            if sig.md5sum() in md5_to_lineage:
                notify('\nin file {},', filename)
                notify('duplicate md5sum: {}; skipping', sig.md5sum())
                continue

            name = sig.name()
            if args.split_identifiers:
                name = name.split(' ')[0].split('.')[0]

            # is this one for which we have a lineage assigned?
            lineage_idx = assignments_idx.get(name)
            if lineage_idx is not None:
                # downsample to specified scaled; this has the side effect of
                # making sure they're all at the same scaled value!
                minhash = sig.minhash.downsample_scaled(args.scaled)

                # connect hashvals to lineage
                for hashval in minhash.get_mins():
                    hashval_to_lineage[hashval].add(lineage_idx)

                # store md5 -> lineage too
                md5_to_lineage[sig.md5sum()] = lineage_idx

    notify(u'\r\033[K', end=u'')
    notify('...found {} genomes with lineage assignments!!',
           len(md5_to_lineage))

    # remove those lineages with no genomes associated
    assigned_lineages = set(md5_to_lineage.values())
    lineage_dict_2 = {}
    for idx in assigned_lineages:
        lineage_dict_2[idx] = lineage_dict[idx]

    notify('{} assigned lineages out of {} distinct lineages in spreadsheet',
           len(lineage_dict_2), len(lineage_dict))
    lineage_dict = lineage_dict_2

    # now, save!
    notify('saving to LCA DB: {}'.format(args.lca_db_out))

    db = lca_utils.LCA_Database()
    db.lineage_dict = lineage_dict
    db.hashval_to_lineage_id = hashval_to_lineage
    db.ksize = args.ksize
    db.scaled = args.scaled
    db.signatures_to_lineage = md5_to_lineage

    db.save(args.lca_db_out)

if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
