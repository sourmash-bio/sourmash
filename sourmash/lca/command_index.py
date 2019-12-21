#! /usr/bin/env python
"""
Build a lowest-common-ancestor database with given taxonomy and genome sigs.
"""
from __future__ import print_function
import sys
import csv
from collections import defaultdict

from .. import sourmash_args, load_signatures
from ..logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import LineagePair
from ..sourmash_args import SourmashArgumentParser


def load_taxonomy_assignments(filename, delimiter=',', start_column=2,
                              use_headers=True, force=False):
    """
    Load a taxonomy assignment spreadsheet into a dictionary.

    The 'assignments' dictionary that's returned maps identifiers to
    lineage tuples.
    """
    mode = 'rt'
    if sys.version_info < (3, ):
        mode = 'rtU'

    # parse spreadsheet!
    fp = open(filename, mode)
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
    n_species = 0
    n_strains = 0
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
                # check duplicates
                if ident in assignments:
                    if assignments[ident] != tuple(lineage):
                        if not force:
                            raise Exception("multiple lineages for identifier {}".format(ident))
                else:
                    assignments[ident] = tuple(lineage)

                    if lineage[-1].rank == 'species':
                        n_species += 1
                    elif lineage[-1].rank == 'strain':
                        n_strains += 1

    fp.close()

    # this is to guard against a bug that happened once and I can't find
    # any more, when building a large GTDB-based database :) --CTB
    if len(assignments) * 0.2 > n_species and len(assignments) > 50:
        if not force:
            raise Exception("error: fewer than 20% of lineages have species-level resolution!? ({} total found)".format(n_species))

    return assignments, num_rows


def generate_report(record_duplicates, record_no_lineage, record_remnants,
                    unused_lineages, unused_identifiers, filename):
    """
    Output a report of anomalies from building the index.
    """
    with open(filename, 'wt') as fp:
        print('Duplicate signatures:', file=fp)
        fp.write("\n".join(record_duplicates))
        fp.write("\n")
        print('----\nUnused identifiers:', file=fp)
        fp.write("\n".join(unused_identifiers))
        fp.write("\n")
        print('----\nNo lineage provided for these identifiers:', file=fp)
        fp.write("\n".join(record_no_lineage))
        fp.write("\n")
        print('----\nNo signatures found for these identifiers:', file=fp)
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
    p = SourmashArgumentParser(prog="sourmash lca index")
    p.add_argument('csv', help='taxonomy spreadsheet')
    p.add_argument('lca_db_out', help='name to save database to')
    p.add_argument('signatures', nargs='+',
                   help='one or more sourmash signatures')
    p.add_argument('--scaled', default=10000, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
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
    p.add_argument('--require-taxonomy', action='store_true',
                   help='ignore signatures with no taxonomy entry')
    args = p.parse_args(args)

    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

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

    # convert identities to numbers.
    ident_to_idx = {}
    idx_to_lid = {}

    lid_to_lineage = {}
    lineage_to_lid = {}

    arg_d = dict(next_index=0, next_lid=0)          # hack to keep from using nonlocal

    def get_ident_index(ident, fail_on_duplicate=False, arg_d=arg_d):
        idx = ident_to_idx.get(ident)
        if fail_on_duplicate:
            assert idx is None     # should be no duplicate identities

        if idx is None:
            idx = arg_d['next_index']
            arg_d['next_index'] += 1

            ident_to_idx[ident] = idx

        return idx

    def get_lineage_id(lineage, arg_d=arg_d):
        # lineage -> id
        lid = lineage_to_lid.get(lineage)
        if lid is None:
            lid = arg_d['next_lid']
            arg_d['next_lid'] += 1
            
            lineage_to_lid[lineage] = lid
            lid_to_lineage[lid] = lineage

        return lid

    for (ident, lineage) in assignments.items():
        idx = get_ident_index(ident, fail_on_duplicate=True)
        lid = get_lineage_id(lineage)
        
        # index -> lineage id
        idx_to_lid[idx] = lid

    notify('{} distinct identities in spreadsheet out of {} rows.',
           len(idx_to_lid), num_rows)
    notify('{} distinct lineages in spreadsheet out of {} rows.',
           len(set(idx_to_lid.values())), num_rows)

    # load signatures, construct index of hashvals to ident
    hashval_to_idx = defaultdict(set)
    md5_to_name = {}
    ident_to_name = {}

#    notify('finding signatures...')
    if args.traverse_directory:
        yield_all_files = False           # only pick up *.sig files?
        if args.force:
            yield_all_files = True
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures,
                                                          yield_all_files=yield_all_files))
    else:
        inp_files = list(args.signatures)

    #
    # main loop, connecting lineage ID to signature.
    #

    n = 0
    total_n = len(inp_files)
    record_duplicates = set()
    record_no_lineage = set()
    record_remnants = set(ident_to_idx.keys())
    record_used_lineages = set()
    record_used_idents = set()
    n_skipped = 0
    for filename in inp_files:
        n += 1
        for sig in load_signatures(filename, ksize=args.ksize):
            notify(u'\r\033[K', end=u'')
            notify('\r... loading signature {} (file {} of {}); skipped {} so far', sig.name()[:30], n, total_n, n_skipped, end='')
            debug(filename, sig.name())

            if sig.md5sum() in md5_to_name:
                debug('WARNING: in file {}, duplicate md5sum: {}; skipping', filename, sig.md5sum())
                record_duplicates.add(filename)
                continue

            ident = sig.name()
            if args.split_identifiers: # hack for NCBI-style names, etc.
                ident = ident.split(' ')[0].split('.')[0]

            # store full name
            ident_to_name[ident] = sig.name()

            # store md5 -> name too
            md5_to_name[sig.md5sum()] = sig.name()
            
            # remove from our list of remnant lineages
            try:
                record_remnants.remove(ident)
            except KeyError:
                # @CTB
                pass

            record_used_idents.add(ident)

            # downsample to specified scaled; this has the side effect of
            # making sure they're all at the same scaled value!
            minhash = sig.minhash.downsample_scaled(args.scaled)

            # connect hashvals to identity (and maybe lineage)
            idx = get_ident_index(ident)
            lid = idx_to_lid.get(idx)

            lineage = None
            if lid is not None:
                lineage = lid_to_lineage.get(lid)

            if lineage is None:
                debug('WARNING: no lineage assignment for {}.', ident)
                record_no_lineage.add(ident)
            else:
                record_used_lineages.add(lineage)

            if lineage is None and args.require_taxonomy:
                debug('(skipping, because --require-taxonomy was specified)')
                n_skipped += 1
                continue

            for hashval in minhash.get_mins():
                hashval_to_idx[hashval].add(idx)

    notify(u'\r\033[K', end=u'')

    if n == 0:
        error('ERROR: no signatures found. ??')
        if args.traverse_directory and not args.force:
            error('(note, with --traverse-directory, you may want to use -f)')
        sys.exit(1)

    if not hashval_to_idx:
        error('ERROR: no hash values found - are there any signatures?')
        sys.exit(1)

    notify('{} assigned lineages out of {} distinct lineages in spreadsheet.',
           len(record_used_lineages), len(set(assignments.values())))
    unused_lineages = set(assignments.values()) - record_used_lineages

    notify('{} identifiers used out of {} distinct identifiers in spreadsheet.',
           len(record_used_idents), len(set(assignments)))

    # remove unused identifiers
    unused_identifiers = set(assignments) - record_used_idents
    for ident in unused_identifiers:
        assert ident not in ident_to_name
        idx = get_ident_index(ident)
        del ident_to_idx[ident]
        if idx in idx_to_lid:
            del idx_to_lid[idx]

    # remove unusued lineages and lids
    for lineage in unused_lineages:
        lid = lineage_to_lid[lineage]
        del lineage_to_lid[lineage]
        del lid_to_lineage[lid]

    # now, save!
    db_outfile = args.lca_db_out
    if not (db_outfile.endswith('.lca.json') or \
                db_outfile.endswith('.lca.json.gz')):
        db_outfile += '.lca.json'
    notify('saving to LCA DB: {}'.format(db_outfile))

    db = lca_utils.LCA_Database()
    db.ident_to_name = ident_to_name
    db.ident_to_idx = ident_to_idx
    db.idx_to_lid = idx_to_lid
    db.lineage_to_lid = lineage_to_lid
    db.lid_to_lineage = lid_to_lineage
    db.hashval_to_idx = hashval_to_idx
    
    db.ksize = int(args.ksize)
    db.scaled = int(args.scaled)

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

        if unused_identifiers:
            notify('WARNING: {} unused identifiers.', len(unused_identifiers))

        if args.report:
            notify("generating a report and saving in '{}'", args.report)
            generate_report(record_duplicates, record_no_lineage,
                            record_remnants, unused_lineages,
                            unused_identifiers, args.report)
        else:
            notify('(You can use --report to generate a detailed report.)')


if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
