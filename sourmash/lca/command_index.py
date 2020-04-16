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
from .lca_utils import LineagePair, LCA_Database
from ..sourmash_args import DEFAULT_LOAD_K


class LCA_Database_Creation(LCA_Database):
    def __init__(self, ksize, scaled):
        super().__init__()
        self.ksize = int(ksize)
        self.scaled = int(scaled)

        self._next_index = 0
        self._next_lid = 0
        self.ident_to_name = {}
        self.ident_to_idx = {}
        self.idx_to_lid = {}
        self.lineage_to_lid = {}
        self.lid_to_lineage = {}
        self.hashval_to_idx = defaultdict(set)

    def build_get_ident_index(self, ident, fail_on_duplicate=False):
        "Get (create if nec) a unique int id, idx, for each identifier."
        idx = self.ident_to_idx.get(ident)
        if fail_on_duplicate:
            assert idx is None     # should be no duplicate identities

        if idx is None:
            idx = self._next_index
            self._next_index += 1

            self.ident_to_idx[ident] = idx

        return idx

    def build_get_lineage_id(self, lineage):
        "Get (create if nec) a unique lineage ID for each LineagePair tuples."
        # does one exist already?
        lid = self.lineage_to_lid.get(lineage)

        # nope - create one. Increment next_lid.
        if lid is None:
            lid = self._next_lid
            self._next_lid += 1

            # build mappings
            self.lineage_to_lid[lineage] = lid
            self.lid_to_lineage[lid] = lineage

        return lid

    def insert_signature(self, ident, sig, require_lineage=False):
        "Add a new signature into the LCA database."
        # store full name
        self.ident_to_name[ident] = sig.name()

        # connect hashvals to identity (and maybe lineage)
        idx = self.build_get_ident_index(ident)
        lid = self.idx_to_lid.get(idx)

        lineage = None
        if lid is not None:
            lineage = self.lid_to_lineage.get(lid)

        if lineage is None and require_lineage:
            return None

        # downsample to specified scaled; this has the side effect of
        # making sure they're all at the same scaled value!
        minhash = sig.minhash.downsample_scaled(self.scaled)

        for hashval in minhash.get_mins():
            self.hashval_to_idx[hashval].add(idx)

        return lineage


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
                        n_species += 1
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
    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    args.scaled = int(args.scaled)

    if args.ksize is None:
        args.ksize = DEFAULT_LOAD_K

    # first, load taxonomy spreadsheet
    delimiter = ','
    if args.tabs:
        delimiter = '\t'
    assignments, num_rows = load_taxonomy_assignments(args.csv,
                                               delimiter=delimiter,
                                               start_column=args.start_column,
                                               use_headers=not args.no_headers,
                                               force=args.force)

    db = LCA_Database_Creation(args.ksize, args.scaled)

    # make integer identifiers for lineages and indices for everything
    # in the spreadsheet.
    for (ident, lineage) in assignments.items():
        # identifiers -> integer indices (idx)
        idx = db.build_get_ident_index(ident, fail_on_duplicate=True)
        # (LineagePairs*) -> integer lineage ids (lids)
        lid = db.build_get_lineage_id(lineage)
        
        # map idx to lid.
        db.idx_to_lid[idx] = lid

    notify('{} distinct identities in spreadsheet out of {} rows.',
           len(db.idx_to_lid), num_rows)
    notify('{} distinct lineages in spreadsheet out of {} rows.',
           len(set(db.idx_to_lid.values())), num_rows)

#    notify('finding signatures...')
    if args.traverse_directory:
        yield_all_files = False           # only pick up *.sig files?
        if args.force:
            yield_all_files = True
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures,
                                                          yield_all_files=yield_all_files))
    else:
        inp_files = list(args.signatures)

    # track duplicates
    md5_to_name = {}

    #
    # main loop, connecting lineage ID to signature.
    #

    n = 0
    total_n = len(inp_files)
    record_duplicates = set()
    record_no_lineage = set()
    record_remnants = set(db.ident_to_idx.keys())
    record_used_lineages = set()
    record_used_idents = set()
    n_skipped = 0
    for filename in inp_files:
        n += 1
        for sig in load_signatures(filename, ksize=args.ksize):
            notify(u'\r\033[K', end=u'')
            notify('\r... loading signature {} (file {} of {}); skipped {} so far', sig.name()[:30], n, total_n, n_skipped, end='')
            debug(filename, sig.name())

            # block off duplicates.
            if sig.md5sum() in md5_to_name:
                debug('WARNING: in file {}, duplicate md5sum: {}; skipping', filename, sig.md5sum())
                record_duplicates.add(filename)
                continue

            md5_to_name[sig.md5sum()] = sig.name()

            # parse identifier, potentially with splitting
            ident = sig.name()
            if args.split_identifiers: # hack for NCBI-style names, etc.
                ident = ident.split(' ')[0].split('.')[0]

            # add the signature into the database.
            lineage = db.insert_signature(ident, sig,
                                       require_lineage=args.require_taxonomy)

            # punt if no lineage and --require-taxonomy
            if lineage is None and args.require_taxonomy:
                debug('(skipping, because --require-taxonomy was specified)')
                n_skipped += 1
                continue

            # remove from our list of remaining lineages
            try:
                record_remnants.remove(ident)
            except KeyError:
                # @CTB
                pass

            # track ident as used
            record_used_idents.add(ident)

            # track lineage info - either no lineage, or this lineage used.
            if lineage is None:
                debug('WARNING: no lineage assignment for {}.', ident)
                record_no_lineage.add(ident)
            else:
                record_used_lineages.add(lineage)

    # end main add signatures loop

    notify(u'\r\033[K', end=u'')

    # check -- did we find any signatures?
    if n == 0:
        error('ERROR: no signatures found. ??')
        if args.traverse_directory and not args.force:
            error('(note, with --traverse-directory, you may want to use -f)')
        sys.exit(1)

    # check -- did the signatures we found have any hashes?
    if not db.hashval_to_idx:
        error('ERROR: no hash values found - are there any signatures?')
        sys.exit(1)

    # summarize:
    notify('{} assigned lineages out of {} distinct lineages in spreadsheet.',
           len(record_used_lineages), len(set(assignments.values())))
    unused_lineages = set(assignments.values()) - record_used_lineages

    notify('{} identifiers used out of {} distinct identifiers in spreadsheet.',
           len(record_used_idents), len(set(assignments)))

    # remove unused identifiers
    unused_identifiers = set(assignments) - record_used_idents
    for ident in unused_identifiers:
        assert ident not in db.ident_to_name
        idx = db.build_get_ident_index(ident)
        del db.ident_to_idx[ident]
        if idx in db.idx_to_lid:
            del db.idx_to_lid[idx]

    # remove unusued lineages and lids
    for lineage in unused_lineages:
        lid = db.lineage_to_lid[lineage]
        del db.lineage_to_lid[lineage]
        del db.lid_to_lineage[lid]

    # now, save!
    db_outfile = args.lca_db_out
    if not (db_outfile.endswith('.lca.json') or \
                db_outfile.endswith('.lca.json.gz')):
        db_outfile += '.lca.json'
    notify('saving to LCA DB: {}'.format(db_outfile))

    db.save(db_outfile)

    ## done!

    # output a record of stuff if requested/available:
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
