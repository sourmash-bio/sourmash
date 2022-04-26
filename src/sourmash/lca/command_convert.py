#! /usr/bin/env python
"""
Convert and/or combine databases.
"""
import sys
import csv
import os

from sourmash import sourmash_args
from sourmash.sourmash_args import DEFAULT_LOAD_K
from sourmash.logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import LineagePair, check_files_exist
from .lca_db import LCA_Database


def convert_main(args):
    """
    main function for converting an LCA database.
    """
    set_quiet(args.quiet, args.debug)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    args.db = [item for sublist in args.db for item in sublist]
    args.scaled = int(args.scaled)

    if args.ksize is None:
        args.ksize = DEFAULT_LOAD_K

    moltype = sourmash_args.calculate_moltype(args, default='DNA')

    if not check_files_exist(*args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    assert len(dblist) == 1
    db = dblist[0]

    db_outfile = db.regularize_filename(args.output_database,
                                        args.database_format)

    if os.path.exists(db_outfile):
        error(f"ERROR: output file {db_outfile} already exists. Not overwriting.")
        sys.exit(-1)

    notify(f"saving LCA DB in format {args.database_format} to file '{format(db_outfile)}'")

    # now, save!
    db.save(db_outfile, format=args.database_format,
            display_progress=not args.quiet)

    ## done!


if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
