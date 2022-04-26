#! /usr/bin/env python
"""
Convert and/or combine databases.
"""
import sys
import csv
import os
from collections import defaultdict

from sourmash.logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import LineagePair
from .lca_db import LCA_Database


def convert_main(args):
    """
    main function for converting an LCA database.
    """
    return

    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    args.scaled = int(args.scaled)

    if args.ksize is None:
        args.ksize = DEFAULT_LOAD_K

    moltype = sourmash_args.calculate_moltype(args, default='DNA')
    picklist = sourmash_args.load_picklist(args)

    db_outfile = args.lca_db_out
    if args.database_format == 'json':
        if not (db_outfile.endswith('.lca.json') or \
                    db_outfile.endswith('.lca.json.gz')):   # logic -> db.save
            db_outfile += '.lca.json'
    else:
        assert args.database_format == 'sql'
        if not db_outfile.endswith('.lca.sql'):
                db_outfile += '.lca.sql'

    if os.path.exists(db_outfile):
        error(f"ERROR: output file {db_outfile} already exists. Not overwriting.")
        sys.exit(-1)

    notify(f'saving to LCA DB: {format(db_outfile)}')

    notify(f'Building LCA database with ksize={args.ksize} scaled={args.scaled} moltype={moltype}.')

    if picklist:
        sourmash_args.report_picklist(args, picklist)

    # now, save!
    db.save(db_outfile, format=args.database_format)

    ## done!


if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
