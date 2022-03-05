"""
Command-line entry point for 'python -m sourmash.sig grep'
"""
import sys
import re

import sourmash
from sourmash import logging, sourmash_args
from sourmash.logging import notify, error
from sourmash.manifest import CollectionManifest
from .__main__ import _extend_signatures_with_from_file

def main(args):
    """
    extract signatures by pattern match.
    """
    logging.set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    pattern = args.pattern
    if args.ignore_case:
        pattern = re.compile(pattern, re.IGNORECASE)
    else:
        pattern = re.compile(pattern)

    # require manifests?
    require_manifest = True
    if args.no_require_manifest:
        require_manifest = False

    # define output
    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    # field retrieval from manifest row
    def get_fields(row):
        if row['name'] is not None:
            yield row['name']
        if row['filename'] is not None:
            yield row['filename']
        yield row['md5']

    # define filter function.
    def filter_fn(row):
        for field in get_fields(row):
            if pattern.search(field):
                return True
        return False

    # start loading!
    total_rows_examined = 0
    for filename in args.signatures:
        idx = sourmash_args.load_file_as_index(filename,
                                               yield_all_files=args.force)

        idx = idx.select(ksize=args.ksize,
                         moltype=moltype,
                         picklist=picklist)

        manifest = sourmash_args.get_manifest(idx, require=require_manifest)

        sub_rows = []
        for row in manifest.rows:
            match = filter_fn(row)

            if match and not args.invert_match:
                sub_rows.append(row)
            elif not match and args.invert_match:
                sub_rows.append(row)

            total_rows_examined += 1

        sub_manifest = CollectionManifest(sub_rows)
        sub_picklist = sub_manifest.to_picklist()

        try:
            idx = idx.select(picklist=sub_picklist)
        except ValueError:
            error("** This input collection doesn't support 'extract' with picklists.")
            error("** EXITING.")
            error("**")
            error("** You can use 'sourmash sig cat' with a picklist,")
            error("** and then pipe the output to 'sourmash sig extract")
            sys.exit(-1)

        for ss in idx.signatures():
            save_sigs.add(ss)

    notify(f"loaded {total_rows_examined} total that matched ksize & molecule type")
    if not save_sigs:
        error("no matching signatures found!")
        sys.exit(-1)

    save_sigs.close()

    notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")

    if picklist:
        sourmash_args.report_picklist(args, picklist)
