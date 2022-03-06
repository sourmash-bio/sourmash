"""
Command-line entry point for 'python -m sourmash.sig grep'
"""
import sys
import re

import sourmash
from sourmash import logging, sourmash_args
from sourmash.logging import notify, error, debug, print_results
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
        debug("sig grep: manifest will not be required")
    else:
        debug("sig grep: manifest required")

    if args.count:
        args.silent = True

    # define output
    if args.silent:
        notify("(no signatures will be output because of --silent/--count).")
        save_sigs = sourmash_args.SaveSignaturesToLocation(None)
    else:
        notify(f"saving matching signatures to '{args.output}'")
        save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
        save_sigs.open()

    csv_obj = None
    if args.csv:
        csv_obj = sourmash_args.FileOutputCSV(args.csv)
        csv_fp = csv_obj.open()
        CollectionManifest.write_csv_header(csv_fp)

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

        # get (and maybe generate) the manifest.
        manifest = idx.manifest
        if manifest is None:
            if require_manifest:
                error(f"ERROR on filename '{filename}'.")
                error("sig grep requires a manifest by default, but no manifest present.")
                error("specify --no-require-manifest to dynamically generate one.")
                sys.exit(-1)
            else:
                manifest = sourmash_args.get_manifest(idx,
                                                      require=False)

        # find all matching rows. CTB: could move this to manifest.py?
        sub_rows = []
        for row in manifest.rows:
            match = filter_fn(row)

            if match and not args.invert_match:
                sub_rows.append(row)
            elif not match and args.invert_match:
                sub_rows.append(row)

            total_rows_examined += 1

        # convert to picklist, grabme.
        sub_manifest = CollectionManifest(sub_rows)
        sub_picklist = sub_manifest.to_picklist()

        if args.csv:
            sub_manifest.write_to_csv(csv_fp)

        if args.count:
            print_results(f"{len(sub_rows)} matches: {filename}")
        elif not args.silent:
            try:
                idx = idx.select(picklist=sub_picklist)
            except ValueError:
                error("** This input collection doesn't support 'grep' with picklists.")
                error("** EXITING.")
                error("**")
                error("** You can use 'sourmash sig cat' with a picklist,")
                error("** and then pipe the output to 'sourmash sig grep -")
                sys.exit(-1)

            # save!
            for ss in idx.signatures():
                save_sigs.add(ss)

    if not args.silent:
        notify(f"loaded {total_rows_examined} total that matched ksize & molecule type")

    if not save_sigs and not args.silent:
        error("no matching signatures found!")
        sys.exit(-1)

    if not args.silent:
        notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")
    save_sigs.close()

    if args.csv:
        notify(f"wrote manifest containing matches to CSV file '{args.csv}'")
        csv_obj.close()

    if picklist:
        sourmash_args.report_picklist(args, picklist)
