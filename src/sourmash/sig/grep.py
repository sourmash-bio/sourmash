"""
Command-line entry point for 'python -m sourmash.sig grep'
"""
import sys
import re

from sourmash import logging, sourmash_args
from sourmash.logging import notify, error, debug, print_results
from sourmash.manifest import CollectionManifest
from .__main__ import _extend_signatures_with_from_file


def main(args):
    """
    extract signatures by pattern match.
    """
    # basic argument parsing
    logging.set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    # build the search pattern
    pattern = args.pattern
    if args.ignore_case:
        pattern = re.compile(pattern, re.IGNORECASE)
    else:
        pattern = re.compile(pattern)

    if args.invert_match:
        search_pattern = lambda vals: all(not pattern.search(val) for val in vals)
    else:
        search_pattern = lambda vals: any(pattern.search(val) for val in vals)

    # require manifests?
    require_manifest = True
    if args.no_require_manifest:
        require_manifest = False
        debug("sig grep: manifest will not be required")
    else:
        debug("sig grep: manifest required")

    # are we doing --count? if so, enforce --silent so no sigs are printed.
    if args.count:
        args.silent = True

    # define output type: signatures, or no?
    if args.silent:
        notify("(no signatures will be saved because of --silent/--count).")
        save_sigs = sourmash_args.SaveSignaturesToLocation(None)
    else:
        notify(f"saving matching signatures to '{args.output}'")
        save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
        save_sigs.open()

    # are we outputting a CSV? if so, initialize that, too.
    csv_obj = None
    if args.csv:
        csv_obj = sourmash_args.FileOutputCSV(args.csv)
        csv_fp = csv_obj.open()
        CollectionManifest.write_csv_header(csv_fp)

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

        # find all matching rows.
        sub_manifest = manifest.filter_on_columns(search_pattern,
                                                  ["name", "filename", "md5"])
        total_rows_examined += len(manifest)

        # write out to CSV, if desired.
        if args.csv:
            sub_manifest.write_to_csv(csv_fp)

        # just print out number of matches?
        if args.count:
            print_results(f"{len(sub_manifest)} matches: {filename}")
        elif not args.silent:
            # nope - do output signatures. convert manifest to picklist, apply.
            sub_picklist = sub_manifest.to_picklist()

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
    # done with the big loop over all indexes!

    if args.silent:
        pass
    else:
        notify(f"loaded {total_rows_examined} total that matched ksize & molecule type")

        if save_sigs:
            notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")
            save_sigs.close()
        else:
            error("no matching signatures found!")
            sys.exit(-1)

    if args.csv:
        notify(f"wrote manifest containing all matches to CSV file '{args.csv}'")
        csv_obj.close()

    if picklist:
        sourmash_args.report_picklist(args, picklist)
