"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import csv
import json
import os
from collections import defaultdict, namedtuple, Counter
import json
import re

import screed
import sourmash
from sourmash.sourmash_args import FileOutput

from sourmash.logging import (set_quiet, error, notify, print_results, debug,
                              debug_literal)
from sourmash import sourmash_args
from sourmash.minhash import _get_max_hash_for_scaled
from sourmash.manifest import CollectionManifest


usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

cat <signature> [<signature> ... ]        - concatenate all signatures
describe <signature> [<signature> ... ]   - show details of signature
downsample <signature> [<signature> ... ] - downsample one or more signatures
extract <signature> [<signature> ... ]    - extract one or more signatures
filter <signature> [<signature> ... ]     - filter k-mers on abundance
flatten <signature> [<signature> ... ]    - remove abundances
intersect <signature> [<signature> ...]   - intersect one or more signatures
manifest <sig/db>                         - build a manifest
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
split <signatures> [<signature> ...]      - split signatures into single files
subtract <signature> <other_sig> [...]    - subtract one or more signatures
import [ ... ]                            - import a mash or other signature
export <signature>                        - export a signature, e.g. to mash
overlap <signature1> <signature2>         - see detailed comparison of sigs
check <locations> --picklist ...          - check picklist against (many) sigs
collect <locations> -o manifest.sqlmf     - collect sigs metadata into manifest

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''


def _check_abundance_compatibility(sig1, sig2):
    if sig1.minhash.track_abundance != sig2.minhash.track_abundance:
        raise ValueError("incompatible signatures: track_abundance is {} in first sig, {} in second".format(sig1.minhash.track_abundance, sig2.minhash.track_abundance))


def _extend_signatures_with_from_file(args, *, target_attr='signatures'):
    # extend input signatures with --from-file
    if args.from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.from_file)

        sigs = list(getattr(args, target_attr))
        sigs.extend(more_files)
        setattr(args, target_attr, sigs)


def _set_num_scaled(mh, num, scaled):
    "set num and scaled values on a MinHash object"
    mh_params = list(mh.__getstate__())
    # Number of hashes is 0th parameter
    mh_params[0] = num
    # Scale is 8th parameter
    mh_params[8] = _get_max_hash_for_scaled(scaled)
    mh.__setstate__(mh_params)
    assert mh.num == num
    assert mh.scaled == scaled


##### actual command line functions


def cat(args):
    """
    concatenate all signatures into one file.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)

    encountered_md5sums = defaultdict(int)   # used by --unique

    # open output for saving sigs
    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    _extend_signatures_with_from_file(args)

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force,
                                                pattern=pattern_search)
    for ss, sigloc in loader:
        md5 = ss.md5sum()
        encountered_md5sums[md5] += 1
        if args.unique and encountered_md5sums[md5] > 1:
            continue

        save_sigs.add(ss)

    notify(f'loaded {len(save_sigs)} signatures total.')
    if picklist:
        sourmash_args.report_picklist(args, picklist)

    save_sigs.close()

    notify(f'output {len(save_sigs)} signatures')

    multiple_md5 = [ 1 for cnt in encountered_md5sums.values() if cnt > 1 ]
    if multiple_md5:
        notify(f'encountered {sum(multiple_md5)} MinHashes multiple times')
        if args.unique:
            notify('...and removed the duplicates, because --unique was specified.')


def split(args):
    """
    split all signatures into individual files
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    output_names = set()
    output_scaled_template = '{md5sum}.k={ksize}.scaled={scaled}.{moltype}.dup={dup}.{basename}.sig'
    output_num_template = '{md5sum}.k={ksize}.num={num}.{moltype}.dup={dup}.{basename}.sig'

    if args.output_dir:
        if not os.path.exists(args.output_dir):
            notify(f'Creating --output-dir {args.output_dir}')
            os.mkdir(args.output_dir)

    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)

    for sig, sigloc in loader:
        # save each file individually --
        md5sum = sig.md5sum()[:8]
        minhash = sig.minhash
        basename = os.path.basename(sig.filename)
        if not basename or basename == '-':
            basename = 'none'

        params = dict(basename=basename,
                      md5sum=md5sum,
                      scaled=minhash.scaled,
                      ksize=minhash.ksize,
                      num=minhash.num,
                      moltype=minhash.moltype)

        if minhash.scaled:
            output_template = output_scaled_template
        else: # num
            assert minhash.num
            output_template = output_num_template

        # figure out if this is duplicate, build unique filename
        n = 0
        params['dup'] = n
        output_name = output_template.format(**params)
        while output_name in output_names:
            params['dup'] = n
            output_name = output_template.format(**params)
            n += 1

        output_names.add(output_name)

        if args.output_dir:
            output_name = os.path.join(args.output_dir, output_name)

        if os.path.exists(output_name):
            notify(f"** overwriting existing file {format(output_name)}")

        # save!
        with sourmash_args.SaveSignaturesToLocation(output_name) as save_sigs:
            save_sigs.add(sig)
            notify(f'writing sig to {output_name}')

    notify(f'loaded and split {len(progress)} signatures total.')
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def describe(args):
    """
    provide basic info on signatures
    """
    set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    _extend_signatures_with_from_file(args)

    # write CSV?
    w = None
    csv_obj = None
    if args.csv:
        csv_obj = sourmash_args.FileOutputCSV(args.csv)
        csv_fp = csv_obj.open()

        w = csv.DictWriter(csv_fp,
                           ['signature_file', 'md5', 'ksize', 'moltype',
                            'num', 'scaled', 'n_hashes', 'seed',
                            'with_abundance', 'name', 'filename', 'license',
                            'sum_hashes'],
                           extrasaction='ignore')
        w.writeheader()

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force,
                                                pattern=pattern_search)

    for sig, location in loader:
        # extract info, write as appropriate.
        signature_file = location
        mh = sig.minhash
        ksize = mh.ksize
        moltype = mh.moltype
        scaled = mh.scaled
        num = mh.num
        seed = mh.seed
        n_hashes = len(mh)
        sum_hashes = sum(mh.hashes.values())
        with_abundance = 0
        if mh.track_abundance:
            with_abundance = 1
        md5 = sig.md5sum()
        name = sig.name
        p_name = name or "** no name **"
        filename = sig.filename
        p_filename = filename or "** no name **"
        license = sig.license

        if w:
            w.writerow(locals())

        print_results('''\
---
signature filename: {location}
signature: {p_name}
source file: {p_filename}
md5: {md5}
k={ksize} molecule={moltype} num={num} scaled={scaled} seed={seed} track_abundance={with_abundance}
size: {n_hashes}
sum hashes: {sum_hashes}
signature license: {license}
''', **locals())

    if csv_obj:
        csv_obj.close()

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def manifest(args):
    """
    build a signature manifest
    """
    set_quiet(args.quiet, args.debug)

    try:
        loader = sourmash_args.load_file_as_index(args.location,
                                                  yield_all_files=args.force)
    except ValueError as exc:
        error(f"Cannot open '{args.location}' as a sourmash signature collection.")
        error("Use -d/--debug for details.")
        sys.exit(-1)

    rebuild = True
    if args.no_rebuild_manifest:
        debug("sig manifest: not forcing rebuild.")
        rebuild = False
    else:
        debug("sig manifest: forcing rebuild.")

    manifest = sourmash_args.get_manifest(loader, require=True,
                                          rebuild=rebuild)

    manifest.write_to_filename(args.output,
                               database_format=args.manifest_format,
                               ok_if_exists=args.force)
    notify(f"manifest contains {len(manifest)} signatures total.")
    notify(f"wrote manifest to '{args.output}' ({args.manifest_format})")


def overlap(args):
    """
    provide detailed comparison of two signatures
    """
    set_quiet(args.quiet)

    moltype = sourmash_args.calculate_moltype(args)

    sig1 = sourmash.load_one_signature(args.signature1, ksize=args.ksize,
                                       select_moltype=moltype)
    sig2 = sourmash.load_one_signature(args.signature2, ksize=args.ksize,
                                       select_moltype=moltype)

    notify(f'loaded one signature each from {args.signature1} and {args.signature2}')

    try:
        similarity = sig1.similarity(sig2)
    except ValueError:
        raise

    cont1 = sig1.contained_by(sig2)
    cont2 = sig2.contained_by(sig1)

    sig1_file = args.signature1
    sig2_file = args.signature2

    name1 = sig1.name
    name2 = sig2.name

    md5_1 = sig1.md5sum()
    md5_2 = sig2.md5sum()

    ksize = sig1.minhash.ksize
    moltype = sig1.minhash.moltype

    num = sig1.minhash.num
    size1 = len(sig1.minhash)
    size2 = len(sig2.minhash)

    scaled = sig1.minhash.scaled

    hashes_1 = set(sig1.minhash.hashes)
    hashes_2 = set(sig2.minhash.hashes)

    num_common = len(hashes_1 & hashes_2)
    disjoint_1 = len(hashes_1 - hashes_2)
    disjoint_2 = len(hashes_2 - hashes_1)
    num_union = len(hashes_1.union(hashes_2))

    print('''\
first signature:
  signature filename: {sig1_file}
  signature: {name1}
  md5: {md5_1}
  k={ksize} molecule={moltype} num={num} scaled={scaled}

second signature:
  signature filename: {sig2_file}
  signature: {name2}
  md5: {md5_2}
  k={ksize} molecule={moltype} num={num} scaled={scaled}

similarity:                  {similarity:.5f}
first contained in second:   {cont1:.5f}
second contained in first:   {cont2:.5f}

number of hashes in first:   {size1}
number of hashes in second:  {size2}

number of hashes in common:  {num_common}
only in first:               {disjoint_1}
only in second:              {disjoint_2}
total (union):               {num_union}
'''.format(**locals()))


def merge(args):
    """
    merge one or more signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    first_sig = None
    mh = None

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)

    for sigobj, sigloc in loader:
        # first signature? initialize a bunch of stuff
        if first_sig is None:
            first_sig = sigobj
            mh = first_sig.minhash.copy_and_clear()

            # forcibly remove abundance?
            if args.flatten:
                mh.track_abundance = False

        try:
            sigobj_mh = sigobj.minhash
            if not args.flatten:
                _check_abundance_compatibility(first_sig, sigobj)
            else:
                sigobj_mh.track_abundance = False

            mh.merge(sigobj_mh)
        except (TypeError, ValueError) as exc:
            error("ERROR when merging signature '{}' ({}) from file {}",
                  sigobj, sigobj.md5sum()[:8], sigloc)
            error(str(exc))
            sys.exit(-1)

    if not len(progress):
        error("no signatures to merge!?")
        sys.exit(-1)

    merged_sigobj = sourmash.SourmashSignature(mh, name=args.name)

    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        save_sigs.add(merged_sigobj)

    notify(f'loaded and merged {len(progress)} signatures')

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.

    This function always removes abundances unless -A specified.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    first_sig = None
    mins = None

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)

    for sigobj, sigloc in loader:
        if first_sig is None:
            first_sig = sigobj
            mins = set(sigobj.minhash.hashes)
        else:
            # check signature compatibility -- if no ksize/moltype specified
            # 'first_sig' may be incompatible with later sigs.
            if not sigobj.minhash.is_compatible(first_sig.minhash):
                error("incompatible minhashes; specify -k and/or molecule type.")
                sys.exit(-1)

        mins.intersection_update(sigobj.minhash.hashes)

    # forcibly turn off track_abundance, unless --abundances-from set.
    intersect_mh = first_sig.minhash.copy_and_clear().flatten()
    intersect_mh.add_many(mins)

    # borrow abundances from a signature?
    if args.abundances_from:
        notify(f'loading signature from {args.abundances_from}, keeping abundances')
        abund_sig = sourmash.load_one_signature(args.abundances_from,
                                                ksize=args.ksize,
                                                select_moltype=moltype)
        if not abund_sig.minhash.track_abundance:
            error("--track-abundance not set on loaded signature?! exiting.")
            sys.exit(-1)

        intersect_mh = intersect_mh.inflate(abund_sig.minhash)

    intersect_sigobj = sourmash.SourmashSignature(intersect_mh)
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        save_sigs.add(intersect_sigobj)

    notify(f'loaded and intersected {len(progress)} signatures')
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def inflate(args):
    """
    inflate one or more other signatures from the first.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    inflate_sig = sourmash_args.load_query_signature(args.signature_from,
                                                     ksize=args.ksize,
                                                     select_moltype=moltype)
    inflate_from_mh = inflate_sig.minhash
    ksize = inflate_from_mh.ksize
    moltype = inflate_from_mh.moltype

    if not inflate_from_mh.track_abundance:
        error(f"ERROR: signature '{inflate_sig.name}' from ")
        error(f"file '{args.signature_from}' has no abundances.")
        sys.exit(-1)

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.other_sigs,
                                                ksize=ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)

    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        for sigobj, sigloc in loader:
            inflated_mh = sigobj.minhash.inflate(inflate_from_mh)
            inflated_sigobj = sourmash.SourmashSignature(inflated_mh,
                                                         name=sigobj.name)

            save_sigs.add(inflated_sigobj)

    if len(progress) == 0:
        error("no signatures to inflate!?")
        sys.exit(-1)

    notify(f'loaded and intersected {len(save_sigs)} signatures')
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def subtract(args):
    """
    subtract one or more signatures from another
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    from_sigfile = args.signature_from
    from_sigobj = sourmash.load_one_signature(from_sigfile, ksize=args.ksize, select_moltype=moltype)

    if args.abundances_from:    # it's ok to work with abund signatures if -A.
        args.flatten = True

    from_mh = from_sigobj.minhash
    if from_mh.track_abundance and not args.flatten:
        error('Cannot use subtract on signatures with abundance tracking, sorry!')
        sys.exit(1)

    subtract_mins = set(from_mh.hashes)

    notify(f'loaded signature from {from_sigfile}...', end='\r')

    progress = sourmash_args.SignatureLoadingProgress()

    for sigfile in args.subtraction_sigs:
        for sigobj in sourmash_args.load_file_as_signatures(sigfile,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress):
            if not sigobj.minhash.is_compatible(from_mh):
                error("incompatible minhashes; specify -k and/or molecule type.")
                sys.exit(-1)

            if sigobj.minhash.track_abundance and not args.flatten:
                error('Cannot use subtract on signatures with abundance tracking, sorry!')
                sys.exit(1)

            subtract_mins -= set(sigobj.minhash.hashes)

            notify(f'loaded and subtracted signatures from {sigfile}...', end='\r')

    if not len(progress):
        error("no signatures to subtract!?")
        sys.exit(-1)

    # build new minhash with new mins
    subtract_mh = from_sigobj.minhash.copy_and_clear().flatten()
    subtract_mh.add_many(subtract_mins)

    # borrow abundances from somewhere?
    if args.abundances_from:
        notify(f'loading signature from {args.abundances_from}, keeping abundances')
        abund_sig = sourmash.load_one_signature(args.abundances_from,
                                                ksize=args.ksize,
                                                select_moltype=moltype)
        if not abund_sig.minhash.track_abundance:
            error("--track-abundance not set on loaded signature?! exiting.")
            sys.exit(-1)

        subtract_mh = subtract_mh.inflate(abund_sig.minhash)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        save_sigs.add(subtract_sigobj)

    notify(f'loaded and subtracted {len(progress)} signatures')


def rename(args):
    """
    rename one or more signatures.
    """
    set_quiet(args.quiet, args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    _extend_signatures_with_from_file(args)

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force,
                                                pattern=pattern_search)

    for sigobj, sigloc in loader:
        sigobj = sigobj.to_mutable()
        sigobj._name = args.name
        save_sigs.add(sigobj)

    save_sigs.close()

    notify(f"set name to '{args.name}' on {len(save_sigs)} signatures")
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def extract(args):
    """
    extract signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    _extend_signatures_with_from_file(args)

    # further filtering on md5 or name?
    filter_fn = None
    if args.md5 is not None or args.name is not None:
        def filter_fn(row):
            # match?
            keep = False
            if args.name:
                name = row['name'] or row['filename']
                if args.name in name:
                    keep = True
            if args.md5 and args.md5 in row['md5']:
                keep = True

            return keep

    # ok! filtering defined, let's go forward
    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    # start loading!
    total_rows_examined = 0
    for filename in args.signatures:
        idx = sourmash_args.load_file_as_index(filename,
                                               yield_all_files=args.force)

        idx = idx.select(ksize=args.ksize, moltype=moltype)

        idx = sourmash_args.apply_picklist_and_pattern(idx, picklist,
                                                       pattern_search)

        manifest = sourmash_args.get_manifest(idx)
        total_rows_examined += len(manifest)

        # do the extra pattern matching on name/md5 that is part of 'extract'.
        # CTB: This should be deprecated and removed at some point, since
        # --include/--exclude now do the same thing.
        if filter_fn and not pattern_search:
            sub_manifest = manifest.filter_rows(filter_fn)
            sub_picklist = sub_manifest.to_picklist()

            try:
                idx = idx.select(picklist=sub_picklist)
            except ValueError:
                error("** This input collection doesn't support 'extract' with picklists or patterns.")
                error("** EXITING.")
                error("**")
                error("** You can use 'sourmash sig cat' with a picklist or pattern,")
                error("** and then pipe the output to 'sourmash sig extract")
                sys.exit(-1)

        for ss in idx.signatures():
            save_sigs.add(ss)

    notify(f"loaded {total_rows_examined} total that matched ksize & molecule type")
    if not save_sigs:
        error("no matching signatures to save!")
        sys.exit(-1)

    save_sigs.close()

    notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def filter(args):
    """
    filter hashes by abundance in all of the signatures
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    for filename in args.signatures:
        siglist = sourmash_args.load_file_as_signatures(filename,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)
        siglist = list(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in str(ss) ]

        for ss in siglist:
            mh = ss.minhash
            if not mh.track_abundance:
                notify(f'ignoring signature {ss} - track_abundance not set.')
                continue

            abunds = mh.hashes
            abunds2 = {}
            for k, v in abunds.items():
                if v >= args.min_abundance:
                    if args.max_abundance is None or \
                       v <= args.max_abundance:
                       abunds2[k] = v

            filtered_mh = mh.copy_and_clear()
            filtered_mh.set_abundances(abunds2)

            ss = ss.to_mutable()
            ss.minhash = filtered_mh

            save_sigs.add(ss)

    save_sigs.close()

    notify(f"loaded {len(progress)} total that matched ksize & molecule type")
    notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")


def flatten(args):
    """
    flatten one or more signatures, removing abundances.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)
    for ss, sigloc in loader:
        # select!
        if args.md5 is not None:
            if args.md5 not in ss.md5sum():
                continue        #  skip

        if args.name is not None:
            if args.name not in ss.name:
                continue        # skip

        ss = ss.to_mutable()
        ss.minhash = ss.minhash.flatten()
        save_sigs.add(ss)

    save_sigs.close()

    notify(f"loaded {len(progress)} total that matched ksize & molecule type")
    notify(f"extracted {len(save_sigs)} signatures from {len(args.signatures)} file(s)")
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def downsample(args):
    """
    downsample num and scaled signatures, and also interconvert.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    if not args.num_hashes and not args.scaled:
        error('ERROR: must specify either --num or --scaled value')
        sys.exit(-1)

    if args.num_hashes and args.scaled:
        error('ERROR: cannot specify both --num and --scaled')
        sys.exit(-1)

    # open output for saving sigs
    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)
    for ss, sigloc in loader:
        sigobj = ss.to_mutable()
        mh = sigobj.minhash

        if args.scaled:
            # downsample scaled to scaled? straightforward.
            if mh.scaled:
                mh_new = mh.downsample(scaled=args.scaled)
            # try to turn a num into a scaled - trickier.
            else:
                # first check: can we?
                max_hash = _get_max_hash_for_scaled(args.scaled)
                mins = mh.hashes
                if max(mins) < max_hash:
                    raise ValueError("this num MinHash does not have enough hashes to convert it into a scaled MinHash.")

                mh_new = mh.copy()
                _set_num_scaled(mh_new, 0, args.scaled)
        elif args.num_hashes:
            # downsample num to num? straightforward.
            if mh.num:
                mh_new = mh.downsample(num=args.num_hashes)
            # try to turn a scaled into a num - trickier.
            else:
                # first check: can we?
                if len(mh) < args.num_hashes:
                    raise ValueError(f"this scaled MinHash has only {len(mh)} hashes")

                mh_new = mh.copy()
                _set_num_scaled(mh_new, args.num_hashes, 0)


        sigobj.minhash = mh_new
        save_sigs.add(sigobj)

    save_sigs.close()

    notify(f"loaded {len(progress)} signatures")
    notify(f"output {len(save_sigs)} downsampled signatures", len(save_sigs))
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def sig_import(args):
    """
    import a signature into sourmash format.
    """
    set_quiet(args.quiet)

    siglist = []
    if args.csv:
        for filename in args.filenames:
            with open(filename, newline='') as csv_fp:
                reader = csv.reader(csv_fp)
                siglist = []
                for row in reader:
                    hashfn = row[0]
                    hashseed = int(row[1])

                    # only support a limited import type, for now ;)
                    assert hashfn == 'murmur64'
                    assert hashseed == 42

                    _, _, ksize, name, hashes = row
                    ksize = int(ksize)

                    hashes = hashes.strip()
                    hashes = list(map(int, hashes.split(' ' )))

                    e = sourmash.MinHash(len(hashes), ksize)
                    e.add_many(hashes)
                    s = sourmash.SourmashSignature(e, filename=name)
                    siglist.append(s)
                    notify(f'loaded signature: {name} {s.md5sum()[:8]}')
    else:
        for filename in args.filenames:
            with open(filename) as fp:
                x = json.loads(fp.read())

            ksize = x['kmer']
            num = x['sketchSize']

            assert x['hashType'] == "MurmurHash3_x64_128"
            assert x['hashBits'] == 64
            assert x['hashSeed'] == 42

            xx = x['sketches'][0]
            hashes = xx['hashes']

            mh = sourmash.MinHash(ksize=ksize, n=num, is_protein=False)
            mh.add_many(hashes)

            s = sourmash.SourmashSignature(mh, filename=filename)
            siglist.append(s)

    notify(f'saving {len(siglist)} signatures to JSON')
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        save_sigs.add_many(siglist)


def export(args):
    """
    export a signature to mash format
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    query = sourmash_args.load_query_signature(args.filename,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               select_md5=args.md5)
    mh = query.minhash

    x = {}
    x['kmer'] = mh.ksize
    x['sketchSize'] = len(mh)

    x['hashType'] = "MurmurHash3_x64_128"
    x['hashBits'] = 64
    x['hashSeed'] = mh.seed

    ll = list(mh.hashes)
    x['sketches'] = [{ 'hashes': ll }]

    with FileOutput(args.output, 'wt') as fp:
        print(json.dumps(x), file=fp)
    notify(f"exported signature {query} ({query.md5sum()[:8]})")


def kmers(args):
    """
    retrieve k-mers and/or sequences contained by the minhashes
    """
    from sourmash.search import format_bp

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    first_sig = None
    query_mh = None


    # start loading!
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_many_signatures(args.signatures,
                                                ksize=args.ksize,
                                                moltype=moltype,
                                                picklist=picklist,
                                                progress=progress,
                                                yield_all_files=args.force,
                                                force=args.force)

    for sigobj, sigloc in loader:
        # first signature? initialize a bunch of stuff
        if first_sig is None:
            first_sig = sigobj
            query_mh = first_sig.minhash.copy_and_clear()

            # remove abundance as it has no purpose here --
            query_mh.track_abundance = False

        try:
            sigobj_mh = sigobj.minhash
            sigobj_mh.track_abundance = False

            query_mh.merge(sigobj_mh)
        except (TypeError, ValueError) as exc:
            error("ERROR when merging signature '{}' ({}) from file {}",
                  sigobj, sigobj.md5sum()[:8], sigloc)
            error(str(exc))
            sys.exit(-1)

    if not len(progress):
        error("no signatures in query!?")
        sys.exit(-1)

    notify(f"loaded and merged {len(progress)} signatures")
    if picklist:
        sourmash_args.report_picklist(args, picklist)

    is_protein = False
    if query_mh.moltype == 'DNA':
        if args.translate:
            error("ERROR: cannot use --translate with DNA sketches.")
            sys.exit(-1)
    else:
        is_protein = True
        if args.translate:      # input sequence is DNA
            is_protein = False

    if not query_mh:
        notify("ERROR: no hashes in query signature!?")
        sys.exit(-1)

    notify("")
    notify(f"merged signature has the following properties:")
    notify(f"k={query_mh.ksize} molecule={query_mh.moltype} num={query_mh.num} scaled={query_mh.scaled} seed={query_mh.seed}")
    notify(f"total hashes in merged signature: {len(query_mh)}")
    notify("")
    notify("now processing sequence files for matches!")

    found_mh = query_mh.copy_and_clear()

    # open outputs...
    save_kmers = None
    kmer_w = None
    if args.save_kmers:
        save_kmers = sourmash_args.FileOutputCSV(args.save_kmers)
        save_kmers.open()
        kmer_w = csv.DictWriter(save_kmers.fp,
                                fieldnames=['sequence_file',
                                            'sequence_name',
                                            'kmer',
                                            'hashval'])
        kmer_w.writeheader()

    save_seqs = None
    if args.save_sequences:
        save_seqs = sourmash_args.FileOutput(args.save_sequences)
        save_seqs.open()

    # figure out protein vs dna
    is_protein = False
    if query_mh.moltype != 'DNA':
        if not args.translate:
            is_protein = True

    n_files_searched = 0
    n_sequences_searched = 0
    n_bp_searched = 0
    n_kmers_found = 0
    n_sequences_found = 0
    n_bp_saved = 0

    progress_threshold = 1e6
    progress_interval = 1e6
    for filename in args.sequences:
        notify(f"opening sequence file '{filename}'")
        n_files_searched += 1

        for record in screed.open(filename):
            seq_mh = query_mh.copy_and_clear()

            # protein? dna?
            if is_protein:
                seq_mh.add_protein(record.sequence)
            else:
                try:
                    seq_mh.add_sequence(record.sequence,
                                        not args.check_sequence)
                except ValueError as exc:
                    seqname = record.name
                    if len(seqname) > 40:
                        seqname = seqname[:37] + '...'
                    notify(f"ERROR in sequence '{seqname}', file '{filename}'")
                    notify(str(exc))
                    if args.force:
                        notify("(continuing)")
                        continue
                    else:
                        sys.exit(-1)

            if seq_mh.intersection(query_mh):
                # match!

                # output matching sequences:
                if save_seqs:
                    save_seqs.fp.write(f">{record.name}\n{record.sequence}\n")
                    n_sequences_found += 1
                    n_bp_saved += len(record.sequence)

                # output matching k-mers:
                if kmer_w:
                    seq = record.sequence
                    kh_iter = seq_mh.kmers_and_hashes(seq, force=False,
                                                      is_protein=is_protein)
                    for kmer, hashval in kh_iter:
                        if hashval in query_mh.hashes:
                            found_mh.add_hash(hashval)
                            n_kmers_found += 1
                            d = dict(sequence_file=filename,
                                     sequence_name=record.name,
                                     kmer=kmer, hashval=hashval)
                            kmer_w.writerow(d)

                # add seq_mh to found_mh
                found_mh += seq_mh.intersection(query_mh)

            # provide progress indicator based on bp...
            n_sequences_searched += 1
            n_bp_searched += len(record.sequence)

            if n_bp_searched >= progress_threshold:
                notify(f"... searched {n_bp_searched} from {n_files_searched} files so far")
                while n_bp_searched >= progress_threshold:
                    progress_threshold += progress_interval

    # END major for loop. Now, clean up!
    if save_kmers:
        save_kmers.close()

    if save_seqs:
        save_seqs.close()

    if not n_sequences_searched:
        notify("ERROR: no sequences searched!?")
        sys.exit(-1)

    # ...and report!
    notify("DONE.")
    notify(f"searched {n_sequences_searched} sequences from {n_files_searched} files, containing a total of {format_bp(n_bp_searched)}.")

    if save_seqs:
        notify(f"matched and saved a total of {n_sequences_found} sequences with {format_bp(n_bp_saved)}.")

    if kmer_w:
        notify(f"matched and saved a total of {n_kmers_found} k-mers.")

    # calculate overlap, even for num minhashes which ordinarily don't
    # permit it, because here we are interested in knowing how many
    # of the expected hashes we found.
    query_hashes = set(query_mh.hashes)
    found_hashes = set(found_mh.hashes)
    cont = len(query_hashes.intersection(found_hashes)) / len(query_hashes)

    notify(f"found {len(found_mh)} distinct matching hashes ({cont*100:.1f}%)")

    if not kmer_w and not save_seqs:
        notify("NOTE: see --save-kmers or --save-sequences for output options.")


_SketchInfo = namedtuple('_SketchInfo', 'ksize, moltype, scaled, num, abund')


def _summarize_manifest(manifest):
    info_d = {}

    # use a namedtuple to track counts of distinct sketch types and n hashes
    total_size = 0
    counter = Counter()
    hashcounts = Counter()
    for row in manifest.rows:
        ski = _SketchInfo(ksize=row['ksize'], moltype=row['moltype'],
                          scaled=row['scaled'], num=row['num'],
                          abund=row['with_abundance'])
        counter[ski] += 1
        hashcounts[ski] += row['n_hashes']
        total_size += row['n_hashes']

    # store in info_d
    info_d['total_hashes'] = total_size
    sketch_info = []
    for ski, count in counter.items():
        sketch_d = dict(ski._asdict())
        sketch_d['count'] = count
        sketch_d['n_hashes'] = hashcounts[ski]
        sketch_info.append(sketch_d)
    info_d['sketch_info'] = sketch_info

    return info_d


# NOTE: also aliased as 'summarize'
def fileinfo(args):
    """
    provide summary information on the given path (collection, index, etc.)
    """
    set_quiet(args.quiet, args.debug)

    text_out = False
    if not args.json_out:
        text_out = True

    # load as index!
    try:
        notify(f"** loading from '{args.path}'")
        idx = sourmash_args.load_file_as_index(args.path,
                                               yield_all_files=args.force)
    except ValueError:
        error(f"Cannot open '{args.path}' as a sourmash signature collection.")
        error("Use -d/--debug for details.")
        sys.exit(-1)

    print_bool = lambda x: "yes" if x else "no"
    print_none = lambda x: "n/a" if x is None else x

    info_d = {}
    info_d['path_filetype'] = type(idx).__name__
    info_d['location'] = "" if not idx.location else idx.location
    info_d['is_database'] = bool(idx.is_database)
    info_d['has_manifest'] = bool(idx.manifest)
    info_d['num_sketches'] = len(idx)

    if text_out:
        print_results(f"path filetype: {info_d['path_filetype']}")
        print_results(f"location: {info_d['location']}")
        print_results(f"is database? {print_bool(info_d['is_database'])}")
        print_results(f"has manifest? {print_bool(info_d['has_manifest'])}")
        print_results(f"num signatures: {info_d['num_sketches']}")

    # also have arg to fileinfo to force recalculation
    notify("** examining manifest...")

    manifest = sourmash_args.get_manifest(idx, rebuild=args.rebuild_manifest,
                                          require=False)

    if manifest is None:
        # actually can't find any file type to trigger this, but leaving it
        # in for future eventualities, I guess?
        notify("** no manifest and cannot be generated; exiting.")
        sys.exit(0)

    info_d.update(_summarize_manifest(manifest))

    if text_out:
        print_results(f"total hashes: {info_d['total_hashes']}")
        print_results("summary of sketches:")

        for ski in info_d['sketch_info']:
            mh_type = f"num={ski['num']}" if ski['num'] else f"scaled={ski['scaled']}"
            mh_abund = ", abund" if ski['abund'] else ""

            sketch_str = f"{ski['count']} sketches with {ski['moltype']}, k={ski['ksize']}, {mh_type}{mh_abund}"

            print_results(f"   {sketch_str: <50} {ski['n_hashes']} total hashes")

    else:
        assert args.json_out
        print(json.dumps(info_d))


def check(args):
    """
    check signature db(s) against a picklist.
    """
    from sourmash.picklist import PickStyle
    set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    _extend_signatures_with_from_file(args)

    if not picklist:
        error("** No picklist provided?! Exiting.")
        sys.exit(-1)

    if picklist.pickstyle == PickStyle.EXCLUDE and args.output_missing:
        error("** ERROR: Cannot use an 'exclude' picklist with '-o/--output-missing'")
        sys.exit(-1)

    # require manifests?
    require_manifest = True
    if args.no_require_manifest:
        require_manifest = False
        debug("sig check: manifest will not be required")
    else:
        debug("sig check: manifest required")

    total_manifest_rows = CollectionManifest([])

    # start loading!
    total_rows_examined = 0
    for filename in args.signatures:
        idx = sourmash_args.load_file_as_index(filename,
                                               yield_all_files=args.force)

        idx = idx.select(ksize=args.ksize, moltype=moltype)

        if idx.manifest is None and require_manifest:
            error(f"ERROR on filename '{filename}'.")
            error("sig check requires a manifest by default, but no manifest present.")
            error("specify --no-require-manifest to dynamically generate one.")
            sys.exit(-1)

        # has manifest, or ok to build (require_manifest=False) - continue!
        new_manifest = sourmash_args.get_manifest(idx, require=True)
        sub_manifest = new_manifest.select_to_manifest(picklist=picklist)
        total_rows_examined += len(new_manifest)

        # rewrite locations so that each signature can be found by filename
        # of its container; this follows `sig collect` logic.
        rows = []
        for row in sub_manifest.rows:
            row['internal_location'] = filename
            total_manifest_rows.add_row(row)

        debug_literal(f"examined {len(new_manifest)} new rows, found {len(sub_manifest)} matching rows")

    notify(f"loaded {total_rows_examined} signatures.")

    sourmash_args.report_picklist(args, picklist)

    # output picklist of non-matching in same format as input picklist
    n_missing = len(picklist.pickset - picklist.found)
    if args.output_missing and n_missing:
        pickfile = picklist.pickfile

        # go through the input file and pick out missing rows.
        n_input = 0
        n_output = 0

        with sourmash_args.FileInputCSV(pickfile) as r:
            with open(args.output_missing, "w", newline='') as outfp:
                w = csv.DictWriter(outfp, fieldnames=r.fieldnames)
                w.writeheader()

                for row in r:
                    n_input += 1
                    if not picklist.matched_csv_row(row):
                        n_output += 1
                        w.writerow(row)
        notify(f"saved {n_output} non-matching rows of {n_input} picklist rows to '{args.output_missing}'")
    elif args.output_missing:
        notify(f"(no remaining picklist entries; not saving to '{args.output_missing}')")

    # save manifest of matching!
    if args.save_manifest_matching and total_manifest_rows:
        mf = total_manifest_rows
        mf.write_to_filename(args.save_manifest_matching,
                             database_format=args.manifest_format)
        notify(f"wrote {len(mf)} matching manifest rows to '{args.save_manifest_matching}'")
    elif args.save_manifest_matching:
        notify(f"(not saving matching manifest to '{args.save_manifest_matching}' because no matches)")

    if args.fail_if_missing and n_missing:
        error("** ERROR: missing values, and --fail-if-missing requested. Exiting.")
        sys.exit(-1)


def collect(args):
    "Collect signature metadata across many locations, save to manifest"
    # TODO:
    # test what happens with directories :)
    set_quiet(False, args.debug)

    if os.path.exists(args.output):
        if args.merge_previous:
            pass
        else:
            error(f"ERROR: '{args.output}' already exists!")
            error(f"ERROR: please remove it, or use --merge-previous to merge")
            sys.exit(-1)
    elif args.merge_previous:
        notify(f"WARNING: --merge-previous specified, but output file '{args.output}' does not already exist?")

    # load previous manifest for --merge-previous. This gets tricky with
    # mismatched manifest types, which we forbid.
    try:
        if args.manifest_format == 'sql':
            # create on-disk manifest
            from sourmash.index.sqlite_index import SqliteCollectionManifest

            if args.merge_previous:
                collected_mf = SqliteCollectionManifest.create_or_open(args.output)
            else:
                collected_mf = SqliteCollectionManifest.create(args.output)
        else:
            # create in-memory manifest that will be saved as CSV
            assert args.manifest_format == 'csv'

            if args.merge_previous and os.path.exists(args.output):
                collected_mf = CollectionManifest.load_from_filename(args.output)
            else:
                collected_mf = CollectionManifest()

            if not isinstance(collected_mf, CollectionManifest):
                raise Exception
    except:
        error(f"ERROR loading '{args.output}' with --merge-previous. Is it of type {args.manifest_format}?")
        sys.exit(-1)

    if args.merge_previous:
        notify(f"merging new locations with {len(collected_mf)} previous rows.")

    # require manifests? yes by default, since generating can be slow.
    require_manifest = True
    if args.no_require_manifest:
        require_manifest = False
        debug("sig check: manifest will not be required")
    else:
        debug("sig check: manifest required")

    n_files = 0

    # load from_file
    _extend_signatures_with_from_file(args, target_attr='locations')

    # convert to abspath
    if args.abspath:
        args.locations = [ os.path.abspath(iloc) for iloc in args.locations ]

    # iterate through, loading all the manifests from all the locations.
    for n_files, loc in enumerate(args.locations):
        notify(f"Loading signature information from {loc}.")

        if n_files % 100 == 0:
            notify(f'... loaded {len(collected_mf)} sigs from {n_files} files')
        idx = sourmash.load_file_as_index(loc)
        if idx.manifest is None and require_manifest:
            error(f"ERROR on location '{loc}'")
            error(f"sig collect requires a manifest by default, but no manifest present.")
            error("specify --no-require-manifest to dynamically generate one.")
            sys.exit(-1)

        mf = sourmash_args.get_manifest(idx)

        rows = []
        for row in mf.rows:
            row['internal_location'] = loc
            collected_mf.add_row(row)

    if args.manifest_format == 'csv':
        collected_mf.write_to_filename(args.output, database_format='csv',
                                       ok_if_exists=args.merge_previous)
    else:
        collected_mf.close()

    notify(f"saved {len(collected_mf)} manifest rows to '{args.output}'")

    return 0


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
