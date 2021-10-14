"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import csv
import json
import os
from collections import defaultdict

import screed
import sourmash
from sourmash.sourmash_args import FileOutput

from sourmash.logging import set_quiet, error, notify, print_results, debug
from sourmash import sourmash_args
from sourmash.minhash import _get_max_hash_for_scaled


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

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''


def _check_abundance_compatibility(sig1, sig2):
    if sig1.minhash.track_abundance != sig2.minhash.track_abundance:
        raise ValueError("incompatible signatures: track_abundance is {} in first sig, {} in second".format(sig1.minhash.track_abundance, sig2.minhash.track_abundance))


def _extend_signatures_with_from_file(args):
    # extend input signatures with --from-file
    if args.from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.from_file)
        args.signatures = list(args.signatures)
        args.signatures.extend(more_files)

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
                                                force=args.force)
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

    if args.outdir:
        if not os.path.exists(args.outdir):
            notify(f'Creating --outdir {args.outdir}')
            os.mkdir(args.outdir)

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

        if args.outdir:
            output_name = os.path.join(args.outdir, output_name)

        if os.path.exists(output_name):
            notify(f"** overwriting existing file {format(output_name)}")

        # save!
        with open(output_name, 'wt') as outfp:
            sourmash.save_signatures([sig], outfp)
            notify(f'writing sig to {output_name}')

    notify(f'loaded and split {len(progress)} signatures total.')
    if picklist:
        sourmash_args.report_picklist(args, picklist)


def describe(args):
    """
    provide basic info on signatures
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    _extend_signatures_with_from_file(args)

    # write CSV?
    w = None
    csv_fp = None
    if args.csv:
        # CTB: might want to switch to sourmash_args.FileOutputCSV here?
        csv_fp = open(args.csv, 'w', newline='')
        w = csv.DictWriter(csv_fp,
                           ['signature_file', 'md5', 'ksize', 'moltype', 'num',
                            'scaled', 'n_hashes', 'seed', 'with_abundance',
                            'name', 'filename', 'license'],
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
                                                force=args.force)

    for sig, location in loader:
        # extract info, write as appropriate.
        mh = sig.minhash
        ksize = mh.ksize
        moltype = mh.moltype
        scaled = mh.scaled
        num = mh.num
        seed = mh.seed
        n_hashes = len(mh)
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
signature license: {license}
''', **locals())

    if csv_fp:
        csv_fp.close()

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def manifest(args):
    """
    build a signature manifest
    """
    from sourmash.index import CollectionManifest

    set_quiet(args.quiet)

    # CTB: might want to switch to sourmash_args.FileOutputCSV here?
    csv_fp = open(args.output, 'w', newline='')

    CollectionManifest.write_csv_header(csv_fp)
    w = csv.DictWriter(csv_fp, fieldnames=CollectionManifest.required_keys)

    try:
        loader = sourmash_args.load_file_as_index(args.location,
                                                  yield_all_files=args.force)
    except Exception as exc:
        error('\nError while reading signatures from {}:'.format(args.location))
        error(str(exc))
        error('(continuing)')
        raise

    n = 0
    # Need to ignore existing manifests here! otherwise circularity...
    try:
        manifest_iter = loader._signatures_with_internal()
    except NotImplementedError:
        error("ERROR: manifests cannot be generated for this file.")
        sys.exit(-1)

    for n, (sig, parent, loc) in enumerate(manifest_iter):
        # extract info, write as appropriate.
        row = CollectionManifest.make_manifest_row(sig, loc,
                                                   include_signature=False)
        w.writerow(row)

    notify(f'built manifest for {n} signatures total.')

    if csv_fp:
        csv_fp.close()


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

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([merged_sigobj], fp=fp)

    notify(f'loaded and merged {len(progress)} signatures')

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.

    This function always removes abundances.
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
            # check signature compatibility --
            if not sigobj.minhash.is_compatible(first_sig.minhash):
                error("incompatible minhashes; specify -k and/or molecule type.")
                sys.exit(-1)

        mins.intersection_update(sigobj.minhash.hashes)

    if len(progress) == 0:
        error("no signatures to merge!?")
        sys.exit(-1)

    # forcibly turn off track_abundance, unless --abundances-from set.
    if not args.abundances_from:
        intersect_mh = first_sig.minhash.copy_and_clear()
        intersect_mh.track_abundance = False
        intersect_mh.add_many(mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)
    else:
        notify(f'loading signature from {args.abundances_from}, keeping abundances')
        abund_sig = sourmash.load_one_signature(args.abundances_from,
                                                ksize=args.ksize,
                                                select_moltype=moltype)
        if not abund_sig.minhash.track_abundance:
            error("--track-abundance not set on loaded signature?! exiting.")
            sys.exit(-1)
        intersect_mh = abund_sig.minhash.copy_and_clear()
        abund_mins = abund_sig.minhash.hashes

        # do one last intersection
        mins.intersection_update(abund_mins)
        abund_mins = { k: abund_mins[k] for k in mins }

        intersect_mh.set_abundances(abund_mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([intersect_sigobj], fp=fp)

    notify(f'loaded and intersected {len(progress)} signatures')
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

    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(subtract_mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([subtract_sigobj], fp=fp)

    notify(f'loaded and subtracted {len(progress)} signatures')


def rename(args):
    """
    rename one or more signatures.
    """
    set_quiet(args.quiet, args.quiet)
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

    for sigobj, sigloc in loader:
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
    _extend_signatures_with_from_file(args)

    # further filtering on md5 or name?
    if args.md5 is not None or args.name is not None:
        def filter_fn(ss):
            # match?
            keep = False
            if args.name and args.name in str(ss):
                keep = True
            if args.md5 and args.md5 in ss.md5sum():
                keep = True

            return keep
    else:
        # whatever comes out of the database is fine
        filter_fn = lambda x: True

    # ok! filtering defined, let's go forward
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
        if filter_fn(ss):
            save_sigs.add(ss)

    notify(f"loaded {len(progress)} total that matched ksize & molecule type")
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
        mh = ss.minhash

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

        ss.minhash = mh_new

        # save!
        save_sigs.add(ss)

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
    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(siglist, fp)


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


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
