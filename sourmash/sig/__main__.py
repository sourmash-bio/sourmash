"""
Command-line entry point for 'python -m sourmash.sig'
"""
from __future__ import print_function, unicode_literals
import sys
import csv
import json

import sourmash
import copy
from sourmash.sourmash_args import FileOutput

from ..logging import set_quiet, error, notify, set_quiet, print_results, debug
from .. import sourmash_args
from .._minhash import get_max_hash_for_scaled

usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

describe <signature> [<signature> ... ]   - show details of signature
downsample <signature> [<signature> ... ] - downsample one or more signatures
extract <signature> [<signature> ... ]    - extract one or more signatures
filter <signature> [<signature> ... ]     - filter k-mers on abundance
flatten <signature> [<signature> ... ]    - remove abundances
intersect <signature> [<signature> ...]   - intersect one or more signatures
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
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


def _set_num_scaled(mh, num, scaled):
    "set num and scaled values on a MinHash object"
    mh_params = list(mh.__getstate__())
    # Number of hashes is 0th parameter
    mh_params[0] = num
    # Scale is 8th parameter
    mh_params[8] = get_max_hash_for_scaled(scaled)
    mh.__setstate__(mh_params)
    assert mh.num == num
    assert mh.scaled == scaled


##### actual command line functions


def describe(args):
    """
    provide basic info on signatures
    """
    set_quiet(args.quiet)

    siglist = []
    for sigfile in args.signatures:
        this_siglist = []
        try:
            this_siglist = sourmash.load_signatures(sigfile, quiet=True, do_raise=True)
            for k in this_siglist:
                siglist.append((k, sigfile))
        except Exception as exc:
            error('\nError while reading signatures from {}:'.format(sigfile))
            error(str(exc))
            error('(continuing)')

        notify('loaded {} signatures from {}...', len(siglist), sigfile,
               end='\r')

    notify('loaded {} signatures total.', len(siglist))

    # write CSV?
    w = None
    csv_fp = None
    if args.csv:
        csv_fp = open(args.csv, 'wt')
        w = csv.DictWriter(csv_fp,
                           ['signature_file', 'md5', 'ksize', 'moltype', 'num',
                            'scaled', 'n_hashes', 'seed', 'with_abundance',
                            'name', 'filename', 'license'],
                           extrasaction='ignore')
        w.writeheader()

    # extract info, write as appropriate.
    for (sig, signature_file) in siglist:
        mh = sig.minhash
        ksize = mh.ksize
        moltype = 'DNA'
        if mh.is_protein:
            if mh.dayhoff:
                moltype = 'dayhoff'
            elif mh.hp:
                moltype = 'hp'
            else:
                moltype = 'protein'
        scaled = mh.scaled
        num = mh.num
        seed = mh.seed
        n_hashes = len(mh)
        with_abundance = 0
        if mh.track_abundance:
            with_abundance = 1
        md5 = sig.md5sum()
        name = sig.name()
        filename = sig.filename
        license = sig.license

        if w:
            w.writerow(locals())

        print_results('''\
---
signature filename: {signature_file}
signature: {name}
source file: {filename}
md5: {md5}
k={ksize} molecule={moltype} num={num} scaled={scaled} seed={seed} track_abundance={with_abundance}
size: {n_hashes}
signature license: {license}
''', **locals())

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

    notify('loaded one signature each from {} and {}', args.signature1,
           args.signature2)

    try:
        similarity = sig1.similarity(sig2)
    except ValueError:
        raise

    cont1 = sig1.contained_by(sig2)
    cont2 = sig2.contained_by(sig1)

    sig1_file = args.signature1
    sig2_file = args.signature2

    name1 = sig1.name()
    name2 = sig2.name()

    md5_1 = sig1.md5sum()
    md5_2 = sig2.md5sum()

    ksize = sig1.minhash.ksize
    moltype = 'DNA'
    if sig1.minhash.is_protein:
        moltype = 'protein'

    num = sig1.minhash.num
    size1 = len(sig1.minhash)
    size2 = len(sig2.minhash)

    scaled = sig1.minhash.scaled

    hashes_1 = set(sig1.minhash.get_mins())
    hashes_2 = set(sig2.minhash.get_mins())

    num_common = len(hashes_1.intersection(hashes_2))
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

    first_sig = None
    mh = None
    total_loaded = 0

    # iterate over all the sigs from all the files.
    for sigfile in args.signatures:
        notify('loading signatures from {}...', sigfile, end='\r')
        this_n = 0
        for sigobj in sourmash.load_signatures(sigfile, ksize=args.ksize,
                                               select_moltype=moltype,
                                               do_raise=True):
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
            except:
                error("ERROR when merging signature '{}' ({}) from file {}",
                      sigobj.name(), sigobj.md5sum()[:8], sigfile)
                raise

            this_n += 1
            total_loaded += 1
        if this_n:
            notify('loaded and merged {} signatures from {}...', this_n, sigfile, end='\r')

    if not total_loaded:
        error("no signatures to merge!?")
        sys.exit(-1)

    merged_sigobj = sourmash.SourmashSignature(mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([merged_sigobj], fp=fp)

    notify('loaded and merged {} signatures', total_loaded)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.

    This function always removes abundances.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    first_sig = None
    mins = None
    total_loaded = 0

    for sigfile in args.signatures:
        for sigobj in sourmash.load_signatures(sigfile, ksize=args.ksize,
                                               select_moltype=moltype,
                                               do_raise=True):
            if first_sig is None:
                first_sig = sigobj
                mins = set(sigobj.minhash.get_mins())

            mins.intersection_update(sigobj.minhash.get_mins())
            total_loaded += 1
        notify('loaded and intersected signatures from {}...', sigfile, end='\r')

    if total_loaded == 0:
        error("no signatures to merge!?")
        sys.exit(-1)

    # forcibly turn off track_abundance, unless --abundances-from set.
    if not args.abundances_from:
        intersect_mh = first_sig.minhash.copy_and_clear()
        intersect_mh.track_abundance = False
        intersect_mh.add_many(mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)
    else:
        notify('loading signature from {}, keeping abundances',
               args.abundances_from)
        abund_sig = sourmash.load_one_signature(args.abundances_from,
                                                ksize=args.ksize,
                                                select_moltype=moltype)
        if not abund_sig.minhash.track_abundance:
            error("--track-abundance not set on loaded signature?! exiting.")
            sys.exit(-1)
        intersect_mh = abund_sig.minhash.copy_and_clear()
        abund_mins = abund_sig.minhash.get_mins(with_abundance=True)

        # do one last intersection
        mins.intersection_update(abund_mins)
        abund_mins = { k: abund_mins[k] for k in mins }

        intersect_mh.set_abundances(abund_mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([intersect_sigobj], fp=fp)

    notify('loaded and intersected {} signatures', total_loaded)


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

    subtract_mins = set(from_mh.get_mins())

    notify('loaded signature from {}...', from_sigfile, end='\r')

    total_loaded = 0
    for sigfile in args.subtraction_sigs:
        for sigobj in sourmash.load_signatures(sigfile, ksize=args.ksize,
                                               select_moltype=moltype,
                                               do_raise=True):

            if sigobj.minhash.track_abundance and not args.flatten:
                error('Cannot use subtract on signatures with abundance tracking, sorry!')
                sys.exit(1)

            subtract_mins -= set(sigobj.minhash.get_mins())

            notify('loaded and subtracted signatures from {}...', sigfile, end='\r')
            total_loaded += 1

    if not total_loaded:
        error("no signatures to subtract!?")
        sys.exit(-1)


    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(subtract_mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([subtract_sigobj], fp=fp)

    notify('loaded and subtracted {} signatures', total_loaded)


def rename(args):
    """
    rename one or more signatures.
    """
    set_quiet(args.quiet, args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    for filename in args.sigfiles:
        debug('loading {}', filename)
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype)

        for sigobj in siglist:
            sigobj._name = args.name
            outlist.append(sigobj)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(outlist, fp=fp)

    notify("set name to '{}' on {} signatures", args.name, len(outlist))


def extract(args):
    """
    extract signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype,
                                           do_raise=True)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

        outlist.extend(siglist)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    if not outlist:
        error("no matching signatures!")
        sys.exit(-1)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(outlist, fp=fp)

    notify("extracted {} signatures from {} file(s)", len(outlist),
           len(args.signatures))


def filter(args):
    """
    filter hashes by abundance in all of the signatures
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype,
                                           do_raise=True)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

        for ss in siglist:
            mh = ss.minhash
            if not mh.track_abundance:
                notify('ignoring signature {} - track_abundance not set.',
                       ss)
                continue

            abunds = mh.get_mins(with_abundance=True)
            abunds2 = {}
            for k, v in abunds.items():
                if v >= args.min_abundance:
                    if args.max_abundance is None or \
                       v <= args.max_abundance:
                       abunds2[k] = v

            filtered_mh = mh.copy_and_clear()
            filtered_mh.set_abundances(abunds2)

            ss.minhash = filtered_mh

        outlist.extend(siglist)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(outlist, fp=fp)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(outlist),
           len(args.signatures))


def flatten(args):
    """
    flatten a signature, removing abundances.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype,
                                           do_raise=True)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

        for ss in siglist:
            flattened_mh = ss.minhash.copy_and_clear()
            flattened_mh.track_abundance = False
            flattened_mh.add_many(ss.minhash.get_mins())

            ss.minhash = flattened_mh

        outlist.extend(siglist)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(outlist, fp=fp)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(outlist),
           len(args.signatures))


def downsample(args):
    """
    downsample a scaled signature.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if not args.num and not args.scaled:
        error('must specify either --num or --scaled value')
        sys.exit(-1)

    if args.num and args.scaled:
        error('cannot specify both --num and --scaled')
        sys.exit(-1)

    output_list = []
    total_loaded = 0
    for sigfile in args.signatures:
        siglist = sourmash.load_signatures(sigfile, ksize=args.ksize, select_moltype=moltype, do_raise=True)

        for sigobj in siglist:
            mh = sigobj.minhash

            notify('loading and downsampling signature from {}...', sigfile, end='\r')
            total_loaded += 1
            if args.scaled:
                if mh.scaled:
                    mh_new = mh.downsample_scaled(args.scaled)
                else:                         # try to turn a num into a scaled
                    # first check: can we?
                    max_hash = get_max_hash_for_scaled(args.scaled)
                    mins = mh.get_mins()
                    if max(mins) < max_hash:
                        raise ValueError("this num MinHash does not have enough hashes to convert it into a scaled MinHash.")

                    mh_new = copy.copy(mh)
                    _set_num_scaled(mh_new, 0, args.scaled)
            elif args.num:
                if mh.num:
                    mh_new = mh.downsample_n(args.num)
                else:                         # try to turn a scaled into a num
                    # first check: can we?
                    if len(mh) < args.num:
                        raise ValueError("this scaled MinHash has only {} hashes")

                    mh_new = copy.copy(mh)
                    _set_num_scaled(mh_new, args.num, 0)

            sigobj.minhash = mh_new

            output_list.append(sigobj)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(output_list, fp=fp)

    notify("loaded and downsampled {} signatures", total_loaded)


def sig_import(args):
    """
    import a signature into sourmash format.
    """
    set_quiet(args.quiet)

    siglist = []
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

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(siglist, fp)


def export(args):
    """
    export a signature to mash format
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    ss = sourmash.load_one_signature(args.filename, ksize=args.ksize,
                                     select_moltype=moltype)
    mh = ss.minhash

    x = {}
    x['kmer'] = mh.ksize
    x['sketchSize'] = len(mh)

    x['hashType'] = "MurmurHash3_x64_128"
    x['hashBits'] = 64
    x['hashSeed'] = mh.seed

    ll = list(mh.get_mins())
    x['sketches'] = [{ 'hashes': ll }]

    with FileOutput(args.output, 'wt') as fp:
        print(json.dumps(x), file=fp)
    notify("exported signature {} ({})", ss.name(), ss.md5sum()[:8])


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
