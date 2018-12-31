"""
Command-line entry point for 'python -m sourmash.sig'
"""
from __future__ import print_function, unicode_literals
import sys
import argparse
import csv
import json

import sourmash
import copy

from ..logging import set_quiet, error, notify, set_quiet
from .. import sourmash_args
from ..sourmash_args import DEFAULT_LOAD_K

usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

info <signature> [<signature> ... ]       - provide basic info on signatures
downsample <signature> [<signature> ... ] - downsample one or more signatures
extract <signature> [<signature> ... ]    - extract one or more signatures
flatten <signature> [<signature> ... ]    - remove abundances
intersect <signature> [<signature> ...]   - intersect one or more signatures
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
subtract <signature> <other_sig> [...]    - subtract one or more signatures
import [ ... ]                            - import a mash or other signature
export <signature>                        - export a signature, e.g. to mash

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''


def _check_abundance_compatibility(sig1, sig2):
    if sig1.minhash.track_abundance != sig2.minhash.track_abundance:
        raise ValueError("incompatible signatures: track_abundance is {} in first sig, {} in second".format(sig1.minhash.track_abundance, sig2.minhash.track_abundance))


def info(args):
    """
    provide basic info on signatures
    """
    p = argparse.ArgumentParser(prog='sourmash signature info')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('--csv', type=argparse.FileType('wt'),
                   help='output information to a CSV file')

    args = p.parse_args(args)
    set_quiet(args.quiet)

    siglist = []
    for sigfile in args.signatures:
        this_siglist = []
        try:
            this_siglist = list(sourmash.load_signatures(sigfile, quiet=True,
                                                         do_raise=True))
        except Exception as e:
            error('Error reading signatures from {}; skipping'.format(sigfile))

        for k in this_siglist:
            siglist.append((k, sigfile))
        notify('loaded {} signatures from {}...', len(this_siglist), sigfile,
               end='\r')

    notify('loaded {} signatures total.', len(siglist))

    # write CSV?
    w = None
    if args.csv:
        w = csv.DictWriter(args.csv,
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
        filename = sig.d.get('filename', '')
        license = sig.d['license']

        if w:
            w.writerow(locals())

        if not args.quiet:
            print('''\
---
signature filename: {signature_file}
signature: {name}
source file: {filename}
md5: {md5}
k={ksize} molecule={moltype} num={num} scaled={scaled} seed={seed} track_abundance={with_abundance}
size: {n_hashes}
signature license: {license}
'''.format(**locals()))


def merge(args):
    """
    merge one or more signatures.
    """
    p = argparse.ArgumentParser(prog='sourmash signature merge')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)

    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # load a first signature; record compatibility of with-abundance.
    first_sigfile = args.signatures[0]
    first_sig = sourmash.load_one_signature(first_sigfile, ksize=args.ksize, select_moltype=moltype)
    notify('loaded signature from {}...', first_sigfile, end='\r')
    total_loaded = 1

    mh = copy.copy(first_sig.minhash)

    # merge each successive 
    for sigfile in args.signatures[1:]:
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)
        try:
            _check_abundance_compatibility(first_sig, sigobj)

            mh.merge(sigobj.minhash)
        except:
            error("ERROR when merging signature '{}' ({}) from file {}",
                  sigobj.name(), sigobj.md5sum()[:8], sigfile)
            raise

        notify('loaded and merged signature from {}...', sigfile, end='\r')
        total_loaded += 1

    merged_sigobj = sourmash.SourmashSignature(mh)

    output_json = sourmash.save_signatures([merged_sigobj], fp=args.output)

    notify('loaded and merged {} signatures', total_loaded)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.

    This function always removes abundances.
    """
    p = argparse.ArgumentParser(prog='sourmash signature intersect')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    first_sigfile = args.signatures[0]
    first_sig = sourmash.load_one_signature(first_sigfile, ksize=args.ksize, select_moltype=moltype)
    notify('loaded signature from {}...', first_sigfile, end='\r')
    total_loaded = 1

    mins = set(first_sig.minhash.get_mins())

    for sigfile in args.signatures[1:]:
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)
        mins.intersection_update(sigobj.minhash.get_mins())
        notify('loaded and intersected signature from {}...', sigfile, end='\r')
        total_loaded += 1

    # forcibly turn off track_abundance
    intersect_mh = first_sig.minhash.copy_and_clear()
    new_mh_params = list(intersect_mh.__getstate__())
    new_mh_params[5] = False
    intersect_mh.__setstate__(new_mh_params)
    assert not intersect_mh.track_abundance
    intersect_mh.add_many(mins)
    intersect_sigobj = sourmash.SourmashSignature(intersect_mh)

    output_json = sourmash.save_signatures([intersect_sigobj], fp=args.output)

    notify('loaded and intersected {} signatures', total_loaded)


def subtract(args):
    """
    subtract one or more signatures from another
    """
    p = argparse.ArgumentParser(prog='sourmash signature subtract')
    p.add_argument('signature_from')
    p.add_argument('subtraction_sigs', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    from_sigfile = args.signature_from
    from_sigobj = sourmash.load_one_signature(from_sigfile, ksize=args.ksize, select_moltype=moltype)

    from_mh = from_sigobj.minhash
    if from_mh.track_abundance:
        error('Cannot use subtract on signatures with abundance tracking, sorry!')
        sys.exit(1)

    subtract_mins = set(from_mh.get_mins())

    notify('loaded signature from {}...', from_sigfile, end='\r')

    total_loaded = 0
    for sigfile in args.subtraction_sigs:
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)

        if sigobj.minhash.track_abundance:
            error('Cannot use subtract on signatures with abundance tracking, sorry!')
            sys.exit(1)

        subtract_mins -= set(sigobj.minhash.get_mins())

        notify('loaded and subtracted signature from {}...', sigfile, end='\r')
        total_loaded += 1

    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(subtract_mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    output_json = sourmash.save_signatures([subtract_sigobj], fp=args.output)

    notify('loaded and subtracted {} signatures', total_loaded)


def rename(args):
    """
    rename a signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature rename')
    p.add_argument('signature')
    p.add_argument('name')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    sigobj = sourmash.load_one_signature(args.signature, ksize=args.ksize, select_moltype=moltype)
    sigobj.d['name'] = args.name
    output_json = sourmash.save_signatures([sigobj], fp=args.output)

    notify("set name to '{}'", args.name)


def extract(args):
    """
    extract a signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature extract')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    p.add_argument('--md5', default=None,
                   help='select signatures whose md5 contains this substring')
    p.add_argument('--name', default=None,
                   help='select signatures whose name contains this substring')

    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

        outlist.extend(siglist)

    output_json = sourmash.save_signatures(outlist, fp=args.output)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(outlist),
           len(args.signatures))


def flatten(args):
    """
    flatten a signature, removing abundances.
    """
    p = argparse.ArgumentParser(prog='sourmash signature flatten')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    p.add_argument('--md5', default=None,
                   help='select signatures whose md5 contains this substring')
    p.add_argument('--name', default=None,
                   help='select signatures whose name contains this substring')

    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_signatures(filename, ksize=args.ksize,
                                           select_moltype=moltype)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

        for ss in siglist:
            flattened_mh = ss.minhash.copy_and_clear()
            new_mh_params = list(flattened_mh.__getstate__())
            new_mh_params[5] = False
            flattened_mh.__setstate__(new_mh_params)
            assert not flattened_mh.track_abundance
            flattened_mh.add_many(ss.minhash.get_mins())

            ss.minhash = flattened_mh

        outlist.extend(siglist)

    output_json = sourmash.save_signatures(outlist, fp=args.output)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(outlist),
           len(args.signatures))


def downsample(args):
    """
    downsample a scaled signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature downsample')
    p.add_argument('signatures', nargs="+")
    p.add_argument('--scaled', type=int, default=0,
                   help='scaled value to downsample to')
    p.add_argument('--num', type=int, default=0,
                   help='num value to downsample to')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
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
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)
        mh = sigobj.minhash

        if args.scaled:
            if mh.scaled:
                mh = mh.downsample_scaled(args.scaled)
            else:
                pass
        elif args.num:
            if mh.num:
                mh = mh.downsample_n(args.num)
            else:
                pass

        sigobj.minhash = mh

        output_list.append(sigobj)
        notify('loaded and downsampled signature from {}...', sigfile, end='\r')
        total_loaded += 1

    output_json = sourmash.save_signatures(output_list, fp=args.output)

    notify("loaded and downsampled {} signatures", total_loaded)


def sig_import(args):
    """
    import a signature into sourmash format.
    """
    p = argparse.ArgumentParser(prog='sourmash signature import')
    p.add_argument('filenames', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    args = p.parse_args(args)
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

    sourmash.save_signatures(siglist, args.output)


def export(args):
    """
    export a signature to mash format
    """
    p = argparse.ArgumentParser(prog='sourmash signature export')
    p.add_argument('filename')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout,
                   help='output signature to this file')
    sourmash_args.add_ksize_arg(p, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(p)
    args = p.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    total_loaded = 0
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

    print(json.dumps(x), file=args.output)
    notify("exported signature {} ({})", ss.name(), ss.md5sum()[:8])


def main(sysv_args):
    set_quiet(False)

    commands = {'merge': merge,
                'intersect': intersect,
                'rename': rename,
                'extract': extract,
                'flatten': flatten,
                'downsample': downsample,
                'subtract': subtract,
                'import': sig_import,
                'export': export,
                'info': info}

    parser = argparse.ArgumentParser(
        description='signature file manipulation utilities', usage=usage)
    parser.add_argument('sig_command', nargs='?')
    args = parser.parse_args(sysv_args[0:1])

    if not args.sig_command:
        print(usage)
        sys.exit(1)

    if args.sig_command not in commands:
        error('Unrecognized command: {}', args.sig_command)
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.sig_command)
    cmd(sysv_args[1:])

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
