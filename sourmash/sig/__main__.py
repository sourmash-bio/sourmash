"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import argparse
import csv

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
intersect <signature> [<signature> ...]   - intersect one or more signatures
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
subtract <signature> <other_sig> [...]    - subtract one or more signatures

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
    p = argparse.ArgumentParser(prog='sourmash signature merge')
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
    p = argparse.ArgumentParser(prog='sourmash signature merge')
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
    subtract_mins = set(from_mh.get_mins())
    notify('loaded signature from {}...', from_sigfile, end='\r')

    total_loaded = 0
    for sigfile in args.subtraction_sigs:
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)
        try:
            _check_abundance_compatibility(from_sigobj, sigobj)
        except:
            error("ERROR when subtracting signature '{}' ({}) from file {}",
                  sigobj.name(), sigobj.md5sum()[:8], sigfile)
            raise

        subtract_mins -= set(sigobj.minhash.get_mins())
        notify('loaded and subtracted signature from {}...', sigfile, end='\r')
        total_loaded += 1

    subtract_mh = from_sigobj.minhash.copy_and_clear()
    if subtract_mh.track_abundance:
        values = {}
        for hashval, abund in from_mh.get_mins(with_abundance=True):
            if hashval in subtract_mins:   # only take hash values in subtract
                values[hashval] = abund
    else:
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


def downsample(args):
    """
    downsample a scaled signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature rename')
    p.add_argument('signatures', nargs="+")
    p.add_argument('--scaled', type=int, default=10000,
                   help='value to downsample to')
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

    output_list = []
    total_loaded = 0
    for sigfile in args.signatures:
        sigobj = sourmash.load_one_signature(sigfile, ksize=args.ksize, select_moltype=moltype)
        sigobj.minhash = sigobj.minhash.downsample_scaled(args.scaled)
        output_list.append(sigobj)
        notify('loaded and downsample signature from {}...', sigfile, end='\r')
        total_loaded += 1

    output_json = sourmash.save_signatures(output_list, fp=args.output)

    notify("loaded and downsampled {} signatures", total_loaded)


def main(sysv_args):
    set_quiet(False)

    commands = {'merge': merge,
                'intersect': intersect,
                'rename': rename,
                'extract': extract,
                'downsample': downsample,
                'subtract': subtract,
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
