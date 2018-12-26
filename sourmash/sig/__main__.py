"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import argparse

import sourmash
import copy

from ..logging import set_quiet, error, notify, set_quiet

usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

downsample <signature> [<signature> ... ] - downsample one or more signatures
extract <signature> [<signature> ... ]    - extract one or more signatures
intersect <signature> [<signature> ...]   - intersect one or more signatures
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
subtract <signature> <other_sig> [...]    - subtract one or more signatures

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''


def merge(args):
    """
    merge one or more signatures.
    """
    p = argparse.ArgumentParser(prog='sourmash signature merge')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    args = p.parse_args(args)
    set_quiet(args.quiet)

    first_sigfile = args.signatures[0]
    first_sig = sourmash.load_one_signature(first_sigfile)
    notify('loaded signature from {}...', first_sigfile, end='\r')
    total_loaded = 1

    mh = copy.copy(first_sig.minhash)

    for sigfile in args.signatures[1:]:
        sigobj = sourmash.load_one_signature(sigfile)
        mh.merge(sigobj.minhash)
        notify('loaded and merged signature from {}...', sigfile, end='\r')
        total_loaded += 1

    merged_sigobj = sourmash.SourmashSignature(mh)

    output_json = sourmash.save_signatures([merged_sigobj])
    print(output_json)

    notify('loaded and merged {} signatures', total_loaded)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.
    """
    p = argparse.ArgumentParser(prog='sourmash signature merge')
    p.add_argument('signatures', nargs='+')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    args = p.parse_args(args)
    set_quiet(args.quiet)

    first_sigfile = args.signatures[0]
    first_sig = sourmash.load_one_signature(first_sigfile)
    notify('loaded signature from {}...', first_sigfile, end='\r')
    total_loaded = 1

    mins = set(first_sig.minhash.get_mins())

    for sigfile in args.signatures[1:]:
        sigobj = sourmash.load_one_signature(sigfile)
        mins.intersection_update(sigobj.minhash.get_mins())
        notify('loaded and intersected signature from {}...', sigfile, end='\r')
        total_loaded += 1

    intersect_mh = first_sig.minhash.copy_and_clear()
    intersect_mh.add_many(mins)
    intersect_sigobj = sourmash.SourmashSignature(intersect_mh)

    output_json = sourmash.save_signatures([intersect_sigobj])
    print(output_json)

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
    args = p.parse_args(args)
    set_quiet(args.quiet)

    from_sigfile = args.signature_from
    from_sigobj = sourmash.load_one_signature(from_sigfile)
    mins = set(from_sigobj.minhash.get_mins())
    notify('loaded signature from {}...', from_sigfile, end='\r')

    total_loaded = 0
    for sigfile in args.subtraction_sigs:
        sigobj = sourmash.load_one_signature(sigfile)
        mins -= set(sigobj.minhash.get_mins())
        notify('loaded and subtracted signature from {}...', sigfile, end='\r')
        total_loaded += 1

    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    output_json = sourmash.save_signatures([subtract_sigobj])
    print(output_json)

    notify('loaded and intersected {} signatures', total_loaded)


def rename(args):
    """
    rename a signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature rename')
    p.add_argument('signature')
    p.add_argument('name')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    args = p.parse_args(args)
    set_quiet(args.quiet)

    sigobj = sourmash.load_one_signature(args.signature)
    sigobj.d['name'] = args.name
    output_json = sourmash.save_signatures([sigobj])
    print(output_json)

    notify("set name to '{}'", args.name)


def extract(args):
    """
    extract a signature.
    """
    p = argparse.ArgumentParser(prog='sourmash signature extract')
    p.add_argument('signature')
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    args = p.parse_args(args)
    set_quiet(args.quiet)

    sigobj = sourmash.load_one_signature(args.signature)
    output_json = sourmash.save_signatures([sigobj])
    print(output_json)

    notify("extracted single signature from '{}'", args.signature)


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
    args = p.parse_args(args)
    set_quiet(args.quiet)

    output_list = []
    total_loaded = 0
    for sigfile in args.signatures:
        sigobj = sourmash.load_one_signature(sigfile)
        sigobj.minhash = sigobj.minhash.downsample_scaled(args.scaled)
        output_list.append(sigobj)
        notify('loaded and downsample signature from {}...', sigfile, end='\r')
        total_loaded += 1

    output_json = sourmash.save_signatures(output_list)
    print(output_json)

    notify("loaded and downsampled {} signatures", total_loaded)


def main(sysv_args):
    set_quiet(False)

    commands = {'merge': merge,
                'intersect': intersect,
                'rename': rename,
                'extract': extract,
                'downsample': downsample,
                'subtract': subtract}

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


if __name__ == '__main__':
    main(sys.argv[1:])
