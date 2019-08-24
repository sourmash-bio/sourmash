"""
sourmash command line.
"""
from __future__ import print_function
import sys
import argparse

from .logging import error, set_quiet

from .commands import (categorize, compare, compute, dump, import_csv,
                       gather, index, sbt_combine, search,
                       plot, watch, info, storage, migrate, multigather)
from .lca import main as lca_main
from .sig import main as sig_main

try:
    from sourmash_utils.__main__ import main as utils_main
    UTILS_AVAILABLE = True
except ImportError:
    UTILS_AVAILABLE = False


usage='''
sourmash <command> [<args>]

** Commands include:

compute <filenames>         Compute MinHash signatures for sequences in files.
compare <filenames.sig>     Compute similarity matrix for multiple signatures.
search <query> <against>    Search a signature against a list of signatures.
plot <matrix>               Plot a distance matrix made by 'compare'.
gather                      Search a metagenome signature for multiple
                                 non-overlapping matches.

** Taxonomic classification utilities:

   Run 'sourmash lca' for the taxonomic classification routines.

** Sequence Bloom Tree (SBT) utilities:

index                   Index a collection of signatures for fast searching.
sbt_combine             Combine multiple SBTs into a new one.
categorize              Identify best matches for many signatures using an SBT.
watch                   Classify a stream of sequences.

** Other commands:

info                        Display sourmash version and other information.
signature                   Sourmash signature manipulation utilities.

Use '-h' to get subcommand-specific help, e.g.

sourmash compute -h

** Documentation is available at https://sourmash.readthedocs.io/
'''


def main():
    set_quiet(False)

    commands = {'search': search, 'compute': compute,
                'compare': compare, 'plot': plot,
                'import_csv': import_csv, 'dump': dump,
                'index': index,
                'categorize': categorize, 'gather': gather,
                'watch': watch,
                'sbt_combine': sbt_combine, 'info': info,
                'storage': storage,
                'lca': lca_main,
                'migrate': migrate,
                'multigather': multigather,
                'sig': sig_main,
                'signature': sig_main}

    loadable_cmds = {'compute'}

    if UTILS_AVAILABLE:
        commands['utils'] = utils_main

    parser = argparse.ArgumentParser(
        description='work with compressed biological sequence representations')
    parser.add_argument('command', nargs='?')
    args = parser.parse_args(sys.argv[1:2])

    if not args.command:
        print(usage)
        sys.exit(1)

    if args.command not in commands:
        error('Unrecognized command')
        print(usage)
        sys.exit(1)

    if args.command in loadable_cmds:
        parser = argparse.ArgumentParser(
            description='work with compressed sequence representations')
        subparsers = parser.add_subparsers()
        subp = subparsers.add_parser(args.command)
        commands[args.command].add_args(subp)
        subp.set_defaults(cmd=commands[args.command])
        args = parser.parse_args()
        return args.cmd(**vars(args))
    else:
        cmd = commands.get(args.command)
        return cmd(sys.argv[2:])


if __name__ == '__main__':
    main()
