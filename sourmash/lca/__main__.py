"""
Command-line entry point for 'python -m sourmash.lca'
"""

import sys
import argparse

from . import classify, index, summarize_main, rankinfo_main, gather_main
from .command_compare_csv import compare_csv
from ..logging import set_quiet, error

usage='''
sourmash lca <command> [<args>] - work with taxonomic information.

** Commands can be:

index <taxonomy.csv> <output_db name> <signature [...]>  - create LCA database
classify --db <db_name [...]> --query <signature [...]>  - classify genomes
gather <signature> <db_name [...]>                       - classify metagenomes
summarize --db <db_name [...]> --query <signature [...]> - summarize mixture
rankinfo <db_name [...]>                                 - database rank info
compare_csv <csv1> <csv2>                                - compare spreadsheets

** Use '-h' to get subcommand-specific help, e.g.

sourmash lca index -h
'''

def main(sysv_args):
    set_quiet(False)

    commands = {'classify': classify,
                'index': index,
                'summarize': summarize_main,
                'rankinfo': rankinfo_main,
                'gather': gather_main,
                'compare_csv': compare_csv}

    parser = argparse.ArgumentParser(
        description='lowest-common ancestor (LCA) utilities', usage=usage)
    parser.add_argument('lca_command', nargs='?')
    args = parser.parse_args(sysv_args[0:1])

    if not args.lca_command:
        print(usage)
        sys.exit(1)

    if args.lca_command not in commands:
        error('Unrecognized command: {}', args.lca_command)
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.lca_command)
    cmd(sysv_args[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
