"""
Command-line entry point for 'python -m sourmash.lca'
"""

import sys
import argparse

from . import classify, index, summarize_main, rankinfo_main
from .command_compare_csv import compare_csv
from ..logging import set_quiet, error

usage='''
sourmash lca <command> [<args>] - work with taxonomic information.

** Commands can be:

index <taxonomy.csv> <output_db name> <signature [...]>  - create LCA database
classify --db <db_name [...]> --query <signature [...]>  - classify genomes
summarize --db <db_name [...]> --query <signature [...]> - summarize mixture
rankinfo <db_name [...]>                                 - database rank info
compare_csv <csv1> <csv2>                                - compare spreadsheets

** Use '-h' to get subcommand-specific help, e.g.

sourmash lca index -h
'''

def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
