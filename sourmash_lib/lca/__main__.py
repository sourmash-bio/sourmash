"""
Command-line entry point for 'python -m sourmash_lib.lca'
"""

import sys
import argparse

from . import classify, index, summarize_main
from ..logging import set_quiet, error

def main(sysv_args):
    set_quiet(False)

    commands = {'classify': classify,
                'index': index,
                'summarize': summarize_main}

    parser = argparse.ArgumentParser(
        description='least-common ancestor (LCA) utilities',
        usage='''sourmash lca <command> [<args>]

Commands can be:

index <taxonomy.csv> <output_db name> <signature [...]>
classify --db <db_name [...]> --query <signature [...]>
summarize --db <db_name [...]> --query <signature [...]>

Use '-h' to get subcommand-specific help, e.g.

sourmash lca index -h
.
''')
    parser.add_argument('lca_command')
    args = parser.parse_args(sysv_args[0:1])
    if args.lca_command not in commands:
        error('Unrecognized command: {}', args.lca_command)
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.lca_command)
    cmd(sysv_args[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
