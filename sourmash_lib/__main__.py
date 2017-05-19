"""
sourmash command line.
"""
from __future__ import print_function
import sys
import argparse

from .logging import notify, error

from .commands import (categorize, compare, compute, dump, import_csv,
                       sbt_gather, sbt_index, sbt_combine, search,
                       plot, watch)


def main():
    commands = {'search': search, 'compute': compute,
                'compare': compare, 'plot': plot,
                'import_csv': import_csv, 'dump': dump,
                'sbt_index': sbt_index,
                'categorize': categorize, 'sbt_gather': sbt_gather,
                'watch': watch,
                'sbt_combine': sbt_combine}
    parser = argparse.ArgumentParser(
        description='work with RNAseq signatures',
        usage='''sourmash <command> [<args>]

Commands can be:

compute <filenames>         Compute signatures for sequences in these files.
compare <filenames.sig>     Compute distance matrix for given signatures.
search <query> <against>    Search for matching signatures.
plot <matrix>               Plot a distance matrix made by 'compare'.

import_csv                  Import signatures from a CSV file.

sbt_index                   Index signatures with a Sequence Bloom Tree.
sbt_combine                 Combine multiple Sequence Bloom Trees into a new one.
sbt_search                  Search a Sequence Bloom Tree.
categorize                  Categorize signatures with a SBT.
sbt_gather                  Search a signature for multiple matches.
watch                       Classify a stream of sequences using a SBT.

Use '-h' to get subcommand-specific help, e.g.

sourmash compute -h
.
''')
    parser.add_argument('command')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands:
        error('Unrecognized command')
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.command)
    notify('# running sourmash subcommand: %s' % args.command,
           file=sys.stderr)
    cmd(sys.argv[2:])
