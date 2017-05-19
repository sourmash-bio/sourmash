"""
sourmash command line.
"""
from __future__ import print_function
import sys
import argparse

from .logging import notify, error

from .commands import (categorize, compare, compute, convert, dump, import_csv,
                       sbt_gather, sbt_index, sbt_combine, sbt_search, search,
                       plot, watch)


def main():
    commands = {'search': search, 'compute': compute,
                'compare': compare, 'plot': plot,
                'import_csv': import_csv, 'dump': dump,
                'sbt_index': sbt_index, 'sbt_search': sbt_search,
                'categorize': categorize, 'sbt_gather': sbt_gather,
                'watch': watch, 'convert': convert,
                'sbt_combine': sbt_combine}
    parser = argparse.ArgumentParser(
        description='work with RNAseq signatures',
        usage='''sourmash <command> [<args>]

Commands can be:

compute <filenames>         Compute MinHash signatures for sequences in files.
compare <filenames.sig>     Compute similarity matrix for multiple signatures.
search <query> <against>    Search a signature against a list of signatures.
plot <matrix>               Plot a distance matrix made by 'compare'.

Sequence Bloom Tree (SBT) utilities:

sbt_index                   Index a collection of signatures with an SBT.
sbt_combine                 Combine multiple SBTs into a new one.
sbt_search                  Search a signature against a SBT.
categorize                  Identify best matches for many signatures using
                              an SBT.
sbt_gather                  Search a metagenome signature for multiple
                              non-overlapping matches in the SBT.
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
