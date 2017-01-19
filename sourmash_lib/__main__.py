"""
sourmash command line.
"""
from __future__ import print_function
import sys
import argparse

from .logging import notify, error

from .commands import (categorize, compare, compute, convert, dump, import_csv,
                       sbt_gather, sbt_index, sbt_search, search, plot, watch)


class SourmashCommands(object):
    def __init__(self):
        self.commands = {'search': search, 'compute': compute,
                         'compare': compare, 'plot': plot,
                         'import_csv': import_csv, 'dump': dump,
                         'sbt_index': sbt_index, 'sbt_search': sbt_search,
                         'categorize': categorize, 'sbt_gather': sbt_gather,
                         'watch': watch, 'convert': convert}
        parser = argparse.ArgumentParser(
            description='work with RNAseq signatures',
            usage='''sourmash <command> [<args>]

Commands can be:

   compute <filenames>         Compute signatures for sequences in these files.
   compare <filenames.sig>     Compute distance matrix for given signatures.
   search <query> <against>    Search for matching signatures.
   plot <matrix>               Plot a distance matrix made by 'compare'.

   import_csv                  Import signatures from a CSV file.
   convert                     Convert signatures from YAML to JSON.

   sbt_index                   Index signatures with a Sequence Bloom Tree.
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
        if args.command not in self.commands:
            error('Unrecognized command')
            parser.print_help()
            sys.exit(1)

        cmd = self.commands.get(args.command)
        notify('# running sourmash subcommand: %s' % args.command,
               file=sys.stderr)
        cmd(sys.argv[2:])


def main():
    SourmashCommands()
    return 0
