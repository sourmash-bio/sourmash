"""'sourmash import_csv' description goes here"""

import sys
from sourmash.logging import notify


def subparser(subparsers):
    subparser = subparsers.add_parser('import_csv')
    subparser.add_argument('mash_csvfile', help='CSV file with mash sketches')
    subparser.add_argument(
        '-o', '--output',
        help='save signature generated from data to this file (default stdout)'
    )


def main(args):
    import sourmash
    notify("** WARNING: 'import_csv' is deprecated as of sourmash 4.0, and will")
    notify("**    be removed in sourmash 5.0; use 'sourmash sig import' instead.")
    notify('')
    return sourmash.commands.import_csv(args)
