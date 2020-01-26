"""'sourmash import_csv' description goes here"""

import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('import_csv')
    subparser.add_argument('mash_csvfile', help='CSV file with mash sketches')
    subparser.add_argument(
        '-o', '--output',
        help='save signature generated from data to this file (default stdout)'
    )


def main(args):
    import sourmash
    return sourmash.commands.import_csv(args)
