"""'sourmash import_csv' description goes here"""

from argparse import FileType
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('import_csv')
    subparser.add_argument('mash_csvfile', help='CSV file with mash sketches')
    subparser.add_argument(
        '-o', '--output', type=FileType('wt'),
        default=sys.stdout,
        help='save signature generated from data here'
    )


def main(args):
    import sourmash
    return sourmash.commands.import_csv(args)
