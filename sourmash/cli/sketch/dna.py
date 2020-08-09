"""create DNA signatures"""

import csv

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('dna')
    subparser.add_argument(
        'filenames', nargs='+', help='file(s) of sequences'
    )


def main(args):
    import sourmash
    return sourmash.command_sketch.dna(args)
