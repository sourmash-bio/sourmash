"""summarize metagenome gather results at rank"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('summarize')
    subparser.add_argument('gather_results', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '-r', '--rank',
        help='Summarize metagenome gather results to this rank and above'
    )


def main(args):
    import sourmash
    return sourmash.tax.__main__.summarize(args)
