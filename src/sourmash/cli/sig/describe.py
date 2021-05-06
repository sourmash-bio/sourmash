"""show details of signature"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('describe')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--csv', metavar='FILE',
        help='output information to a CSV file'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.describe(args)
