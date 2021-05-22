"""classify genomes"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument('gather_results', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    tax_group = subparser.add_mutually_exclusive_group(required=False)
    tax_group.add_argument(
        '-t', '--taxonomy', default='gtdb',
        choices = ['gtdb', 'ncbi'],
        help='Use an included taxonomy (default gtdb)'
    )
    tax_group.add_argument(
        '-u', '--user-taxonomy', metavar='FILE',
        help='Instead, input your own taxonomy file (see docs for formatting instructions)'
    )
    subparser.add_argument(
        '-r', '--rank',
        help='Summarize genome taxonomy at this rank and above'
    )


def main(args):
    import sourmash
    return sourmash.tax.__main__.summarize(args)
