"""output abundance histogram"""

import sys

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('abundhist')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output histogram to this file (in CSV format)'
    )
    subparser.add_argument(
        '--abundances', metavar='FILE',
        help='output hashes and abundances to this file (in CSV format)')
    subparser.add_argument(
        '--md5', default=None,
        help='select signatures whose md5 contains this substring'
    )
    subparser.add_argument(
        '--name', default=None,
        help='select signatures whose name contains this substring'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.abundhist(args)
