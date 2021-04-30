"""merge one or more signatures"""

import sys

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('merge')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '--flatten', action='store_true',
        help='remove abundances from all signatures'
    )
    subparser.add_argument(
        '--name',
        help='rename merged signature'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.merge(args)
