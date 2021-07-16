"""subtract one or more signatures"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg)


def subparser(subparsers):
    subparser = subparsers.add_parser('subtract')
    subparser.add_argument('signature_from')
    subparser.add_argument('subtraction_sigs', nargs='+')
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
        help='remove abundance from signatures before subtracting'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.subtract(args)
