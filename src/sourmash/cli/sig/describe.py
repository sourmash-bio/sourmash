"""show details of signature"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('describe')
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='provide debugging output'
    )
    subparser.add_argument(
        '--csv', metavar='FILE',
        help='output information to a CSV file'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_pattern_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.describe(args)
