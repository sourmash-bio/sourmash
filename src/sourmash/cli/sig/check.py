"""check signature collections against a picklist"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('check')
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output picklist with remaining unmatched entries to this file',
        default='-',
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '--save-manifest-matching',
        help='save a manifest of the matching entries to this file.'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_pattern_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.check(args)
