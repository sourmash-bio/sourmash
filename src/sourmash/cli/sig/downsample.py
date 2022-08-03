"""downsample one or more signatures"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_num_arg)


def subparser(subparsers):
    subparser = subparsers.add_parser('downsample')
    subparser.add_argument('signatures', nargs="*")
    subparser.add_argument(
        '--scaled', type=int, default=0,
        help='scaled value to downsample to'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output signature to this file (default stdout)',
        default='-',
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_num_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.sig.__main__.downsample(args)
