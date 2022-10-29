"""show details of signature"""

usage="""

### `sourmash signature describe` - display detailed information about signatures

Display signature details.

For example,

sourmash sig describe tests/test-data/47.fa.sig

will display:

signature filename: tests/test-data/47.fa.sig
signature: NC_009665.1 Shewanella baltica OS185, complete genome
source file: 47.fa
md5: 09a08691ce52952152f0e866a59f6261
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 5177
signature license: CC0

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('describe', description=__doc__, usage=usage)
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
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_pattern_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.describe(args)
