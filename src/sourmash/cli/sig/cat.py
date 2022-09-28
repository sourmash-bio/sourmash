"""concatenate signature files"""

usage="""

### `sourmash signature cat` - concatenate multiple signatures together

Concatenate signature files.

For example,

sourmash signature cat file1.sig file2.sig -o all.sig

will combine all signatures in `file1.sig` and `file2.sig` and put them
in the file `all.sig`.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    # working on this
    subparser = subparsers.add_parser('cat', description=__doc__, usage=usage)
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '-u', '--unique', action='store_true',
        help='keep only distinct signatures, removing duplicates (based on md5sum)'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_pattern_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.cat(args)
