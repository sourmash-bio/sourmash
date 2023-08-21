"""rename signature"""

usage="""

### `sourmash signature rename` - rename a signature

Rename the display name for one or more signatures - this is the name
output for matches in `compare`, `search`, `gather`, etc.

For example,

sourmash signature rename file1.sig "new name" -o renamed.sig

will place a renamed copy of the hashes in `file1.sig` in the file
`renamed.sig`. If you provide multiple signatures, all will be renamed
to the same name.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('rename', description=__doc__, usage=usage)
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument('name')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='print debugging output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', 
        help='output renamed signature to this file (default stdout)',
        default='-'
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
    add_pattern_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.rename(args)
