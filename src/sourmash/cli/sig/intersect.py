"""intersect two or more signatures"""

usage="""

### `sourmash signature intersect` - intersect two (or more) signatures

Output the intersection of the hash values in multiple signature files.

For example,


sourmash signature intersect file1.sig file2.sig file3.sig -o intersect.sig

will output the intersection of all the hashes in those three files to
`intersect.sig`.

The `intersect` command flattens all signatures, i.e. the abundances
in any signatures will be ignored and the output signature will have
`track_abundance` turned off. See `sourmash signature flatten` for more details.

Note: `intersect` only creates one output file, with one signature in it.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('intersect', description=__doc__, usage=usage)
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '-A', '--abundances-from', metavar='FILE',
        help='intersect with & take abundances from this signature'
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


def main(args):
    import sourmash
    return sourmash.sig.__main__.intersect(args)
