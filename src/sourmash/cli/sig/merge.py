"""merge one or more signatures"""

usage="""

### `sourmash signature merge` - merge two or more signatures into one

Merge two (or more) signatures.

For example,

sourmash signature merge file1.sig file2.sig -o merged.sig

will output the union of all the hashes in `file1.sig` and `file2.sig`
to `merged.sig`.

All of the signatures passed to merge must either have been created
with `-p abund`, or not.  If they have `track_abundance` on,
then the merged signature will have the sum of all abundances across
the individual signatures.  The `--flatten` flag will override this
behavior and allow merging of mixtures by removing all abundances.

Note: `merge` only creates one output file, with one signature in it,
in the JSON `.sig` format.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('merge', description=__doc__, usage=usage)
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
        '--flatten', action='store_true',
        help='remove abundances from all signatures'
    )
    subparser.add_argument(
        '--name',
        help='rename merged signature'
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
    return sourmash.sig.__main__.merge(args)
