"""see detailed comparison of signatures"""

usage="""

### `sourmash signature overlap` - detailed comparison of two signatures' overlap

Display a detailed comparison of two signatures. This calculates the
Jaccard similarity (as in `sourmash compare` or `sourmash search`) and
the Jaccard containment in both directions (as with `--containment`).
It also displays the number of hash values in the union and
intersection of the two signatures, as well as the number of disjoint
hash values in each signature.

This command has two uses - first, it is helpful for understanding how
similarity and containment are calculated, and second, it is useful for
analyzing signatures with very small overlaps, where the similarity
and/or containment might be very close to zero.

For example,

sourmash signature overlap file1.sig file2.sig

will display the detailed comparison of `file1.sig` and `file2.sig`.

"""

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('overlap', description=__doc__, usage=usage)
    subparser.add_argument('signature1')
    subparser.add_argument('signature2')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.overlap(args)
