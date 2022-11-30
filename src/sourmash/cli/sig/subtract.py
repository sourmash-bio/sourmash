"""subtract one or more signatures"""

usage="""

### `sourmash signature subtract` - subtract other signatures from a signature

Subtract all of the hash values from one signature that are in one or more
of the others.

For example,


sourmash signature subtract file1.sig file2.sig file3.sig -o subtracted.sig

will subtract all of the hashes in `file2.sig` and `file3.sig` from
`file1.sig`, and save the new signature to `subtracted.sig`.

To use `subtract` on signatures calculated with
`-p abund`, you must specify `--flatten`.

Note: `subtract` only creates one output file, with one signature in it.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg)


def subparser(subparsers):
    subparser = subparsers.add_parser('subtract', description=__doc__, usage=usage)
    subparser.add_argument('signature_from')
    subparser.add_argument('subtraction_sigs', nargs='+')
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
        help='remove abundance from signatures before subtracting'
    )
    subparser.add_argument(
        '-A', '--abundances-from', metavar='FILE',
        help='intersect with & take abundances from this signature'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.subtract(args)
