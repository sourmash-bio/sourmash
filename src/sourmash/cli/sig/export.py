"""export a signature, e.g. to mash"""

usage="""

### `sourmash signature export` - export signatures to mash.

Export signatures from sourmash format. Currently only supports
mash dump format.

For example,

sourmash signature export filename.sig -o filename.sig.msh.json

"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('export', description=__doc__, usage=usage)
    subparser.add_argument('filename')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.export(args)
