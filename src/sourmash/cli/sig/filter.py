"""filter k-mers on abundance"""

usage="""

### `sourmash signature filter` - remove hashes based on abundance

Filter the hashes in the specified signature(s) by abundance, by either
`-m/--min-abundance` or `-M/--max-abundance` or both. Abundance selection is
inclusive, so `-m 2 -M 5` will select hashes with abundance greater than
or equal to 2, and less than or equal to 5.

For example,

sourmash signature -m 2 *.sig


will output new signatures containing only hashes that occur two or
more times in each signature.

The `filter` command accepts the same selectors as `extract`.

"""

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('filter', description=__doc__, usage=usage)
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output signature to this file (default stdout)',
        default='-'
    )
    subparser.add_argument(
        '--md5', type=str, default=None,
        help='select signatures whose md5 contains this substring'
    )
    subparser.add_argument(
        '--name', type=str, default=None,
        help='select signatures whose name contains this substring'
    )
    subparser.add_argument(
        '-m', '--min-abundance', type=int, default=1,
        help='keep hashes >= this minimum abundance'
    )
    subparser.add_argument(
        '-M', '--max-abundance', type=int, default=None,
        help='keep hashes <= this maximum abundance'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.filter(args)
