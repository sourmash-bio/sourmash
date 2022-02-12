"""filter k-mers on abundance"""

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('filter')
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
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.filter(args)
