"'sourmash multigather' - gather many signatures against multiple databases."

from sourmash.cli.utils import add_ksize_arg, add_moltype_args, add_scaled_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('multigather')
    subparser.add_argument(
        '--query', nargs='*', default=[], action='append',
        help='query signature'
    )
    subparser.add_argument(
        '--query-from-file',
        help='file containing list of signature files to query'
    )
    subparser.add_argument(
        '--db', nargs='+', action='append',
        help='signatures/SBTs to search',
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='threshold (in bp) for reporting results (default=50,000)'
    )
    subparser.add_argument(
        '--ignore-abundance',  action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.commands.multigather(args)
