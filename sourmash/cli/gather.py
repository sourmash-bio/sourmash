"""search a metagenome signature against dbs"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('gather')
    subparser.add_argument('query', help='query signature')
    subparser.add_argument(
        'databases', nargs='+',
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
        '--traverse-directory', action='store_true',
        help='search all signatures underneath directories'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='save the matched signatures from the database to the '
        'specified file'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='threshold (in bp) for reporting results (default=50,000)'
    )
    subparser.add_argument(
        '--output-unassigned', metavar='FILE',
        help='output unassigned portions of the query as a signature to the '
        'specified file'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=0,
        help='downsample query to the specified scaled factor'
    )
    subparser.add_argument(
        '--ignore-abundance',  action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.gather(args)
