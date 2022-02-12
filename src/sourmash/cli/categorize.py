"'sourmash categorize' - query an SBT for bes match, with many signatures."

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('categorize')
    subparser.add_argument('database', help='location of signature collection/database to load')
    subparser.add_argument(
        'queries', nargs='+',
        help='locations of signatures to categorize'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    add_ksize_arg(subparser, 31)
    subparser.add_argument(
        '--threshold', default=0.08, type=float,
        help='minimum threshold for reporting matches; default=0.08'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    add_moltype_args(subparser)

    # TODO: help messages in these
    subparser.add_argument('--csv', help='output summary CSV to this file')
    subparser.add_argument('--load-csv', default=None)


def main(args):
    import sourmash
    return sourmash.commands.categorize(args)
