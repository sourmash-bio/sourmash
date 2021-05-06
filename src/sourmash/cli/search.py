"""search a signature against other signatures"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('search')
    subparser.add_argument(
        'query', help='query signature'
    )
    subparser.add_argument(
        'databases', nargs='+',
        help='signatures/SBTs to search',
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--threshold', metavar='T', default=0.08, type=float,
        help='minimum threshold for reporting matches; default=0.08'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='output matching signatures to the specified file'
    )
    subparser.add_argument(
        '--best-only', action='store_true',
        help='report only the best match (with greater speed)'
    )
    subparser.add_argument(
        '-n', '--num-results', default=3, type=int, metavar='N',
        help='number of results to report'
    )
    subparser.add_argument(
        '--containment', action='store_true',
        help='score based on containment rather than similarity'
    )
    subparser.add_argument(
        '--max-containment', action='store_true',
        help='score based on max containment rather than similarity'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances if present; note: has no effect if '
        '--containment or --max-containment is specified'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=0,
        help='downsample query to this scaled factor (yields greater speed)'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.search(args)
