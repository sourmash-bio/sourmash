"""classify metagenomes"""


def subparser(subparsers):
    subparser = subparsers.add_parser('gather')
    subparser.add_argument('query')
    subparser.add_argument('db', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output')
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--output-unassigned', metavar='FILE',
        help='output unassigned portions of the query as a signature to this '
        'file'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances if present'
    )


def main(args):
    import sourmash
    return sourmash.lca.gather_main(args)
