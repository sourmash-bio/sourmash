"""classify genomes"""


def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument('--db', nargs='+', action='append')
    subparser.add_argument('--query', nargs='+', action='append')
    subparser.add_argument('--threshold', metavar='T', type=int, default=5)
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output CSV to the specified file; by default output to stdout'
    )
    subparser.add_argument('--scaled', type=float)
    subparser.add_argument(
        '--traverse-directory', action='store_true',
        help='load all signatures underneath directories'
    )


def main(args):
    import sourmash
    return sourmash.lca.command_classify.classify(args)
