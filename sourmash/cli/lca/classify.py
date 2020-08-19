"""classify genomes"""


def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument('--db', nargs='+', action='append',
                           help='databases to use to classify')
    subparser.add_argument('--query', nargs='*', default=[], action='append',
                           help='query signatures to classify')
    subparser.add_argument('--query-from-file',
                           help='file containing list of signature files to query')
    subparser.add_argument('--threshold', metavar='T', type=int, default=5)
    subparser.add_argument(
        '--majority', action='store_true',
        help='use majority vote classification instead of lca'
    )
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


def main(args):
    import sourmash
    return sourmash.lca.command_classify.classify(args)
