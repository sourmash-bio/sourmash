"""summarize mixture"""


def subparser(subparsers):
    subparser = subparsers.add_parser('summarize')
    subparser.add_argument('--db', nargs='+', action='append',
                           help='one or more LCA databases to use')
    subparser.add_argument('--query', nargs='*', default=[], action='append',
                           help='one or more signature files to use as queries')
    subparser.add_argument('--query-from-file',
                           help='file containing list of signature files to query')
    subparser.add_argument('--threshold', metavar='T', type=int, default=5,
                           help='minimum number of hashes to require for a match')
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='file to which CSV output will be written'
    )
    subparser.add_argument('--scaled', metavar='FLOAT', type=float,
                           help='scaled value to downsample to')

    subparser.add_argument('--ignore-abundance', action='store_true',
                           help='ignore hash abundances in query signatures do not weight results')

    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )


def main(args):
    import sourmash
    return sourmash.lca.command_summarize.summarize_main(args)
