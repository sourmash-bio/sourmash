"""database rank info"""

def subparser(subparsers):
    subparser = subparsers.add_parser('rankinfo')
    subparser.add_argument('db', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )
    subparser.add_argument('--scaled', metavar='FLOAT', type=float)
    subparser.add_argument(
        '--minimum-num', type=int, default=0,
        help='Minimum number of different lineages a k-mer must be in to be counted'
    )


def main(args):
    import sourmash
    return sourmash.lca.command_rankinfo.rankinfo_main(args)
