"""summarize mixture"""

from argparse import FileType


def subparser(subparsers):
    subparser = subparsers.add_parser('summarize')
    subparser.add_argument('--db', nargs='+', action='append')
    subparser.add_argument('--query', nargs='+', action='append')
    subparser.add_argument('--threshold', metavar='T', type=int, default=5)
    subparser.add_argument(
        '--traverse-directory', action='store_true',
        help='load all signatures underneath directories'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', type=FileType('wt'),
        help='file to which CSV output will be written'
    )
    subparser.add_argument('--scaled', metavar='FLOAT', type=float)
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
