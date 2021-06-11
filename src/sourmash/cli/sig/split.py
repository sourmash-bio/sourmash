"""concatenate signature files"""

import sourmash


def subparser(subparsers):
    "Aggregates argparse for subtract subcommands"
    subparser = subparsers.add_parser('split')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--outdir', help='output signatures to this directory'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.split(args)
