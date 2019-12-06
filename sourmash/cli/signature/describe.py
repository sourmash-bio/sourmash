from argparse import FileType

def subparser(subparsers):
    subparser = subparsers.add_parser('describe')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--csv', metavar='FILE', type=FileType('wt'),
        help='output information to a CSV file'
    )
