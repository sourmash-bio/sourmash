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
