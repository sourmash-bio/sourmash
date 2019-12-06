def subparser(subparsers):
    subparser = subparsers.add_parser('compare')
    subparser.add_argument('csv1', help='taxonomy spreadsheet output by classify')
    subparser.add_argument('csv2', help='custom taxonomy spreadsheet')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )
    subparser.add_argument(
        '-C', '--start-column', metavar='C', default=2, type=int,
        help='column at which taxonomic assignments start; default=2'
    )
    subparser.add_argument(
        '--tabs', action='store_true',
        help='input spreadsheet is tab-delimited; default is commas'
    )
    subparser.add_argument(
        '--no-headers', action='store_true',
        help='no headers present in taxonomy spreadsheet'
    )
    subparser.add_argument('-f', '--force', action='store_true')
