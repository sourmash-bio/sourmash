def subparser(subparsers):
    subparser = subparsers.add_parser('index')
    subparser.add_argument('csv', help='taxonomy spreadsheet')
    subparser.add_argument('lca_db_out', help='output database name')
    subparser.add_argument('signatures', nargs='+', help='one or more sourmash signatures')
