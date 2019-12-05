def subparser(subparsers):
    subparser = subparsers.add_parser('rename')
    subparser.add_argument('signatures', nargs='+', help='list of signatures')
    subparser.add_argument('-o', '--output', metavar='OUT')
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances if present'
    )
