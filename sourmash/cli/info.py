def subparser(subparsers):
    subparser = subparsers.add_parser('info')
    subparser.add_argument(
        '--verbose', action='store_true',
        help='report versions of khmer and screed'
    )
