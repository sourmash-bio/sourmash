def subparser(subparsers):
    subparser = subparsers.add_parser('combine')
    subparser.add_argument('sbt_name', help='name to save SBT into')
    subparser.add_argument(
        'sbts', nargs='+',
        help='SBTs to combine to form a new SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5
    )

    subparser = subparsers.add_parser('sbt_combine')
    subparser.add_argument('sbt_name', help='name to save SBT into')
    subparser.add_argument(
        'sbts', nargs='+',
        help='SBTs to combine to form a new SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5
    )


def alt_subparser(subparsers):
    '''Dirty hack to simultaneously support new and previous interface

    This function can be removed with a major version bump.
    '''
    subparser = subparsers.add_parser('sbt_combine')
    subparser.add_argument('sbt_name', help='name to save SBT into')
    subparser.add_argument(
        'sbts', nargs='+',
        help='SBTs to combine to form a new SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5
    )

    subparser = subparsers.add_parser('sbt_combine')
    subparser.add_argument('sbt_name', help='name to save SBT into')
    subparser.add_argument(
        'sbts', nargs='+',
        help='SBTs to combine to form a new SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5
    )


def main(args):
    print(args)
