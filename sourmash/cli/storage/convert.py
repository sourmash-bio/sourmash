"""'sourmash storage convert' description goes here"""

def subparser(subparsers):
    subparser = subparsers.add_parser('convert')
    subparser.add_argument(
        'sbt', help='name to save SBT into'
    )
    subparser.add_argument(
        '-b', '--backend', type=str,
        help='Backend to convert to'
    )


def main(args):
    import sourmash
    return sourmash.sbt.convert_cmd(args.sbt, args.backend)
