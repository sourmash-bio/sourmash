"""'sourmash migrate' description goes here"""

def subparser(subparsers):
    subparser = subparsers.add_parser('migrate')
    subparser.add_argument('sbt_name', help='name to save SBT into')


def main(args):
    import sourmash
    return sourmash.commands.migrate(args)
