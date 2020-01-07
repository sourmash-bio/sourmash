"""'sourmash dump' description goes here"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('dump')
    subparser.add_argument('filenames', nargs='+')
    add_ksize_arg(subparser, 31)


def main(args):
    import sourmash
    return sourmash.commands.dump(args)
