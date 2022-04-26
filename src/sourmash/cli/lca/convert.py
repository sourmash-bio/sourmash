"""convert (and/or combine) LCA databases into other formats"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('convert')
    subparser.add_argument('--db', nargs='+', action='append',
                           help='databases to combine and convert')
    subparser.add_argument('-o', '--output-database',
                           help="name of output database file")
    subparser.add_argument(
        '--scaled', metavar='S', default=10000, type=float
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debugging output'
    )
    #subparser.add_argument('-f', '--force', action='store_true')
    subparser.add_argument(
        '-F', '--database-format',
        help="format of output database; default is 'json')",
        default='json',
        choices=['json', 'sql'],
    )

    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.lca.command_convert.convert_main(args)
