"""index signatures for rapid search"""

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('index')
    subparser.add_argument('sbt_name', help='name to save SBT into')
    subparser.add_argument(
        'signatures', nargs='+',
        help='signatures to load into SBT'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    add_ksize_arg(subparser, 31)
    subparser.add_argument(
        '-d', '--n_children', metavar='D', type=int, default=2,
        help='number of children for internal nodes; default=2'
    )
    subparser.add_argument(
        '--traverse-directory', action='store_true',
        help='load all signatures underneath any directories'
    )
    subparser.add_argument(
        '--append', action='store_true', default=False,
        help='add signatures to an existing SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5,
        help='Bloom filter size used for internal nodes'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try loading all files with --traverse-directory'
    )
    subparser.add_argument(
        '-s', '--sparseness', metavar='FLOAT', type=float, default=.0,
        help='What percentage of internal nodes will not be saved; ranges '
        'from 0.0 (save all nodes) to 1.0 (no nodes saved)'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=0,
        help='downsample signatures to the specified scaled factor'
    )
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.index(args)
