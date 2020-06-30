"""compare sequence signatures made by compute"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('compare')
    subparser.add_argument(
        'signatures', nargs='*', help='list of signatures to compare',
        default=[]
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true', help='suppress non-error output'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    subparser.add_argument(
        '-o', '--output', metavar='F',
        help='file to which output will be written; default is terminal '
        '(standard output)'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances even if present'
    )
    subparser.add_argument(
        '--containment', action='store_true',
        help='calculate containment instead of similarity'
    )
    subparser.add_argument(
        '--traverse-directory', action='store_true',
        help='compare all signatures underneath directories'
    )
    subparser.add_argument(
        '--from-file',
        help='a file containing a list of signatures file to compare'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='continue past errors in file loading'
    )
    subparser.add_argument(
        '--csv', metavar='F',
        help='write matrix to specified file in CSV format (with column '
        'headers)'
    )
    subparser.add_argument(
        '-p', '--processes', metavar='N', type=int, default=None,
        help='Number of processes to use to calculate similarity')


def main(args):
    import sourmash
    return sourmash.commands.compare(args)
