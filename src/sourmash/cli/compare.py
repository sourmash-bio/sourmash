"""compare sequence signatures made by compute"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args)


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
        '--max-containment', action='store_true',
        help='calculate max containment instead of similarity'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
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
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.compare(args)
