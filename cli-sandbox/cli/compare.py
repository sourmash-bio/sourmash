from argparse import FileType
from . import add_ksize_arg, add_moltype_args, DEFAULT_LOAD_K

def subparser(subparsers):
    subparser = subparsers.add_parser('compare')
    subparser.add_argument(
        'signatures', nargs='+', help='list of signatures to compare'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='file to which output will be written; default is terminal '
        '(standard output)'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances even if present'
    )
    add_ksize_arg(subparser, DEFAULT_LOAD_K)
    add_moltype_args(subparser)
    parser.add_argument(
        '--traverse-directory', action='store_true',
        help='compare all signatures underneath directories'
    )
    parser.add_argument(
        '--csv', type=FileType('w'),
        help='save matrix in CSV format (with column headers)'
    )
    parser.add_argument(
        '-p', '--processes', type=int, default=None,
        help='Number of processes to use to calculate similarity')
    parser.add_argument(
        '-q', '--quiet', action='store_true', help='suppress non-error output'
    )
