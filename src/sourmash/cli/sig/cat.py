"""concatenate signature files"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('cat')
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '-u', '--unique', action='store_true',
        help='keep only distinct signatures, removing duplicates (based on md5sum)'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.cat(args)
