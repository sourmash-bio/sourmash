"""create a manifest for a collection of signatures"""

usage="""

    sourmash sig manifest <filename> -o manifest.csv

This will output a sourmash manifest in CSV format. This manifest
can be used as a picklist with --picklist manifest.csv::manifest.

The manifest will be rebuilt by iterating over the signatures in the
file unless --no-rebuild-manifest is specified; for large
collections, rebuilding the manifest can take a long time!

See also the 'describe' and 'fileinfo' commands under 'sourmash sig'.

"""


def subparser(subparsers):
    subparser = subparsers.add_parser('manifest', usage=usage)
    subparser.add_argument('location')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debug information'
    )
    subparser.add_argument(
        '-o', '--output', '--csv', metavar='FILE',
        help='output information to a CSV file',
        required=True,
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--no-rebuild-manifest', help='use existing manifest if available',
        action='store_true'
    )
    subparser.add_argument(
        '-F', '--manifest-format',
        help="format of manifest output file; default is 'csv')",
        default='csv',
        choices=['csv', 'sql'],
    )

def main(args):
    import sourmash
    return sourmash.sig.__main__.manifest(args)
