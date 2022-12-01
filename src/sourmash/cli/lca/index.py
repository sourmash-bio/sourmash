"""create LCA database"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('index')
    subparser.add_argument('csv', help='taxonomy spreadsheet')
    subparser.add_argument('lca_db_out', help='output database name')
    subparser.add_argument(
        'signatures', nargs='*',
        help='signatures or directory of signatures to index (optional if provided via --from-file)'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
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
    subparser.add_argument(
        '-C', '--start-column', metavar='C', default=2, type=int,
        help='column at which taxonomic assignments start; default=2'
    )
    subparser.add_argument(
        '--tabs', action='store_true',
        help='input spreadsheet is tab-delimited; default is commas'
    )
    subparser.add_argument(
        '--no-headers', action='store_true',
        help='no headers present in taxonomy spreadsheet'
    )
    subparser.add_argument(
        '--split-identifiers', action='store_true',
        help='split names in signatures on whitespace'
    )
    subparser.add_argument(
        '--keep-identifier-versions', action='store_true',
        help='do not remove accession versions'
    )
    subparser.add_argument('-f', '--force', action='store_true')
    subparser.add_argument(
        '--report', help='output a report on anomalies, if any'
    )
    subparser.add_argument(
        '--require-taxonomy', action='store_true',
        help='ignore signatures with no taxonomy entry'
    )
    subparser.add_argument(
        '--fail-on-missing-taxonomy', action='store_true',
        help='fail quickly if taxonomy is not available for an identifier',
    )
    subparser.add_argument(
        '-F', '--database-format',
        help="format of output database; default is 'json')",
        default='json',
        choices=['json', 'sql'],
    )

    add_ksize_arg(subparser, default=31)
    add_moltype_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.lca.command_index.index(args)
