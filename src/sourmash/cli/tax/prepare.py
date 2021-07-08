"""combine multiple taxonomy databases into one."""

usage="""

    sourmash tax prepare --taxonomy-csv <taxonomy_file> [ ... ] -o <output>

The 'tax prepare' command reads in one or more taxonomy databases
and saves them into a new database. It can be used to combine databases
in the desired order, as well as output different database formats.

Please see the 'tax prepare' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-prepare-prepare-and-or-combine-taxonomy-files
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('prepare',
                                      usage=usage)
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv', '--taxonomy', metavar='FILE',
        nargs="+", required=True,
        help='database lineages'
    )
    subparser.add_argument(
        '-o', '--output', required=True,
        help='output file',
    )
    subparser.add_argument(
        '-F', '--database-format',
        help="format of output file; default is 'sql')",
        default='sql',
        choices=['csv', 'sql'],
    )
    subparser.add_argument(
        '--keep-full-identifiers', action='store_true',
        help='do not split identifiers on whitespace'
    )
    subparser.add_argument(
        '--keep-identifier-versions', action='store_true',
        help='after splitting identifiers, do not remove accession versions'
    )
    subparser.add_argument(
        '--fail-on-missing-taxonomy', action='store_true',
        help='fail quickly if taxonomy is not available for an identifier',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in file and taxonomy loading',
    )

def main(args):
    import sourmash
    return sourmash.tax.__main__.prepare(args)
