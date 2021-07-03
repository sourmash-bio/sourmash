"""combine multiple taxonomy databases into one."""

usage="""

    sourmash tax combine_tax --taxonomy-csv [taxonomy-csv(s)] -o <output>

The 'tax combine_tax' command reads in one or more taxonomy databases
and saves them into a new database.

Please see the 'tax combine_tax' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-annotate-annotates-gather-output-with-taxonomy

@CTB fix link
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('combine_tax',
                                      usage=usage)
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv', metavar='FILE',
        nargs="+", required=True,
        help='database lineages'
    )
    subparser.add_argument(
        '-o', '--output', required=True,
        help='output file',
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
    return sourmash.tax.__main__.combine_tax(args)
