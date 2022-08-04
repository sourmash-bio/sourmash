"""search taxonomies and output picklists.."""

usage="""

    sourmash tax grep <term> --taxonomy-csv <taxonomy_file> [ ... ]

@@
The 'tax prepare' command reads in one or more taxonomy databases
and saves them into a new database. It can be used to combine databases
in the desired order, as well as output different database formats.

Please see the 'tax prepare' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-prepare-prepare-and-or-combine-taxonomy-files
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('grep', usage=usage)
    subparser.add_argument('search_term')
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
        '-o', '--output', default='-',
        help='output file (defaults to stdout)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in file and taxonomy loading',
    )

def main(args):
    import sourmash
    return sourmash.tax.__main__.grep(args)
