"""summarize taxonomy/lineage information"""

usage = """

    sourmash tax summarize <taxonomy_file> [ <more files> ... ]

The 'tax summarize' command reads in one or more taxonomy databases
or lineage files (produced by 'tax annotate'), combines them,
and produces a human readable summary.

Please see the 'tax summarize' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#command-line.html#sourmash-tax-summarize-print-summary-information-for-lineage-spreadsheets-or-taxonomy-databases

"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser("summarize", usage=usage)
    subparser.add_argument(
        "-q", "--quiet", action="store_true", help="suppress non-error output"
    )
    subparser.add_argument(
        "taxonomy_files",
        metavar="FILE",
        nargs="+",
        action="extend",
        help="database lineages",
    )
    subparser.add_argument(
        "-o",
        "--output-lineage-information",
        help="output a CSV file containing individual lineage counts",
    )
    subparser.add_argument(
        "--keep-full-identifiers",
        action="store_true",
        help="do not split identifiers on whitespace",
    )
    subparser.add_argument(
        "--keep-identifier-versions",
        action="store_true",
        help="after splitting identifiers, do not remove accession versions",
    )
    subparser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="continue past errors in file and taxonomy loading",
    )
    subparser.add_argument(
        "--lins",
        "--lin-taxonomy",
        action="store_true",
        default=False,
        help="use LIN taxonomy in place of standard taxonomic ranks.",
    )
    subparser.add_argument(
        '--ictv', '--ictv-taxonomy', action='store_true', default=False,
        help="use ICTV taxonomy in place of standard taxonomic ranks.  Note that the taxonomy CSV must contain ICTV ranks."
    )


def main(args):
    return sourmash.tax.__main__.summarize(args)
