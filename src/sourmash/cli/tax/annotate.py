"""annotate gather results with taxonomy information"""

usage="""

    sourmash tax annotate --gather-csv <gather_csv> [ ... ] --taxonomy-csv <taxonomy_csv> [ ... ]

The 'tax annotate' command reads in gather results CSVs and annotates them
 with taxonomic information.

By default, `tax annotate` produces a gather CSV with an additional `lineage`
 column containing the taxonomic information for each database match.

Please see the 'tax annotate' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-annotate-annotates-gather-output-with-taxonomy
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('annotate',
                                      aliases=['annotate'],
                                      usage=usage)
    subparser.add_argument(
        '-g', '--gather-csv', nargs='*', default = [],
        help='CSV output files from sourmash gather'
    )
    subparser.add_argument(
        '--from-file',  metavar='FILE', default=None,
        help='input many gather results as a text file, with one gather CSV per line'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv', '--taxonomy', metavar='FILE',
        nargs="+", required=True,
        help='database lineages CSV'
    )
    subparser.add_argument(
        '-o', '--output-dir', default= "",
        help='directory for output files'
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
    if not args.gather_csv and not args.from_file:
        raise ValueError(f"No gather CSVs found! Please input via `-g` or `--from-file`.")
    return sourmash.tax.__main__.annotate(args)
