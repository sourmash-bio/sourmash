"""summarize metagenome gather results"""

usage="""

    sourmash tax metagenome --gather-csv <gather_csv> [ ... ] --taxonomy-csv <taxonomy-csv> [ ... ]

The 'tax metagenome' command reads in metagenome gather result CSVs and
summarizes by taxonomic lineage.

The default output format consists of four columns,
 `query_name,rank,fraction,lineage`, where `fraction` is the fraction
 of the query matched to that reported rank and lineage. The summarization
 is reported for each taxonomic rank.

Alternatively, you can output results at a specific rank (e.g. species)
in `krona` or `lineage_summary` formats.

Please see the 'tax metagenome' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-metagenome-summarize-metagenome-content-from-gather-results
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('metagenome',
                                      aliases=['summarize'],
                                      usage=usage)
    subparser.add_argument(
        '-g', '--gather-csv', nargs='*', default = [],
        help='CSVs from sourmash gather'
    )
    subparser.add_argument(
        '--from-file',  metavar='FILE', default = None,
        help='input many gather results as a text file, with one gather CSV per line'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output-base', default='-',
        help='base filepath for output file(s) (default stdout)'
    )
    subparser.add_argument(
        '--output-dir', default= "",
        help='directory for output files'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv', '--taxonomy', metavar='FILE',
        nargs='+', required=True,
        help='database lineages CSV'
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
        '--output-format', default=['csv_summary'], nargs='+', choices=["csv_summary", "krona", "lineage_summary"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-r', '--rank', choices=['strain','species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'],
        help='For non-default output formats: Summarize genome taxonomy at this rank and above. Note that the taxonomy CSV must contain lineage information at this rank.'
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in taxonomy database loading',
    )

def main(args):
    import sourmash
    if not args.gather_csv and not args.from_file:
        raise ValueError(f"No gather CSVs found! Please input via `-g` or `--from-file`.")
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    if not args.rank:
        if any(x in ["krona", "lineage_summary"] for x in args.output_format):
            raise ValueError(f"Rank (--rank) is required for krona and lineage_summary output formats.")
    return sourmash.tax.__main__.metagenome(args)
