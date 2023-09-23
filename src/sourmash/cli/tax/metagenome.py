"""summarize metagenome gather results"""

usage="""

    sourmash tax metagenome --gather-csv <gather_csv> [ ... ] --taxonomy-csv <taxonomy-csv> [ ... ]

The 'tax metagenome' command reads in metagenome gather result CSVs and
summarizes by taxonomic lineage.

The default output format consists of four columns,
 'query_name,rank,fraction,lineage', where 'fraction' is the fraction
 of the query matched to that reported rank and lineage. The summarization
 is reported for each taxonomic rank.

Alternatively, you can output results at a specific rank (e.g. species)
in 'krona', 'lineage_summary', and 'human' formats.

Use '-F human' to display human-readable output.

Please see the 'tax metagenome' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-metagenome-summarize-metagenome-content-from-gather-results
"""

import sourmash
from sourmash.logging import notify, print_results, error
from sourmash.cli.utils import add_rank_arg, check_rank, check_tax_outputs



def subparser(subparsers):
    subparser = subparsers.add_parser('metagenome',
                                      usage=usage)
    subparser.add_argument(
        '-g', '--gather-csv', action="extend", nargs='*', default = [],
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
        action="extend", nargs='+', required=True,
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
        '-F', '--output-format', default=[], nargs='*', action="extend",
        choices=["human", "csv_summary", "krona", "lineage_summary", "kreport", "lingroup", "bioboxes"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in taxonomy database loading',
    )
    subparser.add_argument(
        '--lins', '--lin-taxonomy', action='store_true', default=False,
        help="use LIN taxonomy in place of standard taxonomic ranks.  Note that the taxonomy CSV must contain 'lin' lineage information."
    )
    subparser.add_argument(
        '--lingroup', '--lingroups', metavar='FILE', default=None,
        help="CSV containing 'name', 'lin' columns, where 'lin' is the lingroup prefix. Will produce a 'lingroup' report containing taxonomic summarization for each group."
    )
    add_rank_arg(subparser)

def main(args):
    import sourmash
    try:
        if not args.gather_csv and not args.from_file:
            raise ValueError(f"No gather CSVs found! Please input via '-g' or '--from-file'.")
        if args.rank:
            args.rank = check_rank(args)
        args.output_format = check_tax_outputs(args, rank_required = ['krona', 'lineage_summary'], incompatible_with_lins = ['bioboxes', 'kreport'], use_lingroup_format=True)

    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        import sys; sys.exit(-1)

    return sourmash.tax.__main__.metagenome(args)
