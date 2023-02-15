"""classify genomes from gather results"""

usage="""

    sourmash tax genome --gather-csv <gather_csv> [ ... ] --taxonomy-csv <taxonomy-csv> [ ... ]

The 'tax genome' command reads in genome gather result CSVs and reports likely
classification for each query genome.

By default, classification uses a containment threshold of 0.1,
meaning at least 10 percent of the query was covered by matches with
the reported taxonomic rank and lineage.  You can specify an alternate
classification threshold or force classification by taxonomic rank
instead, e.g. at species or genus-level.

The default output format consists of five columns,
 'query_name,status,rank,fraction,lineage', where 'fraction' is the fraction
 of the query matched to the reported rank and lineage. The 'status' column
 provides additional information on the classification, and can be:
  - 'match' - this query was classified
  - 'nomatch'- this query could not be classified
  - 'below_threshold' - this query was classified at the specified rank,
     but the query fraction matched was below the containment threshold

Use '-F human' to display human-readable output instead.

Optionally, you can report classifications in 'krona' format, but note
that this forces classification by rank, rather than containment threshold.

Please see the 'tax genome' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-genome-classify-a-genome-using-gather-results
"""

import argparse
import sourmash
from sourmash.logging import notify, print_results, error
from sourmash.cli.utils import add_tax_threshold_arg

def subparser(subparsers):
    subparser = subparsers.add_parser('genome',
                                      aliases=['classify'],
                                      usage=usage)
    subparser.add_argument(
        '-g', '--gather-csv', action='extend', nargs='*', default = [],
        help='CSVs output by sourmash gather for this sample'
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
        nargs='*', required=True, action='extend',
        help='database lineages CSV'
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
        '-r', '--rank', choices=['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'],
        help='Summarize genome taxonomy at this rank and above. Note that the taxonomy CSV must contain lineage information at this rank.'
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
        '-F', '--output-format', default=[], nargs='*', action='extend',
        choices=["csv_summary", "krona", "human", "lineage_csv", "LINgroup_report"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past survivable errors in loading taxonomy database or gather results',
    )
    subparser.add_argument(
        '--LIN-taxonomy', action='store_true', default=False,
        help='use LIN taxonomy in place of standard taxonomic ranks.  Note that the taxonomy CSV must contain LIN lineage information.'
    )
    subparser.add_argument(
        '--LIN-position', type=int, default=None,
        help='For non-default output formats: summarize taxonomy at this LIN position and above. Replaces "--rank" for standard taxonomy. Note that the taxonomy CSV must contain LIN with information at this position.'
    )
    subparser.add_argument(
        '--LINgroups', metavar='FILE', default=None,
        help='CSV containing LINgroup_name, LINgroup_prefix. Will produce a "LINgroup_report" file containing taxonomic summarization for each LINgroup.'
    )
    add_tax_threshold_arg(subparser, 0.1)


def main(args):
    import sourmash
    try:
        if not args.gather_csv and not args.from_file:
            raise ValueError(f"No gather CSVs found! Please input via '-g' or '--from-file'.")
        # handle LIN options
        if args.LIN_taxonomy:
            if args.LIN_position:
                args.rank = args.LIN_position
            if args.LINgroups:
                if "LINgroup_report" not in args.output_format:
                    args.output_format.append("LINgroup_report")
            elif "LINgroup_report" in args.output_format:
                raise ValueError(f"Must provide LINgroup csv via '--LINgroups' in order to output a LINgroup_report.")
        elif args.LINgroups or "LINgroup_report" in args.output_format:
            raise ValueError(f"Must enable LIN taxonomy via '--LIN-taxonomy' in order to output a LINgroup_report.")

        # handle output formats
        print(args.output_format)
        if not args.rank:
            if any(x in ["krona"] for x in args.output_format):
                raise ValueError(f"Rank (--rank) is required for krona output format.")
        if len(args.output_format) > 1:
            if args.output_base == "-":
                raise ValueError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
        elif not args.output_format:
            # change to "human" for 5.0
            args.output_format = ["csv_summary"]

    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        import sys; sys.exit(-1)

    return sourmash.tax.__main__.genome(args)
