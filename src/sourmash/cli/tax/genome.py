"""classify genomes from gather results"""

usage="""

    sourmash tax genome --gather-csv <gather_csv> [ ... ] --taxonomy-csv <taxonomy-csv> [ ... ]

The 'tax genome' command reads in genome gather result CSVs and reports likely
classification for each query genome.

By default, classification uses a containment threshold of 0.1, meaning at least
10 percent of the query was covered by matches with the reported taxonomic rank and lineage.
You can specify an alternate classification threshold or force classification by
taxonomic rank instead, e.g. at species or genus-level.

The default output format consists of five columns,
 `query_name,status,rank,fraction,lineage`, where `fraction` is the fraction
 of the query matched to the reported rank and lineage. The `status` column
 provides additional information on the classification, and can be:
  - `match` - this query was classified
  - `nomatch`- this query could not be classified
  - `below_threshold` - this query was classified at the specified rank,
     but the query fraction matched was below the containment threshold

Optionally, you can report classifications in `krona` format, but note
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
        '-g', '--gather-csv', nargs='*', default = [],
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
        nargs='+', required=True,
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
        '--output-format', default=['csv_summary'], nargs='+', choices=["csv_summary", "krona"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past survivable errors in loading taxonomy database or gather results',
    )
    add_tax_threshold_arg(subparser, 0.1)


def main(args):
    import sourmash
    if not args.gather_csv and not args.from_file:
        raise ValueError(f"No gather CSVs found! Please input via `-g` or `--from-file`.")
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    if not args.rank:
        if any(x in ["krona"] for x in args.output_format):
            raise ValueError(f"Rank (--rank) is required for krona output format.")
    return sourmash.tax.__main__.genome(args)
