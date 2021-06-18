"""classify genomes from gather results"""

import argparse
import sourmash
from sourmash.logging import notify, print_results, error
from sourmash.cli.utils import add_threshold_arg

def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument(
        '-g', '--gather-csv', nargs='*', default = [],
        help='csvs from sourmash gather'
    )
    subparser.add_argument(
        '--from-file',  metavar='FILE', default=None,
        help='input many gather results as a text file, with one gather csv per line'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv',  metavar='FILE',
        nargs='+', required=True,
        help='database lineages csv'
    )
    subparser.add_argument(
        '-o', '--output-base', default='-',
        help='base filepath for output file(s) (default stdout)'
    )
    subparser.add_argument(
        '-r', '--rank', choices=['strain','species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'],
        help='Summarize genome taxonomy at this rank and above. Note that the taxonomy csv must contain lineage information at this rank.'
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
        '--output-format', default=['summary'], nargs='+', choices=["summary", "krona"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past survivable errors in loading taxonomy database or gather results',
    )
    add_threshold_arg(subparser, 0.1)


def main(args):
    import sourmash
    if not args.gather_csv and not args.from_file:
        raise ValueError(f"No gather csvs found! Please input via `-g` or `--from-file`.")
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    if not args.rank:
        if any(x in ["krona", "lineage_summary"] for x in args.output_format):
            raise ValueError(f"Rank (--rank) is required for krona output format.")
    return sourmash.tax.__main__.classify(args)
