"""summarize metagenome gather results at rank"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('summarize')
    subparser.add_argument('gather_results', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output-base', default='-',
        help='base filepath for output file(s) (default stdout)'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv',  metavar='FILE',
        help='database lineages csv'
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
        '-r', '--rank', choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'], # strain?
        help='For non-default output formats: Summarize genome taxonomy at this rank and above'
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in taxonomy database loading',
    )

def main(args):
    import sourmash
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    return sourmash.tax.__main__.summarize(args)
