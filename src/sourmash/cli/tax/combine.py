"""aggregate summarize metagenome gather results at rank"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('combine')
    subparser.add_argument('summarized_gather_results', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--from-file',  metavar='FILE',
        help='input many gather results as a text file, with one gather csv per line'
    )
    subparser.add_argument(
        '-o', '--output-base', default='-',
        help='basename for output file (default stdout)'
    )
    subparser.add_argument(
        '--output-format', default=['csv'], nargs='+', choices=["csv", "tsv"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-r', '--rank', choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'],
        default='species',
        help='Output combined info for lineages at this rank'
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in file loading',
    )

def main(args):
    import sourmash
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    return sourmash.tax.__main__.combine(args)
