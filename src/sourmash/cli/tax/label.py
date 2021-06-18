"""add taxonomy information to gather results"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('label')
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
        '-t', '--taxonomy-csv', metavar='FILE',
        nargs="+", required=True,
        help='database lineages csv'
    )
    subparser.add_argument(
        '-o', '--output-dir', default= "",
        help='directory for output csv(s)'
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
        raise ValueError(f"No gather csvs found! Please input via `-g` or `--from-file`.")
    return sourmash.tax.__main__.label(args)
