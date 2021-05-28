"""summarize metagenome gather results at rank"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('summarize')
    subparser.add_argument('gather_results')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv',  metavar='FILE',
        help='database lineages csv'
    )
    subparser.add_argument(
        '-r', '--rank', choices=['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'],
        help='Summarize genome taxonomy at this rank and above'
    )
    subparser.add_argument(
        '--split-identifiers', action='store_true',
        help='split names in signatures on whitespace'
    )
    subparser.add_argument(
        '--keep-identifier-versions', action='store_true',
        help='do not remove accession versions'
    )
    subparser.add_argument(
        '--fail-on-missing-taxonomy', action='store_true',
        help='fail quickly if taxonomy is not available for an identifier',
    )

def main(args):
    import sourmash
    return sourmash.tax.__main__.summarize(args)
