"""search taxonomies and output picklists.."""

usage="""

    sourmash tax grep <term> --taxonomy-csv <taxonomy_file> [ ... ]

`sourmash tax grep` searches taxonomies for matching strings,
optionally restricting the string search to a specific taxonomic rank.
It creates new files containing matching taxonomic entries; these new
files can serve as taxonomies and can also be used as picklists.

Please see the 'tax grep' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-grep-subset-taxonomies-and-create-picklists-based-on-taxonomy-string-matches
"""

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('grep', usage=usage)
    subparser.add_argument('pattern')
    subparser.add_argument('-r', '--rank',
                           help="search only this rank",
                           choices=['superkingdom',
                                    'phylum',
                                    'class',
                                    'order',
                                    'family',
                                    'genus',
                                    'species'])
    subparser.add_argument(
        '-v', '--invert-match',
        help="select non-matching lineages",
        action="store_true"
    )
    subparser.add_argument(
        '-i', '--ignore-case',
        help="ignore case distinctions (search lower and upper case both)",
        action="store_true"
    )
    subparser.add_argument(
        '--silent', '--no-picklist-output',
        help="do not output picklist",
        action='store_true',
    )
    subparser.add_argument(
        '-c', '--count',
        help="only output a count of discovered lineages; implies --silent",
        action='store_true'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv', '--taxonomy', metavar='FILE',
        nargs="+", required=True, action="extend",
        help='database lineages'
    )
    subparser.add_argument(
        '-o', '--output', default='-',
        help='output file (defaults to stdout)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past errors in file and taxonomy loading',
    )

def main(args):
    import sourmash
    return sourmash.tax.__main__.grep(args)
