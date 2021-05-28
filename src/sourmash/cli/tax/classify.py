"""classify genomes"""

import sourmash
from sourmash.logging import notify, print_results, error

#https://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __contains__(self, item):
        return self.__eq__(item)

    def __iter__(self):
        yield self

    def __repr__(self):
        return f'[{self.start}, {self.end}]'


def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument('gather_results', nargs='+')
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
        '-r', '--rank', choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'], #strain
        help='Summarize genome taxonomy at this rank and above'
    )
    subparser.add_argument(
        '--containment-threshold', type=float, default=0.1, choices=[Range(0.0, 1.0)],
        help='minimum containment for classification'
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
    return sourmash.tax.__main__.classify(args)
