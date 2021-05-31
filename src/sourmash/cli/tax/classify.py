"""classify genomes"""

import argparse
import sourmash
from sourmash.logging import notify, print_results, error

#https://stackoverflow.com/questions/55324449/how-to-specify-a-minimum-or-maximum-float-value-with-argparse#55410582
# should this go in a different file?
def range_limited_float_type(arg):
    """ Type function for argparse - a float within some predefined bounds """
    min_val = 0
    max_val = 1
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < min_val or f > max_val:
        raise argparse.ArgumentTypeError(f"Argument must be >{str(min_val)} and <{str(max_val)}")
    return f


def subparser(subparsers):
    subparser = subparsers.add_parser('classify')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-t', '--taxonomy-csv',  metavar='FILE',
        help='database lineages csv'
    )
    subparser.add_argument(
        '-g', '--gather-results',  metavar='FILE',
        help='database lineages csv'
    )
    subparser.add_argument(
        '-n', '--query-name', default="",
        help='name of query to be classified'
    )
    subparser.add_argument(
        '--from-csv',  metavar='FILE',
        # to do: if query_name in gather results, can just have textfile of gather_results here
        help='input many gather results as a csv with "name,resultsfile" on each line'
    )
    subparser.add_argument(
        '-o', '--output-base', default='-',
        help='base filepath for output file(s) (default stdout)'
    )
    subparser.add_argument(
        '-r', '--rank', choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'], #strain
        help='Summarize genome taxonomy at this rank and above'
    )
    subparser.add_argument(
        '--containment-threshold', type=range_limited_float_type, default=0.1,
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
    subparser.add_argument(
        '--output-format', default=['summary'], nargs='+', choices=["summary", "krona"],
        help='choose output format(s)',
    )
    subparser.add_argument(
        '-f', '--force', action = 'store_true',
        help='continue past survivable errors in loading taxonomy database or gather results',
    )


def main(args):
    import sourmash
    if len(args.output_format) > 1:
        if args.output_base == "-":
            raise TypeError(f"Writing to stdout is incompatible with multiple output formats {args.output_format}")
    return sourmash.tax.__main__.classify(args)
