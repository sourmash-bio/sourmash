"""create signatures from a CSV file"""

usage="""

    sourmash sketch fromfile @@CTB

@CTB

'sourmash sketch' takes input sequences in FASTA and FASTQ,
uncompressed or gz/bz2 compressed.

Please see the 'sketch' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html
"""

import sourmash
from sourmash.logging import notify, print_results, error

from sourmash import command_sketch


def subparser(subparsers):
    subparser = subparsers.add_parser('fromfile',
                                      usage=usage)
    # @CTB
    subparser.add_argument('csv')
    subparser.add_argument(
        '-p', '--param-string', default=[],
        help='signature parameters to use.', action='append',
    )
    # @CTB
    subparser.add_argument('--already-done', nargs='+', default=[])
    # @CTB
    subparser.add_argument('--force-output-already-exists', action='store_true')
    subparser.add_argument(
        '--license', default='CC0', type=str,
        help='signature license. Currently only CC0 is supported.'
    )
    subparser.add_argument(
        '--check-sequence', action='store_true',
        help='complain if input sequence is invalid (NOTE: only checks DNA)'
    )
    file_args = subparser.add_argument_group('File handling options')
    file_args.add_argument(
        '-o', '--output',
        help='output computed signatures to this file',
        required=True
    )


def main(args):
    import sourmash.command_sketch
    return sourmash.command_sketch.fromfile(args)
