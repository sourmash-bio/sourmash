"""create signatures from a CSV file"""

usage="""

    sourmash sketch fromfile <csv file> --output-signatures <location> -p <...>

The 'sketch fromfile' command takes in a CSV file with list of names
and filenames to be used for building signatures. It is intended for
batch use, when building large collections of signatures.

One or more parameter strings must be specified with '-p'.

One or more existing collections of signatures can be provided via
'--already-done' and already-existing signatures (based on name and
sketch type) will not be recalculated or output.

If a location is provided via '--output-signatures', signatures will be saved
to that location.

Please see the 'sketch' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html
"""

import sourmash
from sourmash.logging import notify, print_results, error

from sourmash import command_sketch


def subparser(subparsers):
    subparser = subparsers.add_parser('fromfile',
                                      usage=usage)
    subparser.add_argument(
        'csvs', nargs='+',
        help="input CSVs providing 'name', 'genome_filename', and 'protein_filename'"
    )
    subparser.add_argument(
        '-p', '--param-string', default=[],
        help='signature parameters to use.', action='append',
    )
    subparser.add_argument(
        '--already-done', nargs='+', default=[],
        help='one or more collections of existing signatures to avoid recalculating'
    )
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
        '-o', '--output-signatures',
        help='output computed signatures to this file',
    )
    file_args.add_argument(
        '--force-output-already-exists', action='store_true',
        help='overwrite/append to --output-signatures location'
    )
    file_args.add_argument(
        '--ignore-missing', action='store_true',
        help='proceed with building possible signatures, even if some input files are missing'
    )
    file_args.add_argument(
        '--output-csv-info',
        help='output information about what signatures need to be generated'
    )
    file_args.add_argument(
        '--output-manifest-matching',
        help='output a manifest file of already-existing signatures'
    )


def main(args):
    import sourmash.command_sketch
    return sourmash.command_sketch.fromfile(args)
