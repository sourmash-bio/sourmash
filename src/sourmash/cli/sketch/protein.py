"""create protein signatures"""

usage="""

    sourmash sketch protein data/*.fna.gz

The 'sketch protein' command reads in protein sequences and outputs protein
sketches.

By default, 'sketch protein' uses the parameter string
'k=10,scaled=200,noabund'.

This corresponds to an amino-acid k-mer size of 10, a scaled factor
of 200, and no abundance tracking of k-mers. You can specify one or
more parameter strings of your own with -p, e.g. 'sourmash sketch
protein -p k=11,noabund -p k=12,scaled=100,abund'. Note that a single `-p` parameter string can contain multiple ksize values, but only a single scaled value or abundance value e.g. -p k=11,k=12,scaled=100,abund.

'sourmash sketch' takes input sequences in FASTA and FASTQ,
uncompressed or gz/bz2 compressed.

Please see the 'sketch' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html
"""

import sourmash
from sourmash.logging import notify, print_results, error

from sourmash import command_sketch
assert command_sketch.DEFAULTS['protein'] == 'k=10,scaled=200,noabund'


def subparser(subparsers):
    subparser = subparsers.add_parser('protein', aliases=['aa', 'prot'],
                                      usage=usage)
    subparser.add_argument(
        '--license', default='CC0', type=str,
        help='signature license. Currently only CC0 is supported.'
    )
    subparser.add_argument(
        '-p', '--param-string', default=[],
        help='signature parameters to use.', action='append',
    )
    
    subparser.add_argument(
        'filenames', nargs='*', help='file(s) of sequences'
    )
    file_args = subparser.add_argument_group('File handling options')
    file_args.add_argument(
        '-f', '--force', action='store_true',
        help='recompute signatures even if the file exists'
    )
    file_args.add_argument(
        '-o', '--output',
        help='output computed signatures to this file'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of sequence files to load'
    )
    file_args.add_argument(
        '--merge', '--name', type=str, default='', metavar="FILE",
        help='merge all input files into one signature file with the '
        'specified name'
    )
    file_args.add_argument(
        '--output-dir', '--outdir',
        help='output computed signatures to this directory',
    )
    file_args.add_argument(
        '--singleton', action='store_true',
        help='compute a signature for each sequence record individually'
    )
    file_args.add_argument(
        '--name-from-first', action='store_true',
        help='name the signature generated from each file after the first '
        'record in the file'
    )
    file_args.add_argument(
        '--randomize', action='store_true',
        help='shuffle the list of input filenames randomly'
    )
    file_args.add_argument(
        '--dayhoff', action='store_true',
        help='compute sketches using the dayhoff alphabet instead'
    )
    file_args.add_argument(
        '--hp', action='store_true',
        help='compute sketches using the dayhoff alphabet instead'
    )


def main(args):
    import sourmash.command_sketch
    return sourmash.command_sketch.protein(args)
