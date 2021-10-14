"""create DNA signatures"""

usage="""

    sourmash sketch dna data/*.fna.gz

The 'sketch dna' command reads in DNA sequences and outputs DNA
sketches.

By default, 'sketch dna' uses the parameter string 'k=31,scaled=1000,noabund'.

This creates sketches with a k-mer size of 31, a scaled factor of
1000, and no abundance tracking of k-mers.  You can specify one or
more parameter strings of your own with -p, e.g.  'sourmash sketch dna
-p k=31,noabund -p k=21,scaled=100,abund'. Note that a single `-p` parameter string can contain multiple ksize values, but only a single scaled value or abundance value, e.g. -p k=21,k=31,abund

'sourmash sketch' takes input sequences in FASTA and FASTQ,
uncompressed or gz/bz2 compressed.

Please see the 'sketch' documentation for more details:
  https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html
"""

import sourmash
from sourmash.logging import notify, print_results, error

from sourmash import command_sketch
assert command_sketch.DEFAULTS['dna'] == 'k=31,scaled=1000,noabund'


def subparser(subparsers):
    subparser = subparsers.add_parser('dna',
                                      aliases=['rna', 'nucleotide', 'nt'],
                                      usage=usage)
    subparser.add_argument(
        '--license', default='CC0', type=str,
        help='signature license. Currently only CC0 is supported.'
    )
    subparser.add_argument(
        '--check-sequence', action='store_true',
        help='complain if input sequence is invalid (NOTE: only checks DNA)'
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
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of sequence files to load'
    )
    file_args.add_argument(
        '-o', '--output',
        help='output computed signatures to this file'
    )
    file_args.add_argument(
        '--merge', '--name', type=str, default='', metavar="FILE",
        help='merge all input files into one signature file with the '
        'specified name'
    )
    file_args.add_argument(
        '--outdir', help='output computed signatures to this directory'
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


def main(args):
    import sourmash.command_sketch
    return sourmash.command_sketch.dna(args)
