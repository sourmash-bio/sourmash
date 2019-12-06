from argparse import FileType
import sourmash
from sourmash._minhash import get_minhash_default_seed
from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('compute')
    subparser.add_argument(
        'filenames', nargs='+', help='file(s) of sequences'
    )

    add_moltype_args(subparser)
    subparser.add_argument(
        '-q', '--quiet', action='store_true', help='suppress non-error output'
    )
    subparser.add_argument(
        '--input-is-protein', action='store_true',
        help='Consume protein sequences - no translation needed.'
    )
    subparser.add_argument(
        '-k', '--ksizes', default='21,31,51',
        help='comma-separated list of k-mer sizes; default=%(default)s'
    )
    subparser.add_argument(
        '-n', '--num-hashes', type=int, default=1500,
        help='number of hashes to use in each sketch; default=%(default)i'
    )
    subparser.add_argument(
        '--check-sequence', action='store_true',
        help='complain if input sequence is invalid'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='recompute signatures even if the file exists'
    )
    subparser.add_argument(
        '-o', '--output', type=FileType('wt'),
        help='output computed signatures to this file'
    )
    subparser.add_argument(
        '--singleton', action='store_true',
        help='compute a signature for each sequence record individually'
    )
    subparser.add_argument(
        '--merge', '--name', type=str, default='', metavar="FILE",
        help='merge all input files into one signature file with the '
        'specified name'
    )
    subparser.add_argument(
        '--name-from-first', action='store_true',
        help='name the signature generated from each file after the first '
        'record in the file'
    )
    subparser.add_argument(
        '--input-is-10x', action='store_true',
        help='input is 10x single cell output folder'
    )
    subparser.add_argument(
        '--count-valid-reads', default=0, type=int,
        help='For 10x input only (i.e input-is-10x flag is True), A barcode '
        'is only considered a valid barcode read and its signature is written '
        'if number of umis are greater than count-valid-reads. It is used to '
        'weed out cell barcodes with few umis that might have been due to '
        'false rna enzyme reactions'
    )
    subparser.add_argument(
        '--write-barcode-meta-csv', type=str,
        help='For 10x input only (i.e input-is-10x flag is True), for each of '
        'the unique barcodes, Write to a given path, number of reads and '
        'number of umis per barcode.'
    )
    subparser.add_argument(
        '-p', '--processes', default=2, type=int,
        help='For 10x input only (i.e input-is-10x flag is True, Number of '
        'processes to use for reading 10x bam file'
    )
    subparser.add_argument(
        '--save-fastas', default="", type=str,
        help='For 10x input only (i.e input-is-10x flag is True), save merged '
        'fastas for all the unique barcodes to {CELL_BARCODE}.fasta in the '
        'absolute path given by this flag, By default, fastas are not saved'
    )
    subparser.add_argument(
        '--line-count', type=int, default=1500,
        help='For 10x input only (i.e input-is-10x flag is True), line count '
        'for each bam shard',
    )
    subparser.add_argument(
        '--track-abundance', action='store_true',
        help='track k-mer abundances in the generated signature'
    )
    subparser.add_argument(
        '--scaled', type=float, default=0,
        help='choose number of hashes as 1 in FRACTION of input k-mers'
    )
    subparser.add_argument(
        '--seed', type=int, default=get_minhash_default_seed(),
        help='seed used by MurmurHash; default=%(default)i'
    )
    subparser.add_argument(
        '--randomize', action='store_true',
        help='shuffle the list of input filenames randomly'
    )
    subparser.add_argument(
        '--license', default='CC0', type=str,
        help='signature license. Currently only CC0 is supported.'
    )
    subparser.add_argument(
        '--rename-10x-barcodes', metavar='FILE',
        help='Tab-separated file mapping 10x barcode name to new name, e.g. '
        'with channel or cell annotation label'
    )
    subparser.add_argument(
        '--barcodes-file', metavar='FILE',
        help='Barcodes file if the input is unfiltered 10x bam file'
    )
