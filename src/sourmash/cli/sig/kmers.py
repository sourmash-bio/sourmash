"""show k-mers/sequences matching the signature hashes"""

usage="""

### `sourmash signature kmers` - extract k-mers and/or sequences that match to signatures

Given one or more compatible sketches and some sequence files, extract
the k-mers and/or sequences corresponding to the hash values in the
sketch. Because the sourmash hash function is one-way, this requires
FASTA or FASTQ sequence files in addition to the sketch.

For example,

sourmash sig kmers --signatures sig1.sig --sequences seqfile.fasta \
    --save-sequences matches.fasta --save-kmers kmer-matches.csv

will search `seqfile.fasta` for matching sequences and k-mers,
and produce two files. The file `matches.fasta` will contain FASTA
sequences that match the hashes in the input signature, while the
file `kmer-matches.csv` provides the matching k-mers and hash values,
together with their originating filename and sequence name.

If the sketch is a protein sketch (protein, dayhoff, or hp), then
the input sequences are assumed to be protein. To search DNA sequences
for translated protein hashes, provide the `--translate` flag to `sig kmers`.

`--save-sequences` and `--save-kmers` are both optional.  If neither are
given, basic statistics on k-mer matching are given.

Please note that `--save-kmers` can be very slow on large files!

The input sketches are the source of the input hashes.  So, for example,
If `--scaled=1` sketches are provided, `sig kmers` can be used to
yield all the k-mers and their matching hashes.  Likewise, if the
sketch is built from the intersection of two other sketches, only
the k-mers and hash values present in both sketches will be used.

Likewise, the input sequences are used for matching; they do not need
to be the same sequences that were used to create the sketches.
Input sequences can be in FASTA or FASTQ format, and either flat text
or compressed with gzip or bzip2; formats are auto-detected.

By default, `sig kmers` ignores bad k-mers (e.g. non-ACGT characters
in DNA). If `--check-sequence` is provided, `sig kmers` will error
exit on the first bad k-mer.  If `--check-sequence --force` is provided,
`sig kmers` will provide error messages (and skip bad sequences), but
will continue processing input sequences.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('kmers', description=__doc__, usage=usage)
    subparser.add_argument('--signatures', nargs='*', default=[])
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)

    subparser.add_argument('--sequences', nargs='+', required=True,
                           help="FASTA/FASTQ/bz2/gz files with sequences")

    subparser.add_argument('--save-kmers',
                           help="save k-mers and hash values to a CSV file")
    subparser.add_argument('--save-sequences',
                           help="save sequences with matching hashes to a FASTA file")
    subparser.add_argument('--translate', action="store_true",
                           help="translate DNA k-mers into amino acids (for protein, dayhoff, and hp sketches)")
    subparser.add_argument(
        '--check-sequence', action='store_true',
        help='complain if input sequence is invalid (NOTE: only checks DNA)'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.kmers(args)
