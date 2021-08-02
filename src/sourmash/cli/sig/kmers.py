"""show k-mers/sequences matching the signature hashes"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('kmers')
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
    add_ksize_arg(subparser, 31)
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
