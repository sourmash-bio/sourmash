from glob import glob
import os


def add_moltype_args(parser):
    parser.add_argument(
        '--protein', dest='protein', action='store_true',
        help='choose a protein signature; by default, a nucleotide signature is used'
    )
    parser.add_argument(
        '--no-protein', dest='protein', action='store_false',
        help='do not choose a protein signature')
    parser.set_defaults(protein=False)

    parser.add_argument(
        '--dayhoff', dest='dayhoff', action='store_true',
        help='build Dayhoff-encoded amino acid signatures'
    )
    parser.add_argument(
        '--no-dayhoff', dest='dayhoff', action='store_false',
        help='do not build Dayhoff-encoded amino acid signatures')
    parser.set_defaults(dayhoff=False)

    parser.add_argument(
        '--hp', '--hydrophobic-polar', dest='hp', action='store_true',
        help='build hydrophobic-polar-encoded amino acid signatures'
    )
    parser.add_argument(
        '--no-hp', '--no-hydrophobic-polar', dest='hp', action='store_false',
        help='do not build hydrophobic-polar-encoded amino acid signatures')
    parser.set_defaults(hp=False)

    parser.add_argument(
        '--dna', '--rna', dest='dna', default=None, action='store_true',
        help='choose a nucleotide signature (default: True)')
    parser.add_argument(
        '--no-dna', '--no-rna', dest='dna', action='store_false',
        help='do not choose a nucleotide signature')
    parser.set_defaults(dna=None)


def add_construct_moltype_args(parser):
    add_moltype_args(parser)
    parser.set_defaults(dna=True)


def add_ksize_arg(parser, default=31):
    parser.add_argument(
        '-k', '--ksize', metavar='K', default=None, type=int,
        help='k-mer size; default={d}'.format(d=default)
    )


def opfilter(path):
    return not path.startswith('__') and path not in ['utils']


def command_list(dirpath):
    paths = glob(os.path.join(dirpath, '*.py'))
    filenames = [os.path.basename(path) for path in paths]
    basenames = [os.path.splitext(path)[0] for path in filenames if not path.startswith('__')]
    basenames = filter(opfilter, basenames)
    return sorted(basenames)
