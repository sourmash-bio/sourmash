from glob import glob
import os
import argparse
from sourmash.logging import notify
from sourmash.sourmash_args import check_scaled_bounds, check_num_bounds


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
        '--dna', '--rna', '--nucleotide', dest='dna', default=None, action='store_true',
        help='choose a nucleotide signature (default: True)')
    parser.add_argument(
        '--no-dna', '--no-rna', '--no-nucleotide', dest='dna', action='store_false',
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

#https://stackoverflow.com/questions/55324449/how-to-specify-a-minimum-or-maximum-float-value-with-argparse#55410582
def range_limited_float_type(arg):
    """ Type function for argparse - a float within some predefined bounds """
    min_val = 0
    max_val = 1
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("\n\tERROR: Must be a floating point number.")
    if f < min_val or f > max_val:
        raise argparse.ArgumentTypeError(f"\n\tERROR: Argument must be >{str(min_val)} and <{str(max_val)}.")
    return f


def add_tax_threshold_arg(parser, default=0.1):
    parser.add_argument(
        '--containment-threshold', default=default, type=range_limited_float_type,
        help=f'minimum containment threshold for classification; default={default}'
    )


def add_picklist_args(parser):
    parser.add_argument(
        '--picklist', default=None,
        help="select signatures based on a picklist, i.e. 'file.csv:colname:coltype'"
    )
    parser.add_argument(
        '--picklist-require-all', default=False, action='store_true',
        help="require that all picklist values be found or else fail"
    )


def opfilter(path):
    return not path.startswith('__') and path not in ['utils']


def command_list(dirpath):
    paths = glob(os.path.join(dirpath, '*.py'))
    filenames = [os.path.basename(path) for path in paths]
    basenames = [os.path.splitext(path)[0] for path in filenames if not path.startswith('__')]
    basenames = filter(opfilter, basenames)
    return sorted(basenames)


def add_scaled_arg(parser, default=None):
    parser.add_argument(
        '--scaled', metavar='FLOAT', type=check_scaled_bounds,
        help='scaled value should be between 100 and 1e6'
    )


def add_num_arg(parser, default=0):
    parser.add_argument(
        '-n', '--num-hashes', '--num', metavar='N', type=check_num_bounds, default=default,
        help='num value should be between 50 and 50000'
    )