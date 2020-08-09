"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import os
import os.path
import sys
import random
import screed
import time

from . import sourmash_args
from .signature import SourmashSignature, save_signatures
from .logging import notify, error, set_quiet
from .utils import RustObject
from ._lowlevel import ffi, lib
from .command_compute import (_compute_individual, _compute_merged,
                              make_minhashes, add_seq, build_siglist,
                              save_siglist)

DEFAULT_COMPUTE_K = '21,31,51'
DEFAULT_LINE_COUNT = 1500


def dna(args):
    """Compute the signature for one or more files.

    CTB: make usable via Python?
    """
    set_quiet(args.quiet)

#    if args.license != 'CC0':
#        error('error: sourmash only supports CC0-licensed signatures. sorry!')
#        sys.exit(-1)

    notify('computing signatures for files: {}', ", ".join(args.filenames))

    # get list of k-mer sizes for which to compute sketches
    #ksizes = args.ksizes
    ksizes = [31]

    notify('Computing signature for ksizes: {}', str(ksizes))
    num_sigs = len(ksizes)

    notify('Computing a total of {} signature(s).', num_sigs)

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    if args.merge and not args.output:
        error("ERROR: must specify -o with --merge")
        sys.exit(-1)

    if args.output and args.outdir:
        error("ERROR: --outdir doesn't make sense with -o/--output")
        sys.exit(-1)

#    if args.track_abundance:
#        notify('Tracking abundance of input k-mers.')

    args.input_is_10x = False # CTB
    args.track_abundance = False # CTB
    args.ksizes = [ 31 ]
    args.seed = 42
    args.protein = False
    args.hp = False
    args.dayhoff = False
    args.dna = True
    args.num_hashes = 0
    args.scaled = 1000
    args.input_is_protein = False
    args.check_sequence = False

    if args.merge:               # single name specified - combine all
        _compute_merged(args)
    else:                        # compute individual signatures
        _compute_individual(args)
