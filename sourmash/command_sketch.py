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
from .signature import SourmashSignature
from .logging import notify, error, set_quiet
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters)


class _signatures_for_sketch_factory(object):
    "Build sigs on demand, based on args input to 'sketch'."
    def __init__(self, args):
        self.args = args

    def __call__(self):
        args = self.args
        params = ComputeParameters(args.ksizes, args.seed, args.protein,
                                   args.dayhoff, args.hp, args.dna,
                                   args.num_hashes,
                                   args.track_abundance, args.scaled)
        sig = SourmashSignature.from_params(params)
        return [sig]


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

    # regular override
    args.track_abundance = False # CTB
    args.ksizes = [ 31 ]
    args.num_hashes = 0
    args.scaled = 1000

    # TBD
    args.input_is_10x = False # CTB
    args.check_sequence = False

    # rare override
    args.seed = 42

    # fixed parameters:
    args.dna = True

    # irrelevant-to-DNA parameters:
    args.input_is_protein = False
    args.protein = False
    args.hp = False
    args.dayhoff = False

    if not args.param_string:
        args.param_string = ['k=31,scaled=1000,noabund']
    print('XXX', args.param_string)

    signatures_factory = _signatures_for_sketch_factory(args)

    if args.merge:               # single name specified - combine all
        _compute_merged(args, signatures_factory)
    else:                        # compute individual signatures
        _compute_individual(args, signatures_factory)
