"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import os
import os.path
import sys
import random
import screed
import time
import collections

from . import sourmash_args
from .signature import SourmashSignature
from .logging import notify, error, set_quiet
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters)


def parse_params_str(params_str):
    d = {}
    pp = params_str.split(',')
    for p in pp:
        if p == 'abund':
            d['track_abundance'] = True
        elif p == 'noabund':
            d['track_abundance'] = False
        elif p.startswith('k='):
            d['ksize'] = int(p[2:])
        elif p.startswith('num='):
            d['num'] = int(p[4:])
        elif p.startswith('scaled='):
            d['scaled'] = int(p[7:])
        elif p.startswith('seed='):
            d['seed'] = int(p[5:])
        else:
            raise ValueError(f"unknown component '{p}' in params string")

    return d


class _signatures_for_sketch_factory(object):
    "Build sigs on demand, based on args input to 'sketch'."
    def __init__(self, params_str_list):
        self.params_list = []
        for params_str in params_str_list:
            d = parse_params_str(params_str)
            self.params_list.append(d)

    def __call__(self):
        x = []
        for d in self.params_list:
            params = ComputeParameters([d.get('ksize', [31])],
                                        d.get('seed', 42),
                                       0, 0, 0, 1,
                                       d.get('num', 0),
                                       d.get('track_abundance', 0),
                                       d.get('scaled', '1000'))
            sig = SourmashSignature.from_params(params)
            x.append(sig)
        return x


def dna(args):
    """Compute the signature for one or more files.

    CTB: make usable via Python?
    """
    set_quiet(args.quiet)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    notify('computing signatures for files: {}', ", ".join(args.filenames))

    if args.merge and not args.output:
        error("ERROR: must specify -o with --merge")
        sys.exit(-1)

    if args.output and args.outdir:
        error("ERROR: --outdir doesn't make sense with -o/--output")
        sys.exit(-1)

    # TBD/FIXME
    args.input_is_10x = False # CTB
    args.check_sequence = False
    args.input_is_protein = False

    # provide good defaults for dna
    if not args.param_string:
        args.param_string = ['k=31,scaled=1000,noabund']

    # get list of k-mer sizes for which to compute sketches
    num_sigs = len(args.param_string)

    notify('Computing a total of {} signature(s).', num_sigs)

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    signatures_factory = _signatures_for_sketch_factory(args.param_string)

    if args.merge:               # single name specified - combine all
        _compute_merged(args, signatures_factory)
    else:                        # compute individual signatures
        _compute_individual(args, signatures_factory)
