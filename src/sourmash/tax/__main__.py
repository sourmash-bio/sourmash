"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import json
import os
from collections import defaultdict

import sourmash
import copy
from sourmash.sourmash_args import FileOutput
from sourmash.lca.lca_utils import pop_to_rank
from sourmash.lca.command_index import load_taxonomy_assignments

from ..sourmash_args import FileOutputCSV

from sourmash.logging import set_quiet, error, notify, set_quiet, print_results, debug
from sourmash import sourmash_args
from sourmash.minhash import _get_max_hash_for_scaled

from . import tax_utils

usage='''
sourmash taxonomy <command> [<args>] - manipulate/work with taxonomy information.
or
sourmash tax <command> [<args>]


** Commands can be:

summarize <gather_results> [<gather_results> ... ]        - summarize taxonomic information for metagenome gather results
classify <gather_results> [<gather_results> ... ]   - taxonomic classification of genomes from gather results

** Use '-h' to get subcommand-specific help, e.g.

sourmash taxonomy summarize -h
'''

##### taxonomy command line functions

def summarize(args):
    """
    summarize taxonomic information for metagenome gather results
    """
    set_quiet(args.quiet)

    # load gather results and taxonomy assignments
    gather_results = tax_utils.load_gather_results(args.gather_results)
    tax_assign, _ = load_taxonomy_assignments(args.taxonomy_csv, use_headers=True, force=False, split_identifiers=args.split_identifiers, keep_identifier_versions = args.keep_identifier_versions)

    # check for match identites not found in lineage spreadsheets
    n_missed, ident_missed = tax_utils.find_missing_identities(gather_results, tax_assign)
    if n_missed:
        print(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
    assert n_missed == 0

    # write output csv
    csv_fp = None
    with FileOutputCSV(args.output) as csv_fp:
        w = csv.writer(csv_fp)
        # actually summarize at rank
        for rank in sourmash.lca.taxlist(include_strain=False): # do we need to do this at all ranks?
            g_at_rank = tax_utils.summarize_gather_at(rank, tax_assign, gather_results)
            for k, v in g_at_rank:
                w.writerow([rank, f'{v:.3f}', sourmash.lca.display_lineage(k)])
    if csv_fp:
        csv_fp.close()


def classify(args):
    """
    taxonomic classification of genomes from gather results
    """
    set_quiet(args.quiet)
    print("entered classify command")

def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
