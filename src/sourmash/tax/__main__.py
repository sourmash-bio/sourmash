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

# some utils
def make_outfile(base, ext):
    if base == "-":
        return base
    return base + ext

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
        notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
    assert n_missed == 0

    # actually summarize at rank
    summarized_gather = {}
    for rank in sourmash.lca.taxlist(include_strain=False):
        summarized_gather[rank] = tax_utils.summarize_gather_at(rank, tax_assign, gather_results)

    # write summarozed output csv
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".summarized.csv")
        header= ["rank", "fraction", "lineage"]
        csv_fp = None
        with FileOutputCSV(summary_outfile) as csv_fp:
            w = csv.writer(csv_fp)
            w.writerow(header)
            for rank, rank_results in summarized_gather.items():
                for sorted_result in rank_results:
                    lin,val = sorted_result
                    w.writerow([rank, f'{val:.3f}', sourmash.lca.display_lineage(lin)])

    # write summarized --> krona output csv
    if "krona" in args.output_format:
        krona_resultslist = tax_utils.format_for_krona(args.rank, summarized_gather)

        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_resultslist, out_fp)


# todo -- fix for new output file format
def classify(args):
    """
    taxonomic classification of genomes from gather results
    """
    ## currently reports a single rank. do we want to optionally report at all ranks? (no, bc summarize does that?)
    set_quiet(args.quiet)

    # load taxonomy assignments
    tax_assign, _ = load_taxonomy_assignments(args.taxonomy_csv, use_headers=True, force=False, split_identifiers=args.split_identifiers, keep_identifier_versions = args.keep_identifier_versions)

    # write output csv
    header= ["rank", "fraction", "lineage"]
    csv_fp = None
    with FileOutputCSV(args.output) as csv_fp:
        w = csv.writer(csv_fp)
        w.writerow(header)

        # load gather results and taxonomy assignments
        for g_result in args.gather_results:
            gather_results = tax_utils.load_gather_results(g_result)

            # check for match identites not found in lineage spreadsheets
            n_missed, ident_missed = tax_utils.find_missing_identities(gather_results, tax_assign)
            if n_missed:
                notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
            assert n_missed == 0


            # if --rank is specified, classify to that rank
            # to do, what to do if don't have gather results at desired rank (e.g. strain)?
            if args.rank:
                # todo: check we have gather results at this rank
                #if not tax_utils.check_taxonomy_exists(tax_assign, args.rank):
                #    notify(f"No taxonomic information at rank {args.rank}: cannot classify at this rank")
                (lineage,containment) = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, best_only=True)
                w.writerow([args.rank, f'{containment:.3f}', sourmash.lca.display_lineage(lineage)])
                if containment <= args.containment_threshold:
                    notify(f"WARNING: classifying at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
            else:
                # classify to the match that passes the containment threshold. To do - do we want to report anything if nothing >= containment threshold?
                for rank in tax_utils.ascending_taxlist(include_strain=False):
                    (lineage,containment) = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, best_only=True)
                    if containment >= args.containment_threshold:
                        w.writerow([rank, f'{containment:.3f}', sourmash.lca.display_lineage(lineage)])
                        break
    if csv_fp:
        csv_fp.close()


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
