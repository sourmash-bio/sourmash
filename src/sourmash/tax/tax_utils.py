"""
Utility functions for taxonomy analysis tools.
"""
import csv
from os.path import exists
from collections import namedtuple, defaultdict, Counter

__all__ = ['get_ident', 'load_gather_results',
           'summarize_gather_at', 'find_missing_identities']

from sourmash.logging import notify, error, debug

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, build_tree, find_lca,
                                    taxlist, count_lca_for_assignments,
                                    zip_lineage, display_lineage,
                                    make_lineage, is_lineage_match,
                                    pop_to_rank)


def get_ident(ident):
    "Hack and slash identifiers."
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident

# load and aggregate all gather results
def load_gather_results(gather_csvs):
    gather_results = []
    for g_csv in gather_csvs:
        with open(g_csv, 'rt') as fp:
            r = csv.DictReader(fp)
            for n, row in enumerate(r):
                gather_results.append(row)
        print(f'loaded {str(n)} gather results from {g_csv}.')
    print(f'loaded {len(gather_results)} gather results in total.')
    return gather_results


# this summarizes at a specific rank.
# want to also have a flexible version that goes up a rank
# if needed for good lca
def summarize_gather_at(rank, tax_assign, gather_results):
    # collect!
    sum_uniq_weighted = defaultdict(float)
    for row in gather_results:
        # move these checks to loading function!
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        lineage = tax_assign[match_ident]
        # actual summarization code
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]

        f_uniq_weighted = row['f_unique_weighted']
        f_uniq_weighted = float(f_uniq_weighted)
        sum_uniq_weighted[lineage] += f_uniq_weighted

    items = list(sum_uniq_weighted.items())
    items.sort(key = lambda x: -x[1])
    return items


def find_missing_identities(gather_results, tax_assign):
    n_missed = 0
    ident_missed= []
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            n_missed += 1
            ident_missed.append(match_ident)

    print(f'of {len(gather_results)}, missed {n_missed} lineage assignments.')
    return n_missed, ident_missed
