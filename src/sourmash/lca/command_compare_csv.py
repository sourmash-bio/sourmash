#! /usr/bin/env python
"""
Compare two taxonomy spreadsheets.
"""
import sys
from collections import defaultdict

from ..logging import notify, error, print_results, set_quiet
from . import lca_utils
from .lca_utils import zip_lineage
from .command_index import load_taxonomy_assignments


def compare_csv(args):
    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    # first, load classify-style spreadsheet
    notify(f'loading classify output from: {args.csv1}')
    assignments0, num_rows0 = load_taxonomy_assignments(args.csv1,
                                                        start_column=3,
                                                        force=args.force)

    notify(f'loaded {len(set(assignments0.values()))} distinct lineages, {num_rows0} rows')
    notify('----')

    # next, load custom taxonomy spreadsheet
    delimiter = ','
    if args.tabs:
        delimiter = '\t'

    notify(f'loading custom spreadsheet from: {args.csv2}')
    assignments, num_rows = load_taxonomy_assignments(args.csv2,
                                               delimiter=delimiter,
                                               start_column=args.start_column,
                                               use_headers=not args.no_headers,
                                               force=args.force)
    notify(f'loaded {len(set(assignments.values()))} distinct lineages, {num_rows} rows')

    # now, compute basic differences:
    missing_1 = set(assignments0.keys()) - set(assignments.keys())
    missing_2 = set(assignments.keys()) - set(assignments0.keys())
    if missing_2:
        notify(f'missing {len(missing_2)} assignments in classify spreadsheet.')
    if missing_1:
        notify(f'missing {len(missing_1)} assignments in custom spreadsheet.')
    if missing_1 or missing_2:
        notify('(these will not be evaluated any further)')
    else:
        notify('note: all IDs are in both spreadsheets!')

    # next, look at differences in lineages
    common = set(assignments0.keys())
    common.intersection_update(assignments.keys())

    n_total = 0
    n_different = 0
    n_compat = 0
    n_incompat = 0
    incompat_rank = defaultdict(int)
    for k in common:
        n_total += 1
        v0 = assignments0[k]
        v1 = assignments[k]
        if v0 != v1:
            n_different += 1
            tree = lca_utils.build_tree([v0])
            lca_utils.build_tree([v1], tree)

            lca, reason = lca_utils.find_lca(tree)
            if reason == 0:               # compatible lineages
                n_compat += 1
                print_results("{},compatible,{}", k, ";".join(zip_lineage(lca)))
            else:
                n_incompat += 1
                print_results("{},incompatible,{}", k, ";".join(zip_lineage(lca)))
                rank = next(iter(lca_utils.taxlist()))
                if lca:
                    rank = lca[-1].rank
                incompat_rank[rank] += 1

    notify(f"{n_total} total assignments, {n_different} differ between spreadsheets.")
    notify(f"{n_compat} are compatible (one lineage is ancestor of another.")
    notify(f"{n_incompat} are incompatible (there is a disagreement in the trees).")

    if n_incompat:
        for rank in lca_utils.taxlist():
            notify(f'{incompat_rank[rank]} incompatible at rank {rank}')
        

if __name__ == '__main__':
    sys.exit(compare_csv(sys.argv[1:]))
