#! /usr/bin/env python
"""
@@
"""
from __future__ import print_function, division
import sys
import argparse
import csv
from collections import Counter, defaultdict, namedtuple

import sourmash_lib
from sourmash_lib import sourmash_args
from sourmash_lib.logging import notify, error, print_results
from sourmash_lib.lca import lca_utils
from sourmash_lib.lca.lca_utils import debug, set_debug

LCAGatherResult = namedtuple('LCAGatherResult',
                             'intersect_bp, f_orig_query, f_unique_to_query, f_unique_weighted, average_abund, lineage')


# pretty-printing code. redundant with ../search.py; fix when refactoring.
def format_bp(bp):
    bp = float(bp)
    if bp < 500:
        return '{:.0f} bp'.format(bp)
    elif bp <= 500e3:
        return '{:.1f} kbp'.format(round(bp / 1e3, 1))
    elif bp < 500e6:
        return '{:.1f} Mbp'.format(round(bp / 1e6, 1))
    elif bp < 500e9:
        return '{:.1f} Gbp'.format(round(bp / 1e9, 1))
    return '???'


def format_lineage(lineage_tup):
    """
    Pretty print lineage.
    """
    # list of ranks present
    present = [ l.rank for l in lineage_tup if l.name ]
    d = dict(lineage_tup) # rank: value

    if 'genus' in present:
        if 'species' in present:
            name = '{} {}'.format(d['genus'], d['species'])
        else:
            name = '{} sp.'.format(d['genus'])
    elif len(present) < 3:
        lineage_str = lca_utils.zip_lineage(lineage_tup, truncate_empty=True)
        lineage_str = "; ".join(lineage_str)
        name = lineage_str + ' - (no further assignment)'
    elif len(present) > 1 and 'superkingdom' in present:
        lowest_rank = present[-1]
        name = '{}; {} {}'.format(d['superkingdom'], lowest_rank,
                                   d[lowest_rank])
    else:
        lineage_str = lca_utils.zip_lineage(lineage_tup, truncate_empty=True)
        lineage_str = "; ".join(lineage_str)
        name = lineage_str

    return name


def gather_signature(query_sig, dblist, threshold):
    """
    Decompose 'query_sig' using the given list of databases.
    """
    notify('loaded query: {}... (k={})', query_sig.name()[:30],
                                         query_sig.minhash.ksize)

    # gather assignments from across all the databases
    query_mins = set(query_sig.minhash.get_mins())
    n_mins = len(query_mins)

    found = []
    while 1:
        # find all of the assignments for the current set of hashes
        assignments = lca_utils.gather_assignments(query_mins, dblist)

        # none? quit.
        if not assignments:
            break

        # count 'em all
        counts = Counter()
        for hashval, assignment_set in assignments.items():
            for assignment in assignment_set:
                counts[assignment] += 1

        # find the most abundant assignment
        top_assignment, top_count = next(iter(counts.most_common()))

        # construct 'result' object
        intersect_bp = top_count * query_sig.minhash.scaled
        result = LCAGatherResult(intersect_bp = intersect_bp,
                                 f_orig_query = top_count / n_mins,
                                 f_unique_to_query=0,
                                 f_unique_weighted=0,
                                 average_abund=0,
                                 lineage=top_assignment)

        if not len(found):
            print_results("")
            print_results("overlap     p_query")
            print_results("---------   -------")

        # output!
        pct_query = '{:.1f}%'.format(result.f_orig_query*100)
        est_bp = format_bp(result.intersect_bp)
        name = format_lineage(result.lineage)

        print_results('{:9}   {:>6}      {}', est_bp, pct_query, name)

        # now, remove from query mins.
        for hashval, assignment_set in assignments.items():
            if top_assignment in assignment_set:
                query_mins.remove(hashval)

        found.append(result)

    ## done.
    
    f_unassigned = len(query_mins) / n_mins
    est_bp = len(query_mins) * query_sig.minhash.scaled

    if found:
        print_results('')
        if f_unassigned:
            print_results('{:.1f}% ({}) have no assignment.', f_unassigned*100,
                          format_bp(est_bp))
        else:
            print_results('Query is completely assigned.')
            print_results('')
    else:
        print_results('')
        print_results('No assignment for est {} of sequence.',
                      format_bp(est_bp))
        print_results('')


def gather_main(args):
    """
    @@
    """
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('db', nargs='+')
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args(args)

    if args.debug:
        set_debug(args.debug)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, None)

    # for each query, gather all the matches across databases
    query_sig = sourmash_args.load_query_signature(args.query, ksize, 'DNA')
    debug('classifying', query_sig.name())

    # make sure we're looking at the same scaled value as database
    query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

    # do the classification, output results
    gather_signature(query_sig, dblist, 0.0)


if __name__ == '__main__':
    sys.exit(gather_main(sys.argv[1:]))
