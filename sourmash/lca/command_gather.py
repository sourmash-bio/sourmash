#! /usr/bin/env python
"""
Execute a greedy search on lineages attached to hashvals in the query.

Mimics `sourmash gather` but provides taxonomic information.
"""
from __future__ import print_function, division
import sys
import argparse
import csv
from collections import Counter, defaultdict, namedtuple

from .. import sourmash_args, save_signatures, SourmashSignature
from ..logging import notify, error, print_results
from . import lca_utils
from .lca_utils import debug, set_debug, check_files_exist
from ..search import format_bp

LCAGatherResult = namedtuple('LCAGatherResult',
                             'intersect_bp, f_unique_to_query, f_unique_weighted, average_abund, lineage, f_match, name, n_equal_matches')


def format_lineage(lineage_tup):
    """
    Pretty print lineage.
    """
    # list of ranks present
    present = [ l.rank for l in lineage_tup if l.name ]
    d = dict(lineage_tup) # rank: value

    if 'genus' in present:
        genus = d['genus']
        if 'strain' in present:
            name = d['strain']
        elif 'species' in present:
            species = d['species']
            if species.startswith(genus + ' ') or \
              species.startswith(genus + '_'):
                name = species
            else:
                name = '{} {}'.format(genus, species)
        else:
            name = '{} sp.'.format(genus)
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


def gather_signature(query_sig, dblist, ignore_abundance):
    """
    Decompose 'query_sig' using the given list of databases.
    """
    notify('loaded query: {}... (k={})', query_sig.name()[:30],
                                         query_sig.minhash.ksize)

    # extract the basic set of mins
    query_mins = set(query_sig.minhash.get_mins())
    n_mins = len(query_mins)

    if query_sig.minhash.track_abundance and not ignore_abundance:
        orig_abunds = query_sig.minhash.get_mins(with_abundance=True)
    else:
        if query_sig.minhash.track_abundance and ignore_abundance:
            notify('** ignoring abundance')
        orig_abunds = { k: 1 for k in query_mins }
    sum_abunds = sum(orig_abunds.values())


    # collect all mentioned lineage_ids -> md5s, from across the databases
    md5_to_lineage = {}
    md5_to_name = {}

    x = set()
    for hashval in query_mins:
        for lca_db in dblist:
            lineage_ids = lca_db.hashval_to_lineage_id.get(hashval, [])
            for lid in lineage_ids:
                md5 = lca_db.lineage_id_to_signature[lid]
                x.add((lca_db, lid, md5))

    for lca_db, lid, md5 in x:
        md5_to_lineage[md5] = lca_db.lineage_dict[lid]
        if lca_db.signature_to_name:
            md5_to_name[md5] = lca_db.signature_to_name[md5]
        else:
            md5_to_name[md5] = ''

    # now! do the gather:
    while 1:
        # find all of the assignments for the current set of hashes
        assignments = defaultdict(set)
        for hashval in query_mins:
            for lca_db in dblist:
                lineage_ids = lca_db.hashval_to_lineage_id.get(hashval, [])
                for lid in lineage_ids:
                    md5 = lca_db.lineage_id_to_signature[lid]
                    signature_size = lca_db.lineage_id_counts[lid]
                    assignments[hashval].add((md5, signature_size))
        # none? quit.
        if not assignments:
            break

        # count the distinct signatures.
        counts = Counter()
        for hashval, md5_set in assignments.items():
            for (md5, sigsize) in md5_set:
                counts[(md5, sigsize)] += 1

        # collect the most abundant assignments
        common_iter = iter(counts.most_common())
        best_list = []
        (md5, sigsize), top_count = next(common_iter)

        best_list.append((md5_to_name[md5], md5, sigsize))
        for (md5, sigsize), count in common_iter:
            if count != top_count:
                break
            best_list.append((md5_to_name[md5], md5, sigsize))

        # sort on name and pick the top (for consistency).
        best_list.sort()
        _, top_md5, top_sigsize = best_list[0]

        equiv_counts = len(best_list) - 1

        # now, remove from query mins.
        intersect_mins = set()
        for hashval, md5_set in assignments.items():
            if (top_md5, top_sigsize) in md5_set:
                query_mins.remove(hashval)
                intersect_mins.add(hashval)

        # should match!
        assert top_count == len(intersect_mins)

        # calculate size of match (# of hashvals belonging to that sig)
        match_size = top_sigsize

        # construct 'result' object
        intersect_bp = top_count * query_sig.minhash.scaled
        f_unique_weighted = sum((orig_abunds[k] for k in intersect_mins)) \
               / sum_abunds
        average_abund = sum((orig_abunds[k] for k in intersect_mins)) \
               / len(intersect_mins)
        f_match = len(intersect_mins) / match_size

        result = LCAGatherResult(intersect_bp = intersect_bp,
                                 f_unique_to_query= top_count / n_mins,
                                 f_unique_weighted=f_unique_weighted,
                                 average_abund=average_abund,
                                 f_match=f_match,
                                 lineage=md5_to_lineage[top_md5],
                                 name=md5_to_name[top_md5],
                                 n_equal_matches=equiv_counts)

        f_unassigned = len(query_mins) / n_mins
        est_bp = len(query_mins) * query_sig.minhash.scaled

        yield result, f_unassigned, est_bp, query_mins

    ## done.


def gather_main(args):
    """
    Do a greedy search for the hash components of a query against an LCA db.

    Here we don't actually do a least-common-ancestor search of any kind; we
    do essentially the same kind of search as we do in `sourmash gather`, with
    the main difference that we are implicitly combining different genomes of
    identical lineages.

    This takes advantage of the structure of the LCA db, where we store the
    full lineage information for each known hash, as opposed to storing only
    the least-common-ancestor information for it.
    """
    p = argparse.ArgumentParser(prog="sourmash lca gather")
    p.add_argument('query')
    p.add_argument('db', nargs='+')
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV containing matches to this file')
    p.add_argument('--output-unassigned', type=argparse.FileType('wt'),
                   help='output unassigned portions of the query as a signature to this file')
    p.add_argument('--ignore-abundance',  action='store_true',
                   help='do NOT use k-mer abundances if present')
    args = p.parse_args(args)

    if args.debug:
        set_debug(args.debug)

    if not check_files_exist(args.query, *args.db):
        sys.exit(-1)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, None)

    # for each query, gather all the matches across databases
    query_sig = sourmash_args.load_query_signature(args.query, ksize, 'DNA')
    debug('classifying', query_sig.name())

    # make sure we're looking at the same scaled value as database
    query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

    # do the classification, output results
    found = []
    for result, f_unassigned, est_bp, remaining_mins in gather_signature(query_sig, dblist, args.ignore_abundance):
        # is this our first time through the loop? print headers, if so.
        if not len(found):
            print_results("")
            print_results("overlap     p_query p_match ")
            print_results("---------   ------- --------")

        # output!
        pct_query = '{:.1f}%'.format(result.f_unique_to_query*100)
        pct_match = '{:.1f}%'.format(result.f_match*100)
        str_bp = format_bp(result.intersect_bp)
        name = format_lineage(result.lineage)

        equal_match_str = ""
        if result.n_equal_matches:
            equal_match_str = " (** {} equal matches)".format(result.n_equal_matches)

        print_results('{:9}   {:>6}  {:>6}      {}{}', str_bp, pct_query,
                      pct_match, name, equal_match_str)

        found.append(result)

    if found:
        print_results('')
        if f_unassigned:
            print_results('{:.1f}% ({}) of hashes have no assignment.', f_unassigned*100,
                          format_bp(est_bp))
        else:
            print_results('Query is completely assigned.')
            print_results('')
    # nothing found.
    else:
        est_bp = len(query_sig.minhash.get_mins()) * query_sig.minhash.scaled
        print_results('')
        print_results('No assignment for est {} of sequence.',
                      format_bp(est_bp))
        print_results('')

    if not found:
        sys.exit(0)

    if args.output:
        fieldnames = ['intersect_bp', 'f_match', 'f_unique_to_query', 'f_unique_weighted',
                      'average_abund', 'name', 'n_equal_matches'] + list(lca_utils.taxlist())

        w = csv.DictWriter(args.output, fieldnames=fieldnames)
        w.writeheader()
        for result in found:
            lineage = result.lineage
            d = dict(result._asdict())
            del d['lineage']

            for (rank, value) in lineage:
                d[rank] = value

            w.writerow(d)

    if args.output_unassigned:
        if not found:
            notify('nothing found - entire query signature unassigned.')
        elif not remaining_mins:
            notify('no unassigned hashes! not saving.')
        else:
            outname = args.output_unassigned.name
            notify('saving unassigned hashes to "{}"', outname)

            e = query_sig.minhash.copy_and_clear()
            e.add_many(remaining_mins)

            save_signatures([ SourmashSignature(e) ], args.output_unassigned)


if __name__ == '__main__':
    sys.exit(gather_main(sys.argv[1:]))
