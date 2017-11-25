#! /usr/bin/env python
"""
Summarize the taxonomic content of the given signatures, combined.
"""
from __future__ import print_function
import sys
import argparse
import csv
from collections import defaultdict, Counter

import sourmash_lib
from sourmash_lib import sourmash_args
from sourmash_lib.logging import notify, error, print_results
from sourmash_lib.lca import lca_utils
from sourmash_lib.lca.lca_utils import debug, set_debug


DEFAULT_THRESHOLD=5


def lca_search_main(args):
    """
    main lca search function.
    """
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('dblist', nargs='+')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='CSV output')
    p.add_argument('--scaled', type=float)
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-k', '--ksize', type=int, default=31)
    args = p.parse_args(args)

    if args.debug:
        set_debug(args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load query
    query_sig = sourmash_lib.load_one_signature(args.query, ksize=args.ksize)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.dblist, args.scaled)
    notify('ksize={} scaled={}', ksize, scaled)
    assert ksize == args.ksize

    # get the correct query downsampling for the database
    query_mh = query_sig.minhash
    query_mh = query_mh.downsample_scaled(scaled)

    # for each hash in query, gather all matches
    counts = Counter()
    for hashval in query_mh.get_mins():
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            if lineages:
                lineages = tuple(sorted(lineages))
                counts[lineages] += 1
    debug(counts)

    total = len(query_mh.get_mins())

    for lineages, count in counts.most_common():
        if count < args.threshold:
            break

        if not lineages:
            display = ['(root)']
        else:
            display = []
            for lineage in lineages:
                if lineage:
                    lineage = lca_utils.zip_lineage(lineage,
                                                    truncate_empty=True)
                    lineage = ';'.join(lineage)
                else:
                    lineage = '(root)'
                    assert 0
                display.append(lineage)

        p = count / total * 100.
        p = '{:.1f}%'.format(p)

        print_results('{:5} {:>5}   {}', p, count, display[0])
        for d in display[1:]:
            print_results('              {}', d)


if __name__ == '__main__':
    sys.exit(lca_search_main(sys.argv[1:]))
