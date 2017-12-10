#! /usr/bin/env python
"""
Extract all hashes that belong to this taxonomic lineage as a signature.
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


def lca_extract_lineage(args):
    """
    main lca search function.
    """
    p = argparse.ArgumentParser()
    p.add_argument('lineage', help="semi-colon separated lineage to extract")
    p.add_argument('dblist', nargs='+')
    p.add_argument('-d', '--debug', action='store_true')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='signature file to output to')
    args = p.parse_args(args)

    if args.debug:
        set_debug(args.debug)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.dblist)
    notify('ksize={} scaled={}', ksize, scaled)

    # lineage foo
    lineage = args.lineage.split(';')
    pairs = zip(lca_utils.taxlist(), lineage)
    pairs = [ lca_utils.LineagePair(a, b) for (a, b) in pairs ]
    prefix_lineage = tuple(pairs)

    print_lineage = lca_utils.zip_lineage(prefix_lineage, truncate_empty=True)
    print_lineage = ";".join(print_lineage)

    notify('Searching for lineages that start with "{}".', print_lineage)

    keep_hashvals = set()
    for db in dblist:
        for hashval, lid_list in db.hashval_to_lineage_id.items():
            for lineage_id in lid_list:
                lineage = db.lineage_dict[lineage_id]
                lineage = lineage[:len(prefix_lineage)]
                if lineage == prefix_lineage:
                    keep_hashvals.add(hashval)
                    break

    notify('{} hashvals found that match.', len(keep_hashvals))

    if args.output:
        mh = sourmash_lib.MinHash(n=0, ksize=ksize, scaled=scaled)
        mh.add_many(keep_hashvals)
        sig = sourmash_lib.SourmashSignature(mh, name='hashes for {}'.format(print_lineage))

        notify('saving sig to {}', args.output.name)
        sourmash_lib.save_signatures([sig], args.output)


if __name__ == '__main__':
    sys.exit(lca_extract_lineage(sys.argv[1:]))
