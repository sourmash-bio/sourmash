#! /usr/bin/env python
"""
Check SBT search by taking every leaf node in a tree and checking to make
sure we can find it.
"""
import argparse
import sourmash
from sourmash.sbtmh import search_minhashes

THRESHOLD=0.08


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sbt')
    args = p.parse_args()

    db = sourmash.sbtmh.load_sbt_index(args.sbt)
    threshold = THRESHOLD

    for leaf in db.leaves():
        query = leaf.data
        matches = db.find(search_minhashes, query, threshold)
        matches = list([ x.data for x in matches ])
        if query not in matches:
            print(query)
            assert 0
                                                 

if __name__ == '__main__':
    main()
