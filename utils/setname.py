#! /usr/bin/env python
from __future__ import print_function
import sys
import os, os.path
import argparse
import csv

import sourmash_lib
from sourmash_lib import signature as sig
from sourmash_lib import fig as sourmash_fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sigfiles', nargs='+')
    parser.add_argument('--name', type=str, default='')
    args = parser.parse_args()

    assert args.name

    for sigfile in args.sigfiles:
        print('setting name on "%s" to "%s"' % (sigfile, args.name))

        with open(sigfile, 'rt') as fp:
            sigs = sig.load_signatures(fp)

        for s in sigs:
            s.d['name'] = args.name

        outputname = os.path.basename(sigfile)
        with open(outputname, 'wt') as fp:
            sig.save_signatures(sigs, fp)


if __name__ == '__main__':
    main()

