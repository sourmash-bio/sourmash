"""
Command-line entry point for 'python -m sourmash.database'
"""
import sys
import csv
import json
import os
from collections import defaultdict

import screed
import sourmash
from sourmash.sourmash_args import FileOutput

from sourmash.logging import set_quiet, error, notify, print_results, debug
from sourmash import sourmash_args
from sourmash.minhash import _get_max_hash_for_scaled
#from sourmash.sig import _extend_signatures_with_from_file
from sourmash.picklist import SignaturePicklist, PickStyle

usage='''
sourmash database <command> [<args>] - manipulate/work with signature files.

** Commands can be:

extract <database> [<database> ... ]     - extract signature(s) from database(s)

** Use '-h' to get subcommand-specific help, e.g.

sourmash database extract -h
'''

##### actual command line functions

def create_picklist_from_args(args):
    # construct a picklist from commandline args
    picklist=None
    if args.name:
        name = str(args.name)
        picklist = SignaturePicklist('name')
        picklist.init([name])
    elif args.identprefix:
        ip = str(args.identprefix)
        picklist = SignaturePicklist("identprefix")
        picklist.init([ip])
    elif args.md5:
        md5 = str(args.md5)
        if len(md5) == 8:
            picklist = SignaturePicklist('md5prefix8')
        else:
            picklist = SignaturePicklist('md5')
        picklist.init([md5])
    return picklist



def extract(args):
    """
    extract signature(s) from database(s)
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    #sourmash.signature._extend_signatures_with_from_file(args)
    picklist = sourmash_args.load_picklist(args)
    if picklist is None:
        picklist = create_picklist_from_args(args)
    
    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()
    
    # start loading! Use manifest + picklist(s)
    progress = sourmash_args.SignatureLoadingProgress()
    loader = sourmash_args.load_file_as_signatures(args.databases,
                                                   select_moltype=moltype,
                                                   ksize=args.ksize,
                                                   picklist=picklist,
                                                   progress=progress)

    for ss in loader:
        #if filter_fn(ss):
        save_sigs.add(ss)

    notify(f"loaded {len(progress)} total that matched ksize & molecule type")
    if not save_sigs:
        error("no matching signatures to save!")
        sys.exit(-1)

    save_sigs.close()

    notify(f"extracted {len(save_sigs)} signatures from {len(args.databases)} file(s)")

    if picklist:
        sourmash_args.report_picklist(args, picklist)

def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.database, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)