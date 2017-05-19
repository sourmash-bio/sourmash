"""
Utility functions for dealing with input args to the sourmash command line.
"""
import sys
import os
from . import signature
from .logging import notify, error

from . import signature as sig
from sourmash_lib.sbt import SBT
from sourmash_lib.sbtmh import SigLeaf


def add_moltype_args(parser, default_dna=None):
    parser.add_argument('--protein', dest='protein', action='store_true')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true')
    parser.add_argument('--no-dna', dest='dna', action='store_false')
    parser.set_defaults(dna=default_dna)


def add_ksize_arg(parser, default):
    parser.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: {d})'.format(d=default))


def get_moltype(sig, require=False):
    if sig.minhash.is_molecule_type('DNA'):
        moltype = 'DNA'
    elif sig.minhash.is_molecule_type('protein'):
        moltype = 'protein'
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig.name()))

    return moltype


def calculate_moltype(args, default=None):
    if args.protein:
        if args.dna is True:
            error('cannot specify both --dna and --protein!')
            sys.exit(-1)
        args.dna = False

    moltype = default
    if args.protein:
        moltype = 'protein'
    elif args.dna:
        moltype = 'DNA'

    return moltype


def load_query_signature(filename, select_ksize, select_moltype):
    print(select_ksize, select_moltype)
    sl = signature.load_signatures(filename,
                                   select_ksize=select_ksize,
                                   select_moltype=select_moltype)
    sl = list(sl)

    if len(sl) != 1:
        error('When loading query from "{}"', filename)
        error('{} signatures matching ksize and molecule type;', len(sl))
        error('need exactly one. Specify --ksize or --dna/--protein.')
        sys.exit(-1)

    return sl[0]


class LoadSingleSignatures(object):
    def __init__(self, filelist,  select_ksize=None, select_moltype=None,
                 ignore_files=set()):
        self.filelist = filelist
        self.select_ksize = select_ksize
        self.select_moltype = select_moltype
        self.ignore_files = ignore_files

        self.skipped_ignore = 0
        self.skipped_nosig = 0
        self.ksizes = set()
        self.moltypes = set()

    def __iter__(self):
        for filename in self.filelist:
            if filename in self.ignore_files:
                self.skipped_ignore += 1
                continue

            sl = signature.load_signatures(filename,
                                           select_ksize=self.select_ksize,
                                           select_moltype=self.select_moltype)
            sl = list(sl)
            if len(sl) == 0:
                self.skipped_nosig += 1
                continue

            for query in sl:
                query_moltype = get_moltype(query)
                query_ksize = query.minhash.ksize

                self.ksizes.add(query_ksize)
                self.moltypes.add(query_moltype)

                yield filename, query, query_moltype, query_ksize

            if len(self.ksizes) > 1 or len(self.moltypes) > 1:
                raise ValueError('multiple k-mer sizes/molecule types present')


def traverse_find_sigs(dirnames):
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig'):
                    fullname = os.path.join(root, name)
                    yield fullname


def load_sbts_and_sigs(filenames, query_ksize, query_moltype):
    databases = []
    for sbt_or_sigfile in filenames:
        try:
            tree = SBT.load(sbt_or_sigfile, leaf_loader=SigLeaf.load)
            databases.append((tree, True))
            notify('loaded SBT {}', sbt_or_sigfile)
        except (ValueError, FileNotFoundError):
            # not an SBT - try as a .sig

            try:
                siglist = sig.load_signatures(sbt_or_sigfile,
                                              select_ksize=query_ksize,
                                              select_moltype=query_moltype)
                siglist = list(siglist)
                databases.append((list(siglist), False))
                notify('loaded {} signatures from {}', len(siglist),
                       sbt_or_sigfile)
            except:
                raise

    return databases
