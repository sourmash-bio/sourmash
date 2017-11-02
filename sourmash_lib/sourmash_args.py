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

DEFAULT_LOAD_K=31


def add_moltype_args(parser):
    parser.add_argument('--protein', dest='protein', action='store_true',
                        help='choose a protein signature (default: False)')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false',
                        help='do not choose a protein signature')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true',
                        help='choose a DNA signature (default: True)')
    parser.add_argument('--no-dna', dest='dna', action='store_false',
                        help='do not choose a DNA signature')
    parser.set_defaults(dna=None)


def add_construct_moltype_args(parser):
    parser.add_argument('--protein', dest='protein', action='store_true',
                        help='build protein signatures (default: False)')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false',
                        help='do not build protein signatures')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true',
                        help='build DNA signatures (default: True)')
    parser.add_argument('--no-dna', dest='dna', action='store_false',
                        help='do not build DNA signatures')
    parser.set_defaults(dna=True)


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


def load_query_signature(filename, ksize, select_moltype):
    sl = signature.load_signatures(filename,
                                   ksize=ksize,
                                   select_moltype=select_moltype)
    sl = list(sl)

    if len(sl) and ksize is None:
        ksizes = set([ ss.minhash.ksize for ss in sl ])
        if len(ksizes) == 1:
            ksize = ksizes.pop()
            sl = [ ss for ss in sl if ss.minhash.ksize == ksize ]
            notify('select query k={} automatically.', ksize)
        elif DEFAULT_LOAD_K in ksizes:
            sl = [ ss for ss in sl if ss.minhash.ksize == DEFAULT_LOAD_K ]
            notify('selecting default query k={}.', DEFAULT_LOAD_K)
        elif ksize:
            notify('selecting specified query k={}', ksize)

    if len(sl) != 1:
        error('When loading query from "{}"', filename)
        error('{} signatures matching ksize and molecule type;', len(sl))
        error('need exactly one. Specify --ksize or --dna/--protein.')
        sys.exit(-1)

    return sl[0]


class LoadSingleSignatures(object):
    def __init__(self, filelist,  ksize=None, select_moltype=None,
                 ignore_files=set()):
        self.filelist = filelist
        self.ksize = ksize
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
                                           ksize=self.ksize,
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


def get_ksize(tree):
    """Walk nodes in `tree` to find out ksize"""
    for node in tree.nodes.values():
        if isinstance(node, SigLeaf):
            return node.data.minhash.ksize


def load_sbts_and_sigs(filenames, query_ksize, query_moltype, traverse=False):
    n_signatures = 0
    n_databases = 0
    databases = []
    for sbt_or_sigfile in filenames:
        if traverse and os.path.isdir(sbt_or_sigfile):
            for sigfile in traverse_find_sigs([sbt_or_sigfile]):
                try:
                    siglist = sig.load_signatures(sigfile,
                                                  ksize=query_ksize,
                                                  select_moltype=query_moltype)
                    siglist = list(siglist)
                    databases.append((list(siglist), sbt_or_sigfile, False))
                    notify('loaded {} signatures from {}', len(siglist),
                           sigfile, end='\r')
                except:                       # ignore errors with traverse
                    continue
            continue
        try:
            tree = SBT.load(sbt_or_sigfile, leaf_loader=SigLeaf.load)
            ksize = get_ksize(tree)
            if ksize != query_ksize:
                error("ksize on tree '{}' is {};", sbt_or_sigfile, ksize)
                error('this is different from query ksize of {}.', query_ksize)
                sys.exit(-1)

            databases.append((tree, sbt_or_sigfile, True))
            notify('loaded SBT {}', sbt_or_sigfile, end='\r')
            n_databases += 1
        except (ValueError, EnvironmentError):
            # not an SBT - try as a .sig

            try:
                siglist = sig.load_signatures(sbt_or_sigfile,
                                              ksize=query_ksize,
                                              select_moltype=query_moltype)
                siglist = list(siglist)
                databases.append((list(siglist), sbt_or_sigfile, False))
                notify('loaded {} signatures from {}', len(siglist),
                       sbt_or_sigfile, end='\r')
                n_signatures += len(siglist)
            except EnvironmentError:
                error("\nfile '{}' does not exist", sbt_or_sigfile)
                sys.exit(-1)
    notify(' '*79, end='\r')
    notify('loaded {} signatures and {} databases total.'.format(n_signatures,
                                                                 n_databases))

    if databases:
        print('')

    return databases
