"""
Utility functions for dealing with input args to the sourmash command line.
"""
import sys
import os
import argparse
from . import signature
from .logging import notify, error

from .index import LinearIndex
from . import signature as sig
from .sbt import SBT
from .sbtmh import SigLeaf
from .lca import lca_utils
import sourmash

DEFAULT_LOAD_K = 31
DEFAULT_N = 500


def get_moltype(sig, require=False):
    if sig.minhash.is_molecule_type('DNA'):
        moltype = 'DNA'
    elif sig.minhash.is_molecule_type('dayhoff'):
        moltype = 'dayhoff'
    elif sig.minhash.is_molecule_type('hp'):
        moltype = 'hp'
    elif sig.minhash.is_molecule_type('protein'):
        moltype = 'protein'
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig.name()))

    return moltype


def calculate_moltype(args, default=None):
    if args.protein:
        if args.dna is True:
            error('cannot specify both --dna/--rna and --protein!')
            sys.exit(-1)
        args.dna = False

    moltype = default
    if args.dna:
        moltype = 'DNA'
    elif args.dayhoff:
        moltype = 'dayhoff'
    elif args.hp:
        moltype = 'hp'
    elif args.protein:
        moltype = 'protein'

    return moltype


def load_query_signature(filename, ksize, select_moltype):
    try:
        sl = signature.load_signatures(filename,
                                       ksize=ksize,
                                       select_moltype=select_moltype,
                                       do_raise=True)
        sl = list(sl)
    except (IOError, ValueError):
        error("Cannot open file '{}'", filename)
        sys.exit(-1)

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
        error('need exactly one. Specify --ksize or --dna, --rna, or --protein.')
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

            if len(self.ksizes) > 1 or len(self.moltypes) > 1:
                raise ValueError('multiple k-mer sizes/molecule types present')

            for query in sl:
                yield filename, query, query_moltype, query_ksize


def traverse_find_sigs(dirnames, yield_all_files=False):
    for dirname in dirnames:
        if (dirname.endswith('.sig') or yield_all_files) and os.path.isfile(dirname):
            yield dirname
            continue

        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig') or yield_all_files:
                    fullname = os.path.join(root, name)
                    yield fullname


def filter_compatible_signatures(query, siglist, force=False):
    for ss in siglist:
        if check_signatures_are_compatible(query, ss):
            yield ss
        else:
            if not force:
                raise ValueError("incompatible signature")


def check_signatures_are_compatible(query, subject):
    # is one scaled, and the other not? cannot do search
    if query.minhash.scaled and not subject.minhash.scaled or \
       not query.minhash.scaled and subject.minhash.scaled:
       error("signature {} and {} are incompatible - cannot compare.",
             query.name(), subject.name())
       if query.minhash.scaled:
           error("{} was calculated with --scaled, {} was not.",
                 query.name(), subject.name())
       if subject.minhash.scaled:
           error("{} was calculated with --scaled, {} was not.",
                 subject.name(), query.name())
       return 0

    return 1

def get_SBT_info(tree):
    # get a minhash from the tree
    leaf = next(iter(tree.leaves()))
    tree_ss = leaf.data
    tree_mh = tree_ss.minhash

    moltype = get_moltype(tree_ss)

    return tree_mh.num, tree_mh.scaled, tree_mh.ksize, moltype

def check_tree_is_compatible(treename, tree, query, is_similarity_query):
    # get a minhash from the tree
    leaf = next(iter(tree.leaves()))
    tree_mh = leaf.data.minhash

    query_mh = query.minhash

    if tree_mh.ksize != query_mh.ksize:
        error("ksize on tree '{}' is {};", treename, tree_mh.ksize)
        error('this is different from query ksize of {}.', query_mh.ksize)
        return 0

    # is one scaled, and the other not? cannot do search.
    if (tree_mh.scaled and not query_mh.scaled) or \
       (query_mh.scaled and not tree_mh.scaled):
        error("for tree '{}', tree and query are incompatible for search.",
              treename)
        if tree_mh.scaled:
            error("tree was calculated with scaled, query was not.")
        else:
            error("query was calculated with scaled, tree was not.")
        return 0

    # are the scaled values incompatible? cannot downsample tree for similarity
    if tree_mh.scaled and tree_mh.scaled < query_mh.scaled and \
      is_similarity_query:
        error("for tree '{}', scaled value is smaller than query.", treename)
        error("tree scaled: {}; query scaled: {}. Cannot do similarity search.",
              tree_mh.scaled, query_mh.scaled)
        return 0

    return 1


class SearchDatabaseLoader(object):
    def __init__(self, filenames, is_similarity_query, traverse):
        self.filenames = filenames
        self.is_similarity_query = is_similarity_query
        self.traverse = traverse
        self.databases = []

        self.ksizes = set()
        self.moltypes = set()
        self.num_db = set()
        self.scaled_db = set()

    def load_list_of_signatures(self, filename):
        # are we collecting signatures from a directory/path?
        if not (self.traverse and os.path.isdir(filename)):
            return False

        for sigfile in traverse_find_sigs([filename]):
            try:
                siglist = sig.load_signatures(sigfile)
                #siglist = filter_compatible_signatures(query, siglist, 1)
                linear = LinearIndex(siglist, filename=sigfile)
                self.databases.append((linear, filename, 'linear'))
                notify('loaded {} signatures from {}', len(linear),
                       sigfile, end='\r')
                n_signatures += len(linear)
            except Exception:        # ignore errors with traverse
                pass

        return True

    def load_SBT(self, filename):
        # are we loading an SBT?
        try:
            tree = SBT.load(filename, leaf_loader=SigLeaf.load)

            is_num, is_scaled, ksize, moltype = get_SBT_info(tree)
            if is_num:
                self.num_db.add(filename)
            if is_scaled:
                self.scaled_db.add(filename)
            self.ksizes.add(ksize)
            self.moltypes.add(moltype)

            self.databases.append((tree, filename, 'SBT'))
            notify('loaded SBT {}', filename, end='\r')
        except (ValueError, EnvironmentError):
            # not an SBT
            return False

        return True

    def load_LCA(self, filename):
        try:
            lca_db = lca_utils.LCA_Database()
            lca_db.load(filename)

            # LCA databases are always made from --scaled signatures
            self.scaled_db.add(filename)
            self.ksizes.add(lca_db.ksize)
            # LCA databases are always DNA
            self.moltypes.add('DNA')

            notify('loaded LCA {}', filename, end='\r')

            self.databases.append((lca_db, filename, 'LCA'))
        except (ValueError, TypeError, EnvironmentError):
            # not an LCA database
            return False

        return True

    def load_signature_file(self, filename):
        try:
            siglist = sig.load_signatures(filename)
            siglist = list(siglist)
            if len(siglist) == 0:         # file not found, or parse error?
                raise ValueError

            #siglist = filter_compatible_signatures(query, siglist, False)
            linear = LinearIndex(siglist, filename=filename)
            self.databases.append((linear, filename, 'linear'))

            notify('loaded {} signatures from {}', len(linear), filename,
                   end='\r')
        except (EnvironmentError, ValueError):
            return False

        return True

    def load_all(self):
        for filename in self.filenames:
            notify('loading from {}...', filename, end='\r')

            loaded = self.load_list_of_signatures(filename)
            if not loaded:
                loaded = self.load_SBT(filename)
            if not loaded:
                loaded = self.load_LCA(filename)
            if not loaded:
                loaded = self.load_signature_file(filename)

            if not loaded:
                error("\nCannot load signatures from file '{}'", filename)
                sys.exit(-1)

            if self.num_db and self.scaled_db:
                if filename in self.num_db:
                    assert self.scaled_db
                    assert len(self.num_db) == 1
                    error('file {} is incompatible with the other databases',
                          filename)
                    error('the signatures in it were calculated with --num')
                    error('the other databases have --scaled signatures')
                elif filename in self.scaled_db:
                    assert self.num_db
                    assert len(self.scaled_db) == 1
                    error('file {} is incompatible with the other databases',
                          filename)
                    error('the signatures in it were calculated with --scaled')
                    error('the other databases have --num signatures')
                else:
                    assert 0

    def filter_signatures(self):
        assert (self.scaled_db and not self.num_db) or \
               (self.num_db and not self.scaled_db)
        assert len(self.ksizes) == 1
        assert len(self.moltypes) == 1

        ksize = list(self.ksizes)[0]
        moltype = list(self.moltypes)[0]
        is_scaled = bool(self.scaled_db)

        for (obj, filename, objtype) in self.databases:
            if objtype == 'linear':
                siglist = []
                for sig in obj.signatures():
                    if sig.minhash.ksize != ksize:
                        continue
                    if get_moltype(sig) != moltype:
                        continue
                    if sig.minhash.scaled and not is_scaled:
                        continue
                    if not sig.minhash.scaled and is_scaled:
                        continue

                    siglist.append(sig)

                if not siglist:
                    error("No compatible signatures in {}", filename)
                    sys.exit(-1)

                # @CTB we shouldn't reach deeply into the guts here :)
                obj._signatures = siglist

    def summarize_files(self):
        n_signatures = 0
        n_databases = 0

        for (obj, filename, objtype) in self.databases:
            if objtype == 'linear':
                n_signatures += len(obj)
            else:
                assert objtype in ('SBT', 'LCA')
                n_databases += 1

        return n_signatures, n_databases




def load_dbs_and_sigs(filenames, query, is_similarity_query, traverse=False):
    """
    Load one or more SBTs, LCAs, and/or signatures.

    Check for compatibility with query.
    """
    loader = SearchDatabaseLoader(filenames, is_similarity_query, traverse)

    # this loads all the databases and signatures without filtering
    loader.load_all()

    ksize = query.minhash.ksize
    scaled = query.minhash.scaled
    moltype = get_moltype(query)

    # did we load any databases? if so, check against query values. @CTB
    if loader.ksizes:
        assert loader.moltypes
        assert loader.scaled_db or loader.num_db

        if ksize not in loader.ksizes:
            assert 0
        if moltype not in loader.moltypes:
            assert 0
        if scaled and not loader.scaled_db:
            assert 0
        if not scaled and not loader.num_db:
            assert 0
    else:                                 # no databases loaded. @CTB
        loader.ksizes.add(ksize)
        loader.moltypes.add(moltype)
        if scaled:
            loader.scaled_db.add('foo')
        else:
            loader.num_db.add('foo')

    # now, see if we can filter the signatures down to match.
    loader.filter_signatures()

    n_signatures, n_databases = loader.summarize_files()

    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.',
                    n_signatures, n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        sys.exit(-1)

    return loader.databases


class FileOutput(object):
    """A context manager for file outputs that handles sys.stdout gracefully.

    Usage:

       with FileOutput(filename, mode) as fp:
          ...

    does what you'd expect, but it handles the situation where 'filename'
    is '-' or None. This makes it nicely compatible with argparse usage,
    e.g.

    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    ...
    with FileOutput(args.output, 'wt') as fp:
       ...

    will properly handle no argument or '-' as sys.stdout.
    """
    def __init__(self, filename, mode='wt'):
        self.filename = filename
        self.mode = mode
        self.fp = None

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        self.fp = open(self.filename, self.mode)
        return self.fp

    def __enter__(self):
        return self.open()

    def __exit__(self, type, value, traceback):
        # do we need to handle exceptions here?
        if self.fp:
            self.fp.close()

        return False
