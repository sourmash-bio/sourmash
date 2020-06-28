"""
Utility functions for dealing with input args to the sourmash command line.
"""
import sys
import os
import argparse
import itertools
from enum import Enum

from sourmash import load_sbt_index
from sourmash.lca.lca_db import load_single_database

from . import signature
from .logging import notify, error

from .index import LinearIndex
from . import signature as sig
from .sbt import SBT
from .sbtmh import SigLeaf
from .lca import LCA_Database
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
    moltype = default

    n = 0
    if args.dna:
        moltype = 'DNA'
        n += 1
    if args.dayhoff:
        moltype = 'dayhoff'
        n += 1
    if args.hp:
        moltype = 'hp'
        n += 1
    if args.protein:
        moltype = 'protein'
        n += 1

    if n > 1:
        error("cannot specify more than one of --dna/--rna/--protein/--hp/--dayhoff")
        sys.exit(-1)

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


def load_dbs_and_sigs(filenames, query, is_similarity_query, traverse=False):
    """
    Load one or more SBTs, LCAs, and/or signatures.

    Check for compatibility with query.
    """
    query_ksize = query.minhash.ksize
    query_moltype = get_moltype(query)

    n_signatures = 0
    n_databases = 0
    databases = []
    for sbt_or_sigfile in filenames:
        notify('loading from {}...', sbt_or_sigfile, end='\r')

        db, dbtype = load_database(sbt_or_sigfile)

        # are we collecting signatures from a directory/path?
        if 0 and traverse and os.path.isdir(sbt_or_sigfile):
            for sigfile in traverse_find_sigs([sbt_or_sigfile]):
                try:
                    siglist = sig.load_signatures(sigfile,
                                                  ksize=query_ksize,
                                                  select_moltype=query_moltype)
                    siglist = filter_compatible_signatures(query, siglist, 1)
                    linear = LinearIndex(siglist, filename=sigfile)
                    databases.append((linear, sbt_or_sigfile, False))
                    notify('loaded {} signatures from {}', len(linear),
                           sigfile, end='\r')
                    n_signatures += len(linear)
                except Exception:                       # ignore errors with traverse
                    pass

            # done! jump to beginning of main 'for' loop
            continue

        if dbtype == DatabaseType.SBT:
            tree = db
            if not check_tree_is_compatible(sbt_or_sigfile, tree, query,
                                            is_similarity_query):
                sys.exit(-1)

            databases.append((tree, sbt_or_sigfile, 'SBT'))
            notify('loaded SBT {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            # done! jump to beginning of main 'for' loop
            continue

        if dbtype == DatabaseType.LCA:
            lca_db = db

            assert query_ksize == lca_db.ksize
            query_scaled = query.minhash.scaled

            notify('loaded LCA {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            databases.append((lca_db, sbt_or_sigfile, 'LCA'))

            continue

        if dbtype == DatabaseType.SIGLIST:
            siglist = _select_sigs(db, moltype=query_moltype, ksize=query_ksize)
            siglist = list(siglist)
            if len(siglist) == 0:         # file not found, or parse error?
                raise ValueError

            siglist = filter_compatible_signatures(query, siglist, False)
            linear = LinearIndex(siglist, filename=sbt_or_sigfile)
            databases.append((linear, sbt_or_sigfile, 'signature'))

            notify('loaded {} signatures from {}', len(linear),
                   sbt_or_sigfile, end='\r')
            n_signatures += len(linear)

            continue

    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.', n_signatures, 
                                                               n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        sys.exit(-1)

    if databases:
        print('')

    return databases


class DatabaseType(Enum):
    SIGLIST = 1
    SBT = 2
    LCA = 3


def load_database(filename):
    """Load file as a database - list of signatures, LCA, SBT, etc.

    Return (db, dbtype), where dbtype is a DatabaseType enum.

    This will (eventually) supersede load_dbs_and_sigs.

    TODO:
    - add traversal behavior + force load for directories.
    - add stdin for reading signatures?
    - maybe add file lists?
    """
    loaded = False
    dbtype = None
    try:
        # CTB: could make this a generator, with some trickery; but for
        # now, just force into list.
        with open(filename, 'rt') as fp:
            db = sourmash.load_signatures(fp, quiet=True, do_raise=True)
            db = list(db)

        loaded = True
        dbtype = DatabaseType.SIGLIST
    except Exception as exc:
        pass

    if not loaded:                    # try load as SBT
        try:
            db = load_sbt_index(filename)
            loaded = True
            dbtype = DatabaseType.SBT
        except:
            pass

    if not loaded:                    # try load as LCA
        try:
            db, _, _ = load_single_database(filename)
            loaded = True
            dbtype = DatabaseType.LCA
        except:
            pass

    if not loaded:
        error('\nError while reading signatures from {}.'.format(filename))
        sys.exit(-1)

    return db, dbtype


# note: dup from index.py internal function.
def _select_sigs(siglist, ksize, moltype):
    for ss in siglist:
        if (ksize is None or ss.minhash.ksize == ksize) and \
           (moltype is None or ss.minhash.moltype == moltype):
           yield ss


def load_file_as_signatures(filename, select_moltype=None, ksize=None):
    """Load 'filename' as a collection of signatures. Return an iterable.

    If it's an LCA or SBT, call the .signatures() method on it.

    Applies selector function if select_moltype, ksize are given.
    """
    db, dbtype = load_database(filename)
    if dbtype in (DatabaseType.LCA, DatabaseType.SBT):
        db = db.select(moltype=select_moltype, ksize=ksize)
        return db.signatures()
    elif dbtype == DatabaseType.SIGLIST:
        return list(_select_sigs(db, moltype=select_moltype, ksize=ksize))


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
