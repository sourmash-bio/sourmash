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
import sourmash.exceptions

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


def load_query_signature(filename, ksize, select_moltype, select_md5=None):
    try:
        sl = load_file_as_signatures(filename, ksize=ksize,
                                     select_moltype=select_moltype)
        sl = list(sl)
    except (OSError, ValueError):
        error("Cannot open file '{}'", filename)
        sys.exit(-1)

    if len(sl) and select_md5:
        for sig in sl:
            sig_md5 = sig.md5sum()
            if sig_md5.startswith(select_md5.lower()):
                sl = [sig]
                break

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


def traverse_find_sigs(filenames, yield_all_files=False):
    for filename in filenames:
        if os.path.isfile(filename) and \
                  (filename.endswith('.sig') or yield_all_files):
            yield filename
            continue

        # filename is a directory --
        dirname = filename

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


def check_lca_db_is_compatible(filename, db, query):
    query_mh = query.minhash
    if db.ksize != query_mh.ksize:
        error("ksize on db '{}' is {};", filename, db.ksize)
        error('this is different from query ksize of {}.', query_mh.ksize)
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
    for filename in filenames:
        notify('loading from {}...', filename, end='\r')

        try:
            db, dbtype = _load_database(filename, traverse=traverse)
        except IOError as e:
            notify(str(e))
            sys.exit(-1)

        # are we collecting signatures from a directory/path?
        # NOTE: error messages about loading will now be attributed to
        # directory, not individual file.
        if traverse and os.path.isdir(filename):
            assert dbtype == DatabaseType.SIGLIST

            siglist = _select_sigs(db, moltype=query_moltype, ksize=query_ksize)
            siglist = filter_compatible_signatures(query, siglist, 1)
            linear = LinearIndex(siglist, filename=filename)
            databases.append((linear, filename, False))

            n_signatures += len(linear)

        # SBT
        elif dbtype == DatabaseType.SBT:
            if not check_tree_is_compatible(filename, db, query,
                                            is_similarity_query):
                sys.exit(-1)

            databases.append((db, filename, 'SBT'))
            notify('loaded SBT {}', filename, end='\r')
            n_databases += 1

        # LCA
        elif dbtype == DatabaseType.LCA:
            if not check_lca_db_is_compatible(filename, db, query):
                sys.exit(-1)
            query_scaled = query.minhash.scaled

            notify('loaded LCA {}', filename, end='\r')
            n_databases += 1

            databases.append((db, filename, 'LCA'))

        # signature file
        elif dbtype == DatabaseType.SIGLIST:
            siglist = _select_sigs(db, moltype=query_moltype, ksize=query_ksize)
            siglist = list(siglist)
            if not siglist:          # file not found, or parse error?
                raise ValueError

            siglist = filter_compatible_signatures(query, siglist, False)
            linear = LinearIndex(siglist, filename=filename)
            databases.append((linear, filename, 'signature'))

            notify('loaded {} signatures from {}', len(linear),
                   filename, end='\r')
            n_signatures += len(linear)

        # unknown!?
        else:
            raise Exception("unknown dbtype {}".format(dbtype))

        # END for loop


    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.', n_signatures, 
                                                               n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        notify('** ERROR: no signatures or databases loaded?')
        sys.exit(-1)

    if databases:
        print('')

    return databases


class DatabaseType(Enum):
    SIGLIST = 1
    SBT = 2
    LCA = 3


def _load_database(filename, traverse=False, traverse_yield_all=False):
    """Load file as a database - list of signatures, LCA, SBT, etc.

    Return (db, dbtype), where dbtype is a DatabaseType enum.

    This is an internal function used by other functions in sourmash_args.
    """
    loaded = False
    dbtype = None

    # special case stdin
    if not loaded and filename == '-':
        db = sourmash.load_signatures(sys.stdin, quiet=True, do_raise=True)
        db = list(db)
        loaded = True
        dbtype = DatabaseType.SIGLIST

    # load signatures from directory
    if not loaded and os.path.isdir(filename) and traverse:
        all_sigs = []
        for thisfile in traverse_find_sigs([filename], traverse_yield_all):
            try:
                with open(thisfile, 'rt') as fp:
                    x = sourmash.load_signatures(fp, quiet=True, do_raise=True)
                    siglist = list(x)
                    all_sigs.extend(siglist)
            except (IOError, sourmash.exceptions.SourmashError):
                if traverse_yield_all:
                    continue
                else:
                    raise

        loaded=True
        db = all_sigs
        dbtype = DatabaseType.SIGLIST

    # load signatures from single file
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
        raise OSError("Error while reading signatures from '{}'.".format(filename))

    return db, dbtype


# note: dup from index.py internal function.
def _select_sigs(siglist, ksize, moltype):
    for ss in siglist:
        if (ksize is None or ss.minhash.ksize == ksize) and \
           (moltype is None or ss.minhash.moltype == moltype):
           yield ss


def load_file_as_index(filename, traverse=True, yield_all_files=False):
    "Load 'filename' as an Index class; generic database loader."
    db, dbtype = _load_database(filename, traverse,
                                traverse_yield_all=yield_all_files)
    if dbtype in (DatabaseType.LCA, DatabaseType.SBT):
        return db                         # already an index!
    elif dbtype == DatabaseType.SIGLIST:
        # turn siglist into a LinearIndex
        idx = LinearIndex(db, filename)
        return idx
    else:
        assert 0                          # unknown enum!?


def load_file_as_signatures(filename, select_moltype=None, ksize=None,
                            traverse=False, yield_all_files=False):
    """Load 'filename' as a collection of signatures. Return an iterable.

    If it's an LCA or SBT, call the .signatures() method on it.

    Applies selector function if select_moltype, ksize are given.
    """
    db, dbtype = _load_database(filename, traverse=traverse,
                                traverse_yield_all=yield_all_files)

    if dbtype in (DatabaseType.LCA, DatabaseType.SBT):
        db = db.select(moltype=select_moltype, ksize=ksize)
        return db.signatures()
    elif dbtype == DatabaseType.SIGLIST:
        return list(_select_sigs(db, moltype=select_moltype, ksize=ksize))
    else:
        assert 0                          # unknown enum!?


def load_file_list_of_signatures(filename):
    with open(filename, 'rt') as fp:
        file_list = [ x.rstrip('\r\n') for x in fp ]

    return file_list


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
