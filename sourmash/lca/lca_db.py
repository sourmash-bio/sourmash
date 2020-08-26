"LCA database class and utilities."

from __future__ import print_function, division
import json
import gzip
from collections import OrderedDict, defaultdict, Counter
import functools

import sourmash
from sourmash._minhash import get_max_hash_for_scaled, to_bytes
from sourmash.logging import notify, error, debug
from sourmash.index import Index
from sourmash.exceptions import SourmashError

from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str

def cached_property(fun):
    """A memoize decorator for class properties."""
    @functools.wraps(fun)
    def get(self):
        try:
            return self._cache[fun]
        except AttributeError:
            self._cache = {}
        except KeyError:
            pass
        ret = self._cache[fun] = fun(self)
        return ret
    return property(get)


class LCA_Database(RustObject, Index):
    """
    An in-memory database that indexes signatures by hash, and provides
    optional taxonomic lineage classification.

    Follows the `Index` API for `insert`, `search`, `gather`, and `signatures`.

    Identifiers `ident` must be unique, and are taken by default as the
    entire signature name upon insertion. This can be overridden with
    the `ident` keyword argument in `insert`.

    Integer `idx` indices can be used as keys in dictionary attributes:
    * `idx_to_lid`, to get an (optional) lineage index.
    * `idx_to_ident`, to retrieve the unique string identifier for that `idx`.

    Integer `lid` indices can be used as keys in dictionary attributes:
    * `lid_to_idx`, to get a set of `idx` with that lineage.
    * `lid_to_lineage`, to get a lineage for that `lid`.

    `lineage_to_lid` is a dictionary with tuples of LineagePair as keys,
    `lid` as values.

    `ident_to_name` is a dictionary from unique str identifer to a name.

    `ident_to_idx` is a dictionary from unique str identifer to integer `idx`.

    `hashval_to_idx` is a dictionary from individual hash values to sets of
    `idx`.
    """
    def __init__(self, ksize, scaled, moltype="DNA", filename=""):
        self._objptr = lib.lcadb_new_with_params(int(ksize), int(scaled), to_bytes(filename), to_bytes(moltype))

    @cached_property
    def ksize(self):
        return self._methodcall(lib.lcadb_ksize)

    @cached_property
    def scaled(self):
        return self._methodcall(lib.lcadb_scaled)

    @cached_property
    def moltype(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_moltype, size)

        size = size[0]
        buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
        moltype = ffi.string(buf, size).decode('utf-8')

        return moltype

    @cached_property
    def filename(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_filename, size)

        size = size[0]
        filename = None
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            filename = ffi.string(buf, size).decode('utf-8')

        return filename

    # --------------------------------------
    # new helper functions... keep?

    # TODO: hopefully finish this function: to get rid of dependency on last exposed internal properties.
    def _make_assignments_helper(self, min_num):
        return self._methodcall(lib.make_assignments_helper, int(min_num))

    def _invalidate_cache(self):
        if hasattr(self, '_cache'):
            del self._cache

    def _hashval_to_idx_len(self):
        return self._methodcall(lib.lcadb_hashval_to_idx_len)

    def _get_match_size(self, best_idx):
        return self._methodcall(lib.lcadb_get_match_size, int(best_idx))

    def _best_name(self, best_idx):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_best_name, int(best_idx), size)
        size = size[0]
        buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)

        return ffi.string(buf, size).decode('utf-8')

    def _get_lineage_from_idx(self, idx):
        from .lca_utils import taxlist, LineagePair

        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_get_lineage_from_idx, int(idx), size)

        size = size[0]
        buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
        lineage_vec = json.loads(ffi.string(buf, size).decode('utf-8'))

        result = []
        for pair in lineage_vec:
            result.append(LineagePair(rank=pair[0], name=pair[1]))

        return tuple(result)

    def _get_idx_from_hashval(self, hashval):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_get_idx_from_hashval, int(hashval), size)

        size = size[0]

        result = []
        for i in range(size):
            result.append(rawbuff[i])

        return result

    def _get_hashvals(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_get_hashvals, size)
        size = size[0]

        result = []
        for i in range(size):
            result.append(rawbuff[i])

        return result

    # --------------------------------------


    def insert(self, sig, ident=None, lineage=None):
        from .lca_utils import taxlist, LineagePair
        """Add a new signature into the LCA database.

        Takes optional arguments 'ident' and 'lineage'.

        'ident' must be a unique string identifer across this database;
        if not specified, the signature name (sig.name()) is used.

        'lineage', if specified, must contain a tuple of LineagePair objects.
        """
        minhash = sig.minhash

        if minhash.ksize != self.ksize:
            raise ValueError("cannot insert signature with ksize {} into DB (ksize {})".format(minhash.ksize, self.ksize))

        if minhash.moltype != self.moltype:
            raise ValueError("cannot insert signature with moltype {} into DB (moltype {})".format(minhash.moltype, self.moltype))

        if ident is None:
            ident = sig.name()

        self._invalidate_cache()

        # trying to put lineage into a json
        if lineage:
            try:
                if type(lineage[0]) is not LineagePair:
                    raise ValueError('lineage cannot be used as a key?!')
            except:
                raise ValueError('lineage cannot be used as a key?!')

            lineage_dict = dict(lineage)
            lineage_json = json.dumps(lineage_dict)
        else:
            lineage_json = ""

        try:
            return self._methodcall(lib.lcadb_insert, sig._objptr, to_bytes(ident), to_bytes(lineage_json))
        except SourmashError as e:
            raise ValueError(str(e))


    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def signatures(self):
        "Return all of the signatures in this LCA database."
        from sourmash import SourmashSignature

        size = ffi.new("uintptr_t *")

        sigs_ptr = self._methodcall(lib.lcadb_signatures, size)

        size = size[0]
        sigs = []
        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            sigs.append(sig)

        for sig in sigs:
            yield sig

    def select(self, ksize=None, moltype=None):
        "Selector interface - make sure this database matches requirements."
        ok = True
        if ksize is not None and self.ksize != ksize:
            ok = False
        if moltype is not None and moltype != self.moltype:
            ok = False

        if ok:
            return self

        raise ValueError("cannot select LCA on ksize {} / moltype {}".format(ksize, moltype))

    @classmethod
    def load(cls, db_name):
        "Load LCA_Database from a JSON file."
        from .lca_utils import taxlist, LineagePair

        dbs_ptr = rustcall(lib.lcadb_load_db, to_bytes(db_name))

        db = LCA_Database._from_objptr(dbs_ptr)

        return db

    def save(self, db_name):
        "Save LCA_Database to a JSON file."
        self._methodcall(lib.lcadb_save, to_bytes(db_name))

    def search(self, query, *args, **kwargs):
        from sourmash import SourmashSignature
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments:
          * do_containment: default False. If True, use Jaccard containment.
          * best_only: default False. If True, allow optimizations that
            may. May discard matches better than threshold, but first match
            is guaranteed to be best.
          * ignore_abundance: default False. If True, and query signature
            and database support k-mer abundances, ignore those abundances.

        Note, the "best only" hint is ignored by LCA_Database
        """
        if not query.minhash:
            return []

        # check arguments
        if 'threshold' not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs['threshold']
        do_containment = kwargs.get('do_containment', False)
        ignore_abundance = kwargs.get('ignore_abundance', False)

        size = ffi.new("uintptr_t *")
        search_results = self._methodcall(lib.lcadb_search, query._objptr, threshold, do_containment, ignore_abundance, size)
        size = size[0]

        results = []
        for i in range(size):
            # set filename
            name_size = search_results[i].name_size
            filename = None
            if name_size:
                rawbuff = search_results[i].filename
                buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, name_size), name_size)
                filename = ffi.string(buf, name_size).decode('utf-8')

            sig = SourmashSignature._from_objptr(search_results[i].sig)

            tup = (search_results[i].score, sig, filename)
            results.append(tup)

        return results
        

    def gather(self, query, *args, **kwargs):
        from sourmash import SourmashSignature
        "Return the match with the best Jaccard containment in the database."
        if not query.minhash:
            return []

        results = []
        threshold_bp = kwargs.get('threshold_bp', 0.0)

        size = ffi.new("uintptr_t *")
        gather_results = self._methodcall(lib.lcadb_gather, query._objptr, threshold_bp, size)
        size = size[0]

        for i in range(size):
            # set filename
            name_size = gather_results[i].name_size
            filename = None
            if name_size > 0:
                rawbuff = gather_results[i].filename
                buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, name_size), name_size)
                filename = ffi.string(buf, name_size).decode('utf-8')

            sig = SourmashSignature._from_objptr(gather_results[i].sig)

            tup = (gather_results[i].score, sig, filename)
            results.append(tup)

        results.sort(key=lambda x: -x[0])
        return results

    def find(self, search_fn, *args, **kwargs):
        """Not implemented; 'find' cannot be implemented efficiently on
        an LCA database."""
        raise NotImplementedError

    def downsample_scaled(self, scaled):
        """
        Downsample to the provided scaled value, i.e. eliminate all hashes
        that don't fall in the required range.

        This applies to this database in place.
        """
        if scaled == self.scaled:
            return
        elif scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        self._invalidate_cache()

        self._methodcall(lib.lcadb_downsample_scaled, scaled)

    def get_lineage_assignments(self, hashval):
        from .lca_utils import LineagePair
        """
        Get a list of lineages for this hashval.
        """
        size = ffi.new("uintptr_t *")
        rawbuf = self._methodcall(lib.lcadb_get_lineage_assignments, int(hashval), size)
        size = size[0]

        buf = ffi.gc(rawbuf, lambda o: lib.nodegraph_buffer_free(o, size), size)
        result = ffi.string(buf, size)
        result = json.loads(result)

        lineagelist = []
        for v in result:
            lineage = []
            for vv in v:
                lineage.append(LineagePair(rank=vv[0], name=vv[1]))
            
            lineagelist.append(tuple(lineage))

        return lineagelist


def load_single_database(filename, verbose=False):
    "Load a single LCA database; return (db, ksize, scaled)"
    dblist, ksize, scaled = load_databases([filename], verbose=verbose)
    return dblist[0], ksize, scaled


def load_databases(filenames, scaled=None, verbose=True):
    "Load multiple LCA databases; return (dblist, ksize, scaled)"
    ksize_vals = set()
    scaled_vals = set()
    moltype_vals = set()
    dblist = []

    # load all the databases
    for db_name in filenames:
        if verbose:
            notify(u'\r\033[K', end=u'')
            notify('... loading database {}'.format(db_name), end='\r')

        lca_db = LCA_Database.load(db_name)

        ksize_vals.add(lca_db.ksize)
        if len(ksize_vals) > 1:
            raise Exception('multiple ksizes, quitting')

        if scaled and scaled > lca_db.scaled:
            lca_db.downsample_scaled(scaled)
        scaled_vals.add(lca_db.scaled)

        moltype_vals.add(lca_db.moltype)
        if len(moltype_vals) > 1:
            raise Exception('multiple moltypes, quitting')

        dblist.append(lca_db)

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()
    moltype = moltype_vals.pop()

    if verbose:
        notify(u'\r\033[K', end=u'')
        notify('loaded {} LCA databases. ksize={}, scaled={} moltype={}',
               len(dblist), ksize, scaled, moltype)

    return dblist, ksize, scaled
