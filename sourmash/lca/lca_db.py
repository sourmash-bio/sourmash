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

        # self.ksize = ksize
        # self.scaled = scaled
        # self.filename = ""
        # self.moltype = moltype

        # self._next_index = 0
        # self._next_lid = 0
        # self.ident_to_name = {}
        # self.ident_to_idx = {}
        # self.idx_to_lid = {}
        # self.lineage_to_lid = {}
        # self.lid_to_lineage = {}
        # self.hashval_to_idx = defaultdict(set)

        self._objptr = lib.lcadb_new_with_params(int(ksize), int(scaled), to_bytes(filename), to_bytes(moltype))

    @property
    def ksize(self):
        return self._methodcall(lib.lcadb_ksize)

    @property
    def scaled(self):
        return self._methodcall(lib.lcadb_scaled)

    @property
    def moltype(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_moltype, size)

        size = size[0]
        buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
        moltype = ffi.string(buf, size).decode('utf-8')

        return moltype

    @property
    def filename(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_filename, size)

        size = size[0]
        filename = None
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            filename = ffi.string(buf, size).decode('utf-8')

        return filename

    @property
    def ident_to_name(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_ident_to_name, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            ident_to_name = json.loads(ffi.string(buf, size).decode('utf-8'))

        return ident_to_name

    @property
    def ident_to_idx(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_ident_to_idx, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            ident_to_idx = json.loads(ffi.string(buf, size).decode('utf-8'))

        result = {}
        for ident, idx in ident_to_idx.items():
            result[ident] = int(idx)

        return result

    @property
    def idx_to_lid(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_idx_to_lid, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            idx_to_lid = json.loads(ffi.string(buf, size).decode('utf-8'))

        result = {}
        for idx, lid in idx_to_lid.items():
            result[int(idx)] = int(lid)
        
        return result

    @property
    def lineage_to_lid(self):
        from .lca_utils import taxlist, LineagePair

        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_lid_to_lineage, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            lid_to_lineage = json.loads(ffi.string(buf, size).decode('utf-8'))

        # turn lineage into a list of tuples
        result = {}
        for (lid, lineage) in lid_to_lineage.items():
            pairslist = []
            for v in lineage:
                pairslist.append(LineagePair(rank=v[0], name=v[1]))

            result[tuple(pairslist)] = int(lid)

        return result

    @property
    def lid_to_lineage(self):
        from .lca_utils import taxlist, LineagePair

        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_lid_to_lineage, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            temp = ffi.string(buf, size).decode('utf-8')
            lid_to_lineage = json.loads(temp)

        # turn lineage into a list of tuples
        result = {}
        for (lid, lineage) in lid_to_lineage.items():
            pairslist = []
            for v in lineage:
                pairslist.append(LineagePair(rank=v[0], name=v[1]))

            result[int(lid)] = tuple(pairslist)

        return result

    @property
    def hashval_to_idx(self):
        size = ffi.new('uintptr_t *')
        rawbuff = self._methodcall(lib.lcadb_hashval_to_idx, size)

        size = size[0]
        if size:
            buf = ffi.gc(rawbuff, lambda o: lib.nodegraph_buffer_free(o, size), size)
            hashval_to_idx = json.loads(ffi.string(buf, size).decode('utf-8'))

        # convert json's u64 string to int
        result = {}
        for hashval, idxlist in hashval_to_idx.items():
            result[int(hashval)] = idxlist

        return result

    def _invalidate_cache(self):
        if hasattr(self, '_cache'):
            del self._cache

    def _hashval_to_idx_len(self):
        return self._methodcall(lib.lcadb_hashval_to_idx_len)

    def _get_ident_index(self, ident, fail_on_duplicate=False):
        "Get (create if nec) a unique int id, idx, for each identifier."
        idx = self.ident_to_idx.get(ident)
        if fail_on_duplicate:
            assert idx is None     # should be no duplicate identities

        if idx is None:
            idx = self._next_index
            self._next_index += 1

            self.ident_to_idx[ident] = idx

        return idx

    def _get_lineage_id(self, lineage):
        "Get (create if nec) a unique lineage ID for each LineagePair tuples."
        # does one exist already?
        lid = self.lineage_to_lid.get(lineage)

        # nope - create one. Increment next_lid.
        if lid is None:
            lid = self._next_lid
            self._next_lid += 1

            # build mappings
            self.lineage_to_lid[lineage] = lid
            self.lid_to_lineage[lid] = lineage

        return lid

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

        # downsample to specified scaled; this has the side effect of
        # making sure they're all at the same scaled value!
        try:
            minhash = minhash.downsample_scaled(self.scaled)
        except ValueError:
            raise ValueError("cannot downsample signature; is it a scaled signature?")

        if ident is None:
            ident = sig.name()

        # if ident in self.ident_to_name:
        #     raise ValueError("signature {} is already in this LCA db.".format(ident))

        # before adding, invalide any caching from @cached_property
        self._invalidate_cache()

        # # store full name
        # self.ident_to_name[ident] = sig.name()

        # # identifier -> integer index (idx)
        # idx = self._get_ident_index(ident, fail_on_duplicate=True)
        # if lineage:
        #     try:
        #         lineage = tuple(lineage)

        #         # (LineagePairs*) -> integer lineage ids (lids)
        #         lid = self._get_lineage_id(lineage)

        #         # map idx to lid as well.
        #         self.idx_to_lid[idx] = lid
        #     except TypeError:
        #         raise ValueError('lineage cannot be used as a key?!')

        # for hashval in minhash.get_mins():
        #     self.hashval_to_idx[hashval].add(idx)


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

        return self._methodcall(lib.lcadb_insert, sig._objptr, to_bytes(ident), to_bytes(lineage_json))


    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def signatures(self):
        "Return all of the signatures in this LCA database."
        from sourmash import SourmashSignature
        # print("\n", self._signatures, "\n")
        # for v in self._signatures.values():
        #     print(v)
        #     yield v
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

        # xopen = open
        # if db_name.endswith('.gz'):
        #     xopen = gzip.open

        # with xopen(db_name, 'rt') as fp:
        #     load_d = {}
        #     try:
        #         load_d = json.load(fp)
        #     except json.decoder.JSONDecodeError:
        #         pass

        #     if not load_d:
        #         raise ValueError("cannot parse database file '{}' as JSON; invalid format.")

        #     version = None
        #     db_type = None
        #     try:
        #         version = load_d.get('version')
        #         db_type = load_d.get('type')
        #     except AttributeError:
        #         pass

        #     if db_type != 'sourmash_lca':
        #         raise ValueError("database file '{}' is not an LCA db.".format(db_name))

        #     version = float(version)
        #     if version < 2.0 or 'lid_to_lineage' not in load_d:
        #         raise ValueError("Error! This is an old-style LCA DB. You'll need to rebuild or download a newer one.")

        #     ksize = int(load_d['ksize'])
        #     scaled = int(load_d['scaled'])
        #     moltype = load_d.get('moltype', 'DNA')

        #     db = cls(ksize, scaled, moltype)

        #     # convert lineage_dict to proper lineages (tuples of LineagePairs)
        #     lid_to_lineage_2 = load_d['lid_to_lineage']
        #     lid_to_lineage = {}
        #     lineage_to_lid = {}
        #     for k, v in lid_to_lineage_2.items():
        #         v = dict(v)
        #         vv = []
        #         for rank in taxlist():
        #             name = v.get(rank, '')
        #             vv.append(LineagePair(rank, name))

        #         vv = tuple(vv)
        #         lid_to_lineage[int(k)] = vv
        #         lineage_to_lid[vv] = int(k)
        #     db.lid_to_lineage = lid_to_lineage
        #     db.lineage_to_lid = lineage_to_lid

        #     # convert hashval -> lineage index keys to integers (looks like
        #     # JSON doesn't have a 64 bit type so stores them as strings)
        #     hashval_to_idx_2 = load_d['hashval_to_idx']
        #     hashval_to_idx = {}

        #     for k, v in hashval_to_idx_2.items():
        #         hashval_to_idx[int(k)] = v
        #     db.hashval_to_idx = hashval_to_idx

        #     db.ident_to_name = load_d['ident_to_name']
        #     db.ident_to_idx = load_d['ident_to_idx']

        #     db.idx_to_lid = {}
        #     for k, v in load_d['idx_to_lid'].items():
        #         db.idx_to_lid[int(k)] = v

        # print("\n\n\n", type(db_ptr), "\n\n\n")

        dbs_ptr = rustcall(lib.lcadb_load_db, to_bytes(db_name))

        db = LCA_Database._from_objptr(dbs_ptr)

        return db

    def save(self, db_name):
        "Save LCA_Database to a JSON file."
        
        # xopen = open
        # if db_name.endswith('.gz'):
        #     xopen = gzip.open

        # with xopen(db_name, 'wt') as fp:
        #     # use an OrderedDict to preserve output order
        #     save_d = OrderedDict()
        #     save_d['version'] = '2.1'
        #     save_d['type'] = 'sourmash_lca'
        #     save_d['license'] = 'CC0'
        #     save_d['ksize'] = self.ksize
        #     save_d['scaled'] = self.scaled
        #     save_d['moltype'] = self.moltype

        #     # convert lineage internals from tuples to dictionaries
        #     d = OrderedDict()
        #     for k, v in self.lid_to_lineage.items():
        #         d[k] = dict([ (vv.rank, vv.name) for vv in v ])
        #     save_d['lid_to_lineage'] = d

        #     # convert values from sets to lists, so that JSON knows how to save
        #     save_d['hashval_to_idx'] = \
        #        dict((k, list(v)) for (k, v) in self.hashval_to_idx.items())

        #     save_d['ident_to_name'] = self.ident_to_name
        #     save_d['ident_to_idx'] = self.ident_to_idx
        #     save_d['idx_to_lid'] = self.idx_to_lid
        #     save_d['lid_to_lineage'] = self.lid_to_lineage
            
        #     # print("\n\nPYTHON:\n", save_d, "\n\n")
        #     json.dump(save_d, fp)
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
        # mh = query.minhash
        # if ignore_abundance:
        #     mh.track_abundance = False

        # # find all the matches, then sort & return.
        # results = []
        # for x in self._find_signatures(mh, threshold, do_containment):
        #     (score, match, filename) = x
        #     results.append((score, match, filename))

        # results.sort(key=lambda x: -x[0])
        # return results
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
        # threshold = threshold_bp / (len(query.minhash) * self.scaled)

        # # grab first match, if any, and return that; since _find_signatures
        # # is a generator, this will truncate further searches.
        # for x in self._find_signatures(query.minhash, threshold,
        #                                containment=True, ignore_scaled=True):
        #     (score, match, filename) = x
        #     if score:
        #         results.append((score, match, filename))
        #     break

        # return results
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

        # max_hash = get_max_hash_for_scaled(scaled)

        # # filter out all hashes over max_hash in value.
        # new_hashvals = {}
        # for k, v in self.hashval_to_idx.items():
        #     if k < max_hash:
        #         new_hashvals[k] = v
        # self.hashval_to_idx = new_hashvals
        # self.scaled = scaled

    def get_lineage_assignments(self, hashval):
        from .lca_utils import LineagePair
        """
        Get a list of lineages for this hashval.
        """
        # x = []

        # idx_list = self.hashval_to_idx.get(hashval, [])
        # for idx in idx_list:
        #     lid = self.idx_to_lid.get(idx, None)
        #     if lid is not None:
        #         lineage = self.lid_to_lineage[lid]
        #         x.append(lineage)

        # return x
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

    @cached_property
    def _signatures(self):
        "Create a _signatures member dictionary that contains {idx: sigobj}."
        from sourmash import MinHash, SourmashSignature

        is_protein = False
        is_hp = False
        is_dayhoff = False
        if self.moltype == 'protein':
            is_protein = True
        elif self.moltype == 'hp':
            is_hp = True
        elif self.moltype == 'dayhoff':
            is_dayhoff = True
        minhash = MinHash(n=0, ksize=self.ksize, scaled=self.scaled,
                          is_protein=is_protein, hp=is_hp, dayhoff=is_dayhoff)

        debug('creating signatures for LCA DB...')
        mhd = defaultdict(minhash.copy_and_clear)
        temp_vals = defaultdict(list)

        # invert the hashval_to_idx dictionary
        for (hashval, idlist) in self.hashval_to_idx.items():
            for idx in idlist:
                temp_hashes = temp_vals[idx]
                temp_hashes.append(hashval)

                # 50 is an arbitrary number. If you really want
                # to micro-optimize, list is resized and grow in this pattern:
                # 0, 4, 8, 16, 25, 35, 46, 58, 72, 88, ...
                # (from https://github.com/python/cpython/blob/b2b4a51f7463a0392456f7772f33223e57fa4ccc/Objects/listobject.c#L57)
                if len(temp_hashes) > 50:
                    mhd[idx].add_many(temp_hashes)

                    # Sigh, python 2... when it goes away,
                    # we can do `temp_hashes.clear()` instead.
                    del temp_vals[idx]

        # We loop temp_vals again to add any remainder hashes
        # (each list of hashes is smaller than 50 items)
        for sig, vals in temp_vals.items():
            mhd[sig].add_many(vals)

        sigd = {}
        for idx, mh in mhd.items():
            ident = self.idx_to_ident[idx]
            name = self.ident_to_name[ident]
            sigd[idx] = SourmashSignature(mh, name=name)

        debug('=> {} signatures!', len(sigd))
        return sigd

    def _find_signatures(self, minhash, threshold, containment=False,
                       ignore_scaled=False):
        """
        Do a Jaccard similarity or containment search, yield results.

        This is essentially a fast implementation of find that collects all
        the signatures with overlapping hash values. Note that similarity
        searches (containment=False) will not be returned in sorted order.
        """
        # make sure we're looking at the same scaled value as database
        if self.scaled > minhash.scaled:
            minhash = minhash.downsample_scaled(self.scaled)
        elif self.scaled < minhash.scaled and not ignore_scaled:
            # note that containment can be calculated w/o matching scaled.
            raise ValueError("lca db scaled is {} vs query {}; must downsample".format(self.scaled, minhash.scaled))

        query_mins = set(minhash.get_mins())

        # collect matching hashes for the query:
        c = Counter()
        for hashval in query_mins:
            idx_list = self.hashval_to_idx.get(hashval, [])
            for idx in idx_list:
                c[idx] += 1

        debug('number of matching signatures for hashes: {}', len(c))

        # for each match, in order of largest overlap,
        for idx, count in c.most_common():
            # pull in the hashes. This reconstructs & caches all input
            # minhashes, which is kinda memory intensive...!
            # NOTE: one future low-mem optimization could be to support doing
            # this piecemeal by iterating across all the hashes, instead.
            match_sig = self._signatures[idx]
            match_mh = match_sig.minhash
            match_size = len(match_mh)

            # calculate the containment or similarity
            if containment:
                score = count / len(query_mins)
            else:
                # query_mins is size of query signature
                # match_size is size of match signature
                # count is overlap
                score = count / (len(query_mins) + match_size - count)

            # ...and return.
            if score >= threshold:
                yield score, match_sig, self.filename

    @cached_property
    def lid_to_idx(self):
        d = defaultdict(set)
        for idx, lid in self.idx_to_lid.items():
            d[lid].add(idx)
        return d

    @cached_property
    def idx_to_ident(self):
        d = defaultdict(set)
        for ident, idx in self.ident_to_idx.items():
            assert idx not in d
            d[idx] = ident
        return d


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
