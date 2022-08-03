"LCA database class and utilities."
import os
import json
import gzip
from collections import OrderedDict, defaultdict, Counter
import functools

import sourmash
from sourmash.minhash import _get_max_hash_for_scaled
from sourmash.logging import notify, error, debug
from sourmash.index import Index, IndexSearchResult
from sourmash.picklist import passes_all_picklists


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


class LCA_Database(Index):
    """
    An in-memory database that indexes signatures by hash, and provides
    optional taxonomic lineage classification.

    Follows the `Index` API for `insert`, `search`, `gather`, and `signatures`.

    Identifiers `ident` must be unique, and are taken by default as the
    entire signature name upon insertion. This can be overridden with
    the `ident` keyword argument in `insert`.

    Integer `idx` indices can be used as keys in dictionary attributes:
    * `_idx_to_lid`, to get an (optional) lineage index.
    * `_idx_to_ident`, to retrieve the unique string identifier for that `idx`.

    Integer `lid` indices can be used as keys in dictionary attributes:
    * `_lid_to_idx`, to get a set of `idx` with that lineage.
    * `_lid_to_lineage`, to get a lineage for that `lid`.

    `_lineage_to_lid` is a dictionary with tuples of LineagePair as keys,
    `lid` as values.

    `_ident_to_name` is a dictionary from unique str identifer to a name.

    `_ident_to_idx` is a dictionary from unique str identifer to integer `idx`.

    `_hashval_to_idx` is a dictionary from individual hash values to sets of
    `idx`.
    """
    is_database = True

    # we set manifest to None to avoid implication of fast on-disk access to
    # sketches. This may be revisited later.
    manifest = None

    def __init__(self, ksize, scaled, moltype='DNA'):
        self.ksize = int(ksize)
        self.scaled = int(scaled)
        self.filename = None
        self.moltype = moltype

        self._next_index = 0
        self._next_lid = 0
        self._ident_to_name = {}
        self._ident_to_idx = {}
        self._idx_to_lid = {}
        self._lineage_to_lid = {}
        self._lid_to_lineage = {}
        self._hashval_to_idx = defaultdict(set)
        self.picklists = []

    @property
    def location(self):
        """Return source filename.

        Part of the Index protocol.
        """
        return self.filename

    def __len__(self):
        """Return number of sketches.

        Part of the Index protocol.
        """
        return self._next_index

    def _invalidate_cache(self):
        """Force rebuild of signatures after an 'insert'.

        Internal method.
        """
        if hasattr(self, '_cache'):
            del self._cache

    def _get_ident_index(self, ident, fail_on_duplicate=False):
        """Get (create if necessary) a unique int idx, for each identifier.

        Internal method.
        """
        idx = self._ident_to_idx.get(ident)
        if fail_on_duplicate:
            assert idx is None     # should be no duplicate identities

        if idx is None:
            idx = self._next_index
            self._next_index += 1

            self._ident_to_idx[ident] = idx

        return idx

    def _get_lineage_id(self, lineage):
        """Get (create if necessary) a unique lineage ID for each
        LineagePair tuples."

        Internal method of this class.
        """
        # does one exist already?
        lid = self._lineage_to_lid.get(lineage)

        # nope - create one. Increment next_lid.
        if lid is None:
            lid = self._next_lid
            self._next_lid += 1

            # build mappings
            self._lineage_to_lid[lineage] = lid
            self._lid_to_lineage[lid] = lineage

        return lid

    def insert(self, sig, ident=None, lineage=None):
        """Add a new signature into the LCA database.

        Takes optional arguments 'ident' and 'lineage'.

        'ident' must be a unique string identifer across this database;
        if not specified, the signature name (sig.name) is used.

        'lineage', if specified, must contain a tuple of LineagePair objects.

        Method unique to this class.
        """
        minhash = sig.minhash

        if minhash.ksize != self.ksize:
            raise ValueError("cannot insert signature with ksize {} into DB (ksize {})".format(minhash.ksize, self.ksize))

        if minhash.moltype != self.moltype:
            raise ValueError("cannot insert signature with moltype {} into DB (moltype {})".format(minhash.moltype, self.moltype))

        # downsample to specified scaled; this has the side effect of
        # making sure they're all at the same scaled value!
        try:
            minhash = minhash.downsample(scaled=self.scaled)
        except ValueError:
            raise ValueError("cannot downsample signature; is it a scaled signature?")

        if not ident:
            ident = str(sig)

        if ident in self._ident_to_name:
            raise ValueError("signature '{}' is already in this LCA db.".format(ident))

        # before adding, invalide any caching from @cached_property
        self._invalidate_cache()

        # store full name
        self._ident_to_name[ident] = sig.name

        # identifier -> integer index (idx)
        idx = self._get_ident_index(ident, fail_on_duplicate=True)
        if lineage:
            try:
                lineage = tuple(lineage)

                # (LineagePairs*) -> integer lineage ids (lids)
                lid = self._get_lineage_id(lineage)

                # map idx to lid as well.
                self._idx_to_lid[idx] = lid
            except TypeError:
                raise ValueError('lineage cannot be used as a key?!')

        for hashval in minhash.hashes:
            self._hashval_to_idx[hashval].add(idx)

        return len(minhash)

    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def signatures(self):
        """Return all of the signatures in this LCA database.

        Part of the Index protocol.
        """
        from sourmash import SourmashSignature

        if self.picklists:
            pl = self.picklists
            for v in self._signatures.values():
                if passes_all_picklists(v, pl):
                    yield v
        else:
            for v in self._signatures.values():
                yield v

    def _signatures_with_internal(self):
        """Return all of the signatures in this LCA database.

        Part of the Index protocol; used for buulding manifests.
        """

        for idx, ss in self._signatures.items():
            yield ss, idx

    def select(self, ksize=None, moltype=None, num=0, scaled=0, abund=None,
               containment=False, picklist=None):
        """Select a subset of signatures to search.

        As with SBTs, queries with higher scaled values than the database
        can still be used for containment search, but not for similarity
        search. See SBT.select(...) for details, and _find_signatures for
        implementation.

        Will always raise ValueError if a requirement cannot be met.
        """
        if num:
            raise ValueError("cannot use 'num' MinHashes to search LCA database")

        if scaled > self.scaled and not containment:
            raise ValueError(f"cannot use scaled={scaled} on this database (scaled={self.scaled})")

        if ksize is not None and self.ksize != ksize:
            raise ValueError(f"ksize on this database is {self.ksize}; this is different from requested ksize of {ksize}")
        if moltype is not None and moltype != self.moltype:
            raise ValueError(f"moltype on this database is {self.moltype}; this is different from requested moltype of {moltype}")

        if abund:
            raise ValueError("LCA databases do not support sketches with abund=True")

        if picklist is not None:
            self.picklists.append(picklist)
            if len(self.picklists) > 1:
                raise ValueError("we do not (yet) support multiple picklists for LCA databases")

        return self

    @classmethod
    def load(cls, db_name):
        """Load LCA_Database from a JSON file.

        Method specific to this class.
        """
        from .lca_utils import taxlist, LineagePair

        if not os.path.isfile(db_name):
            raise ValueError(f"'{db_name}' is not a file and cannot be loaded as an LCA database")

        try:
            from sourmash.index.sqlite_index import LCA_SqliteDatabase
            return LCA_SqliteDatabase.load(db_name)
        except ValueError:
            pass

        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
            try:
                first_ch = fp.read(1)
            except ValueError:
                first_ch = 'X'
            if not first_ch or first_ch[0] != '{':
                raise ValueError(f"'{db_name}' is not an LCA database file.")

            fp.seek(0)

            load_d = {}
            try:
                load_d = json.load(fp)
            except json.decoder.JSONDecodeError:
                pass

            if not load_d:
                raise ValueError("cannot parse database file '{}' as JSON; invalid format.")

            version = None
            db_type = None
            try:
                version = load_d.get('version')
                db_type = load_d.get('type')
            except AttributeError:
                pass

            if db_type != 'sourmash_lca':
                raise ValueError("database file '{}' is not an LCA db.".format(db_name))

            version = float(version)
            if version < 2.0 or 'lid_to_lineage' not in load_d:
                raise ValueError("Error! This is an old-style LCA DB. You'll need to rebuild or download a newer one.")

            ksize = int(load_d['ksize'])
            scaled = int(load_d['scaled'])
            moltype = load_d.get('moltype', 'DNA')
            if moltype != 'DNA':
                assert ksize % 3 == 0
                ksize = int(ksize / 3)

            db = cls(ksize, scaled, moltype)

            # convert lineage_dict to proper lineages (tuples of LineagePairs)
            lid_to_lineage_2 = load_d['lid_to_lineage']
            lid_to_lineage = {}
            lineage_to_lid = {}
            for k, v in lid_to_lineage_2.items():
                v = dict(v)
                vv = []
                for rank in taxlist():
                    name = v.get(rank, '')
                    vv.append(LineagePair(rank, name))

                vv = tuple(vv)
                lid_to_lineage[int(k)] = vv
                lineage_to_lid[vv] = int(k)
            db._lid_to_lineage = lid_to_lineage
            db._lineage_to_lid = lineage_to_lid

            # convert hashval -> lineage index keys to integers (looks like
            # JSON doesn't have a 64 bit type so stores them as strings)
            hashval_to_idx_2 = load_d['hashval_to_idx']
            hashval_to_idx = {}

            for k, v in hashval_to_idx_2.items():
                hashval_to_idx[int(k)] = v
            db._hashval_to_idx = hashval_to_idx

            db._ident_to_name = load_d['ident_to_name']
            db._ident_to_idx = load_d['ident_to_idx']

            db._idx_to_lid = {}
            for k, v in load_d['idx_to_lid'].items():
                db._idx_to_lid[int(k)] = v

        if db._ident_to_idx:
            db._next_index = max(db._ident_to_idx.values()) + 1
        else:
            db._next_index = 0
        if db._idx_to_lid:
            db._next_lid = max(db._idx_to_lid.values()) + 1
        else:
            db._next_lid = 0

        db.filename = db_name

        return db

    def save(self, db_name, *, format='json'):
        if format == 'sql':
            self.save_to_sql(db_name)
        else:
            assert format == 'json'
            self.save_to_json(db_name)

    def save_to_json(self, db_name):
        """Save LCA_Database to a JSON file.

        Method specific to this class.
        """
        if os.path.exists(db_name):
            raise ValueError(f"LCA database {db_name} already exists; not overwriting or appending")

        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '2.1'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'

            if self.moltype != 'DNA':
                ksize = self.ksize*3
            else:
                ksize = self.ksize
            save_d['ksize'] = ksize
            save_d['scaled'] = self.scaled
            save_d['moltype'] = self.moltype

            # convert lineage internals from tuples to dictionaries
            d = OrderedDict()
            for k, v in self._lid_to_lineage.items():
                d[k] = dict([ (vv.rank, vv.name) for vv in v ])
            save_d['lid_to_lineage'] = d

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_to_idx'] = \
               dict((k, list(v)) for (k, v) in self._hashval_to_idx.items())

            save_d['ident_to_name'] = self._ident_to_name
            save_d['ident_to_idx'] = self._ident_to_idx
            save_d['idx_to_lid'] = self._idx_to_lid
            save_d['lid_to_lineage'] = self._lid_to_lineage
            
            json.dump(save_d, fp)

    def save_to_sql(self, dbname):
        "Save this LCA_Database into an LCA_SqliteDatabase"
        from sourmash.index.sqlite_index import LCA_SqliteDatabase
        from sourmash.tax.tax_utils import LineageDB

        if os.path.exists(dbname):
            raise ValueError(f"LCA database {dbname} already exists; not overwriting or appending")

        # create a new in-memory lineage db...
        assignments = {}
        available_ranks = set() # track ranks, too
        for ident, idx in self._ident_to_idx.items():
            lid = self._idx_to_lid.get(idx)
            if lid is not None:
                lineage = self._lid_to_lineage[lid]
                assignments[ident] = lineage
                for pair in lineage:
                    available_ranks.add(pair.rank)

        ldb = LineageDB(assignments, available_ranks)

        # ...and pass over to create, using 'self' as index.
        LCA_SqliteDatabase.create(dbname, self, ldb)

    def downsample_scaled(self, scaled):
        """
        Downsample to the provided scaled value, i.e. eliminate all hashes
        that don't fall in the required range.

        This applies to this database in place.

        Method specific to LCA databases.
        """
        if scaled == self.scaled:
            return
        elif scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        self._invalidate_cache()

        max_hash = _get_max_hash_for_scaled(scaled)

        # filter out all hashes over max_hash in value.
        new_hashvals = defaultdict(set)
        for k, v in self._hashval_to_idx.items():
            if k < max_hash:
                new_hashvals[k] = v
        self._hashval_to_idx = new_hashvals
        self.scaled = scaled

    @property
    def hashvals(self):
        "Return all hashvals stored in this database."
        return self._hashval_to_idx.keys()

    def get_lineage_assignments(self, hashval, min_num=None):
        """Get a list of lineages for this hashval.

        Method specific to LCA Databases.
        """
        x = []

        idx_list = self._hashval_to_idx.get(hashval, [])

        if min_num and len(idx_list) < min_num:
            return []

        for idx in idx_list:
            lid = self._idx_to_lid.get(idx, None)
            if lid is not None:
                lineage = self._lid_to_lineage[lid]
                x.append(lineage)

        return x

    def get_identifiers_for_hashval(self, hashval):
        """
        Get a list of identifiers for signatures containing this hashval
        """
        idx_list = self._hashval_to_idx.get(hashval, [])

        for idx in idx_list:
            yield self._idx_to_ident[idx]

    @cached_property
    def _signatures(self):
        """Create a _signatures member dictionary that contains {idx: sigobj}.

        Internal method of this class.
        """
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
        for (hashval, idlist) in self._hashval_to_idx.items():
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
            ident = self._idx_to_ident[idx]
            name = self._ident_to_name[ident]
            ss = SourmashSignature(mh, name=name)
            ss.into_frozen()

            if passes_all_picklists(ss, self.picklists):
                sigd[idx] = ss

        debug('=> {} signatures!', len(sigd))
        return sigd

    def find(self, search_fn, query, **kwargs):
        """
        Do a Jaccard similarity or containment search, yield results.

        Here 'search_fn' should be an instance of 'JaccardSearch'.

        As with SBTs, queries with higher scaled values than the database
        can still be used for containment search, but not for similarity
        search. See SBT.select(...) for details.

        Part of the Index protocol.
        """
        search_fn.check_is_compatible(query)

        # make sure we're looking at the same scaled value as database
        query_mh = query.minhash
        query_scaled = query_mh.scaled
        if self.scaled > query_scaled:
            query_mh = query_mh.downsample(scaled=self.scaled)
            query_scaled = query_mh.scaled
            prepare_subject = lambda x: x # identity
        else:
            prepare_subject = lambda subj: subj.downsample(scaled=query_scaled)

        # collect matching hashes for the query:
        c = Counter()
        query_hashes = set(query_mh.hashes)
        for hashval in query_hashes:
            idx_list = self._hashval_to_idx.get(hashval, [])
            for idx in idx_list:
                c[idx] += 1

        debug('number of matching signatures for hashes: {}', len(c))

        # for each match, in order of largest overlap,
        for idx, count in c.most_common():
            # pull in the hashes. This reconstructs & caches all input
            # minhashes, which is kinda memory intensive...!
            # NOTE: one future low-mem optimization could be to support doing
            # this piecemeal by iterating across all the hashes, instead.

            subj = self._signatures.get(idx)
            if subj is None:    # must be because of a picklist exclusion
                assert self.picklists
                continue

            subj_mh = prepare_subject(subj.minhash)

            # all numbers calculated after downsampling --
            query_size = len(query_mh)
            subj_size = len(subj_mh)
            shared_size = query_mh.count_common(subj_mh)
            total_size = len(query_mh + subj_mh)

            score = search_fn.score_fn(query_size, shared_size, subj_size,
                                       total_size)

            # CTB note to self: even with JaccardSearchBestOnly, this will
            # still iterate over & score all signatures. We should come
            # up with a protocol by which the JaccardSearch object can
            # signal that it is done, or something.
            # For example, see test_lca_jaccard_ordering, where
            # for containment we could be done early, but for Jaccard we
            # cannot.
            if search_fn.passes(score):
                if search_fn.collect(score, subj):
                    if passes_all_picklists(subj, self.picklists):
                        yield IndexSearchResult(score, subj, self.location)

    @cached_property
    def _lid_to_idx(self):
        """Connect lineage id lid (int) to idx set (set of ints).""

        Method specific to LCA databases.
        """
        d = defaultdict(set)
        for idx, lid in self._idx_to_lid.items():
            d[lid].add(idx)
        return d

    @cached_property
    def _idx_to_ident(self):
        """Connect idx (int) to ident (str).

        Method specific to LCA databases.
        """
        d = defaultdict(set)
        for ident, idx in self._ident_to_idx.items():
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
            notify(f'... loading database {format(db_name)}', end='\r')

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
        notify(f'loaded {len(dblist)} LCA databases. ksize={ksize}, scaled={scaled} moltype={moltype}')

    return dblist, ksize, scaled
