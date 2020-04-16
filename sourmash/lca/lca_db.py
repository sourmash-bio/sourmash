"LCA database class and utilities."

from __future__ import print_function, division
import json
import gzip
from collections import OrderedDict, defaultdict, Counter
import functools

from sourmash._minhash import get_max_hash_for_scaled
from sourmash.logging import notify, error, debug
from sourmash.index import Index


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
    Wrapper class for taxonomic database.

    obj.ident_to_idx: key 'identifier' to 'idx'
    obj.idx_to_lid: key 'idx' to 'lid'
    obj.lid_to_lineage: key 'lid' to tuple of LineagePair objects
    obj.hashval_to_idx: key 'hashval' => set('idx')
    """
    def __init__(self, ksize, scaled):
        self.ksize = int(ksize)
        self.scaled = int(scaled)
        self.filename = None

        self._next_index = 0
        self._next_lid = 0
        self.ident_to_name = {}
        self.ident_to_idx = {}
        self.idx_to_lid = {}
        self.lineage_to_lid = {}
        self.lid_to_lineage = {}
        self.hashval_to_idx = defaultdict(set)

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

    def insert_signature(self, sig, ident=None, lineage=None): # @CTB -> insert
        """Add a new signature into the LCA database.

        Takes optional arguments 'ident' and 'lineage'.

        'ident' must be a unique string identifer across this database;
        if not specified, the signature name (sig.name()) is used.

        'lineage', if specified, must contain a tuple of LineagePair objects.
        """
        if ident is None:
            ident = sig.name()

        if ident in self.ident_to_name:
            raise ValueError("signature {} is already in this LCA db.".format(ident))

        # store full name
        self.ident_to_name[ident] = sig.name()

        # identifier -> integer index (idx)
        idx = self._get_ident_index(ident, fail_on_duplicate=True)
        if lineage:
            # (LineagePairs*) -> integer lineage ids (lids)
            lid = self._get_lineage_id(lineage)

            # map idx to lid as well.
            self.idx_to_lid[idx] = lid

        # downsample to specified scaled; this has the side effect of
        # making sure they're all at the same scaled value!
        minhash = sig.minhash.downsample_scaled(self.scaled)

        for hashval in minhash.get_mins():
            self.hashval_to_idx[hashval].add(idx)

        return lineage

    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def signatures(self):
        from .. import SourmashSignature
        for v in self._signatures.values():
            yield SourmashSignature(v)

    @classmethod
    def load(cls, db_name):
        from .lca_utils import taxlist, LineagePair

        "Load from a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
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

            if version != '2.0' or 'lid_to_lineage' not in load_d:
                raise ValueError("Error! This is an old-style LCA DB. You'll need to build or download a newer one.")

            ksize = int(load_d['ksize'])
            scaled = int(load_d['scaled'])

            db = cls(ksize, scaled)

            # convert lineage_dict to proper lineages (tuples of LineagePairs)
            lid_to_lineage_2 = load_d['lid_to_lineage']
            lid_to_lineage = {}
            for k, v in lid_to_lineage_2.items():
                v = dict(v)
                vv = []
                for rank in taxlist():
                    name = v.get(rank, '')
                    vv.append(LineagePair(rank, name))

                lid_to_lineage[int(k)] = tuple(vv)
            db.lid_to_lineage = lid_to_lineage

            # convert hashval -> lineage index keys to integers (looks like
            # JSON doesn't have a 64 bit type so stores them as strings)
            hashval_to_idx_2 = load_d['hashval_to_idx']
            hashval_to_idx = {}

            for k, v in hashval_to_idx_2.items():
                hashval_to_idx[int(k)] = v
            db.hashval_to_idx = hashval_to_idx

            db.ident_to_name = load_d['ident_to_name']
            db.ident_to_idx = load_d['ident_to_idx']

            db.idx_to_lid = {}
            for k, v in load_d['idx_to_lid'].items():
                db.idx_to_lid[int(k)] = v

        db.filename = db_name

        return db

    def save(self, db_name):
        "Save to a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '2.0'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'
            save_d['ksize'] = self.ksize
            save_d['scaled'] = self.scaled

            # convert lineage internals from tuples to dictionaries
            d = OrderedDict()
            for k, v in self.lid_to_lineage.items():
                d[k] = dict([ (vv.rank, vv.name) for vv in v ])
            save_d['lid_to_lineage'] = d

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_to_idx'] = \
               dict((k, list(v)) for (k, v) in self.hashval_to_idx.items())

            save_d['ident_to_name'] = self.ident_to_name
            save_d['ident_to_idx'] = self.ident_to_idx
            save_d['idx_to_lid'] = self.idx_to_lid
            save_d['lid_to_lineage'] = self.lid_to_lineage
            
            json.dump(save_d, fp)

    def search(self, query, *args, **kwargs):
        # check arguments
        if 'threshold' not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs['threshold']
        do_containment = kwargs.get('do_containment', False)
        ignore_abundance = kwargs.get('ignore_abundance', False)
        mh = query.minhash
        if ignore_abundance:
            mh.track_abundance = False

        results = []
        for x in self.find_signatures(mh, threshold, do_containment):
            (score, match, filename) = x
            results.append((score, match, filename))

        results.sort(key=lambda x: -x[0])
        return results

    def gather(self, query, *args, **kwargs):
        if not query.minhash:
            return []

        results = []
        threshold_bp = kwargs.get('threshold_bp', 0.0)
        threshold = threshold_bp / (len(query.minhash) * self.scaled)
        for x in self.find_signatures(query.minhash, threshold,
                                      containment=True, ignore_scaled=True):
            (score, match, filename) = x
            if score:
                results.append((score, match, filename))
                break

        return results

    def insert(self, node):
        raise NotImplementedError

    def find(self, search_fn, *args, **kwargs):
        raise NotImplementedError

    def downsample_scaled(self, scaled):
        """
        Downsample to the provided scaled value, i.e. eliminate all hashes
        that don't fall in the required range.

        NOTE: we probably need to invalidate some of the dynamically
        calculated members of this object, like _signatures, when we do this.
        But we aren't going to right now.
        """
        if scaled == self.scaled:
            return
        elif scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        max_hash = get_max_hash_for_scaled(scaled)
        new_hashvals = {}
        for k, v in self.hashval_to_idx.items():
            if k < max_hash:
                new_hashvals[k] = v
        self.hashval_to_idx = new_hashvals
        self.scaled = scaled

    def get_lineage_assignments(self, hashval):
        """
        Get a list of lineages for this hashval.
        """
        x = []

        idx_list = self.hashval_to_idx.get(hashval, [])
        for idx in idx_list:
            lid = self.idx_to_lid.get(idx, None)
            if lid is not None:
                lineage = self.lid_to_lineage[lid]
                x.append(lineage)

        return x

    @cached_property
    def _signatures(self):
        "Create a _signatures member dictionary that contains {idx: minhash}."
        from .. import MinHash

        minhash = MinHash(n=0, ksize=self.ksize, scaled=self.scaled)

        debug('creating signatures for LCA DB...')
        sigd = defaultdict(minhash.copy_and_clear)
        temp_vals = defaultdict(list)

        for (k, v) in self.hashval_to_idx.items():
            for vv in v:
                temp_hashes = temp_vals[vv]
                temp_hashes.append(k)

                # 50 is an arbitrary number. If you really want
                # to micro-optimize, list is resized and grow in this pattern:
                # 0, 4, 8, 16, 25, 35, 46, 58, 72, 88, ...
                # (from https://github.com/python/cpython/blob/b2b4a51f7463a0392456f7772f33223e57fa4ccc/Objects/listobject.c#L57)
                if len(temp_hashes) > 50:
                    sigd[vv].add_many(temp_hashes)

                    # Sigh, python 2... when it goes away,
                    # we can do `temp_hashes.clear()` instead.
                    del temp_vals[vv]

        # We loop temp_vals again to add any remainder hashes
        # (each list of hashes is smaller than 50 items)
        for sig, vals in temp_vals.items():
            sigd[sig].add_many(vals)

        debug('=> {} signatures!', len(sigd))
        return sigd

    def find_signatures(self, minhash, threshold, containment=False,
                       ignore_scaled=False):
        """
        Do a Jaccard similarity or containment search.
        """
        # make sure we're looking at the same scaled value as database
        if self.scaled > minhash.scaled:
            minhash = minhash.downsample_scaled(self.scaled)
        elif self.scaled < minhash.scaled and not ignore_scaled:
            # note that containment can be calculated w/o matching scaled.
            raise ValueError("lca db scaled is {} vs query {}; must downsample".format(self.scaled, minhash.scaled))

        query_mins = set(minhash.get_mins())

        # collect matching hashes:
        c = Counter()
        for hashval in query_mins:
            idx_list = self.hashval_to_idx.get(hashval, [])
            for idx in idx_list:
                c[idx] += 1

        debug('number of matching signatures for hashes: {}', len(c))

        for idx, count in c.items():
            ident = self.idx_to_ident[idx]
            name = self.ident_to_name[ident]

            match_mh = self._signatures[idx]
            match_size = len(match_mh)

            debug('count: {}; query_mins: {}; match size: {}',
                  count, len(query_mins), match_size)

            if containment:
                score = count / len(query_mins)
            else:
                score = count / (len(query_mins) + match_size - count)

            debug('score: {} (containment? {}), threshold: {}',
                  score, containment, threshold)

            if score >= threshold:
                from .. import SourmashSignature
                match_sig = SourmashSignature(match_mh, name=name)

                yield score, match_sig, self.filename

    @cached_property
    def lineage_to_lids(self):
        d = defaultdict(set)
        for lid, lineage in self.lid_to_lineage.items():
            d[lineage].add(lid)
        return d

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

        dblist.append(lca_db)

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()

    if verbose:
        notify(u'\r\033[K', end=u'')
        notify('loaded {} LCA databases. ksize={}, scaled={}', len(dblist),
               ksize, scaled)

    return dblist, ksize, scaled


