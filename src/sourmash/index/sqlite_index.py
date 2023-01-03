"""sqlite3 based Index, CollectionManifest, and LCA_Database
implementations.

These classes support a variety of flexible and fast on-disk storage,
search, and retrieval functions.

SqliteIndex stores full scaled signatures; sketches are stored as
reverse-indexed collections of hashes. Search is optimized via the
reverse index. Num and abund sketches are not supported. All scaled
values must be the same upon insertion. Multiple moltypes _are_
supported.

SqliteCollectionManifest provides a full implementation of the
manifest API. It can store details for all signature types. When used
as part of a SqliteIndex database, it does not support independent
insertion.

LCA_SqliteDatabase builds on top of SqliteIndex and LineageDB_Sqlite
(in the tax submodule) to provide a full on-disk implementation of
LCA_Database.

Using these classes
-------------------

These classes are fully integrated into sourmash loading.

Internally, use `sqlite_index.load_sqlite_index(...)` to load a specific
file; this will return the appropriate SqliteIndex, StandaloneManifestIndex,
or LCA_Database object.

Use `CollectionManifest.load_from_filename(...)` to load the manifest
directly as a manifest object.

Implementation Details
----------------------

SqliteIndex:

* Hashes with values above MAX_SQLITE_INT=2**63-1 are transformed into
  signed long longs upon insertion, and then back into ulong longs upon
  retrieval.

* Hash overlap is calculated via a SELECT.

* SqliteIndex relies on SqliteCollectionManifest for manifest functionality,
  including signature selection and picklists.

SqliteCollectionManifest:

* each object maintains info about whether it is being "managed" by a
  SqliteIndex class or not. If it is, `_insert_row(...)` cannot be
  called directly.

* `select(...)` operates directly with SQL queries, except for
  picklist selection, which involves inspect each manifest row in
  Python. In addition to being (much) simpler, this ends up being
  faster in some important real world situations, even for millions of
  rows!

* filter_on_rows and filter_on_columns also both operate in Python,
  not SQL.

* for this reason, the `locations()` method returns a superset of
  locations.  This is potentially very significant if you do a select
  with a picklist that ignores most sketches - the `locations()`
  method will ignore the picklist.

Limitations:

* all of these classes share a single connection object, and it could
  get confusing quickly if you simultaneously insert and query. We suggest
  separating creation and insertion. That having been said, these databases
  should work fine for many simultaneous queries; just don't write :).

"""
import time
import os
import sqlite3
from collections import defaultdict
import itertools

from bitstring import BitArray

from sourmash.index import Index
from sourmash.exceptions import IndexNotSupported
from sourmash import MinHash, SourmashSignature
from sourmash.index import IndexSearchResult, StandaloneManifestIndex
from sourmash.picklist import SignaturePicklist
from sourmash.logging import debug_literal
from sourmash import sqlite_utils

from sourmash.lca.lca_db import cached_property
from sourmash.manifest import BaseCollectionManifest

# converters for unsigned 64-bit ints: if over MAX_SQLITE_INT,
# convert to signed int.

MAX_SQLITE_INT = 2 ** 63 - 1
convert_hash_to = lambda x: BitArray(uint=x, length=64).int if x > MAX_SQLITE_INT else x
convert_hash_from = lambda x: BitArray(int=x, length=64).uint if x < 0 else x


def load_sqlite_index(filename, *, request_manifest=False):
    """Load a SqliteIndex, SqliteCollectionManifest, or LCA_SqliteDatabase.

    This is the main top-level API for loading an Index-like object. The logic
    is roughly:

    * does this database have both index and lineage tables? If so,
      return an LCA_SqliteDatabase.
    * if it only has an index, return a SqliteIndex.
    * if it only has a manifest, return a StandaloneManifestIndex.

    If you would like only a manifest, specify 'request_manifest=True'.
    """
    conn = sqlite_utils.open_sqlite_db(filename)

    if conn is None:
        debug_literal("load_sqlite_index: conn is None.")
        return

    c = conn.cursor()
    internal_d = sqlite_utils.get_sourmash_internal(c)

    is_index = False
    is_manifest = False
    is_lca_db = False

    if 'SqliteIndex' in internal_d:
        v = internal_d['SqliteIndex']
        if v != '1.0':
            raise IndexNotSupported
        is_index = True
        debug_literal("load_sqlite_index: it's an index!")

    if is_index and 'SqliteLineage' in internal_d:
        v = internal_d['SqliteLineage']
        if v != '1.0':
            raise IndexNotSupported

        is_lca_db = True
        debug_literal("load_sqlite_index: it's got a lineage table!")

    if 'SqliteManifest' in internal_d:
        v = internal_d['SqliteManifest']
        if v != '1.0':
            raise IndexNotSupported
        is_manifest = True
        debug_literal(f"load_sqlite_index: it's a manifest! request_manifest: {request_manifest}")

    # every Index is a Manifest!
    if is_index or is_lca_db:
        assert is_manifest

    idx = None
    if is_index and not request_manifest:
        conn.close()

        if is_lca_db:
            debug_literal("load_sqlite_index: returning LCA_SqliteDatabase")
            idx = LCA_SqliteDatabase.load(filename)
        else:
            debug_literal("load_sqlite_index: returning SqliteIndex")
            idx = SqliteIndex(filename)
    elif is_manifest:
        managed_by_index=False
        if is_index:
            assert request_manifest
            managed_by_index=True

        prefix = os.path.dirname(filename)
        mf = SqliteCollectionManifest(conn, managed_by_index=managed_by_index)
        idx = StandaloneManifestIndex(mf, filename, prefix=prefix)
        debug_literal("load_sqlite_index: returning StandaloneManifestIndex")

    return idx


class SqliteIndex(Index):
    is_database = True
    
    # NOTE: we do not need _signatures_with_internal for this class
    # because it supplies a manifest directly :tada:.

    def __init__(self, dbfile, *, sqlite_manifest=None, conn=None):
        "Constructor. 'dbfile' should be valid filename or ':memory:'."
        self.dbfile = dbfile

        # no connection? connect and/or create!
        if conn is None:
            conn = self._open(dbfile)

        # build me a SQLite manifest class to use for selection.
        if sqlite_manifest is None:
            sqlite_manifest = SqliteCollectionManifest(conn,
                                                       managed_by_index=True)
        self.manifest = sqlite_manifest
        self.conn = conn

        # set 'scaled'.
        c = self.conn.cursor()
        c.execute("SELECT DISTINCT scaled FROM sourmash_sketches")
        scaled_vals = c.fetchall()
        if len(scaled_vals) > 1:
            raise ValueError("this database has multiple scaled values, which is not currently allowed")

        if scaled_vals:
            self.scaled = scaled_vals[0][0]
        else:
            self.scaled = None

    @classmethod
    def _open(cls, dbfile, *, empty_ok=True):
        "Connect to existing SQLite database or create new."
        try:
            conn = sqlite3.connect(dbfile)
            c = conn.cursor()

            c.execute("PRAGMA cache_size=10000000")
            c.execute("PRAGMA synchronous = OFF")
            c.execute("PRAGMA journal_mode = MEMORY")
            c.execute("PRAGMA temp_store = MEMORY")

            if not empty_ok:
                c.execute("SELECT * FROM sourmash_hashes LIMIT 1")
                c.fetchone()
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            raise ValueError(f"cannot open '{dbfile}' as SqliteIndex database")

        return conn

    @classmethod
    def load(self, dbfile):
        "Load an existing SqliteIndex from dbfile."
        return SqliteIndex(dbfile)

    @classmethod
    def create(cls, dbfile, *, append=False):
        "Create a new SqliteIndex in dbfile."
        conn = cls._open(dbfile, empty_ok=True)
        cls._create_tables(conn.cursor(), ignore_exists=append)
        conn.commit()

        return cls(dbfile, conn=conn)

    @classmethod
    def _create_tables(cls, c, *, ignore_exists=False):
        "Create sqlite tables for SqliteIndex"
        try:
            sqlite_utils.add_sourmash_internal(c, 'SqliteIndex', '1.0')
            SqliteCollectionManifest._create_tables(c)

            c.execute("""
            CREATE TABLE IF NOT EXISTS sourmash_hashes (
               hashval INTEGER NOT NULL,
               sketch_id INTEGER NOT NULL,
               FOREIGN KEY (sketch_id) REFERENCES sourmash_sketches (id)
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS sourmash_hashval_idx ON sourmash_hashes (
               hashval,
               sketch_id
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS sourmash_hashval_idx2 ON sourmash_hashes (
               hashval
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS sourmash_sketch_idx ON sourmash_hashes (
               sketch_id
            )
            """
            )
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            if not ignore_exists:
                raise ValueError("cannot create SqliteIndex tables")

        return c

    def cursor(self):
        return self.conn.cursor()

    def close(self):
        self.conn.close()

    def commit(self):
        self.conn.commit()

    def __len__(self):
        return len(self.manifest)

    def insert(self, ss, *, cursor=None, commit=True):
        """
        Insert a signature into the sqlite database.

        If a cursor object is supplied, use that cursor instead of
        generating a new one.

        If 'commit' is True, commit after add; otherwise, do not.
        """
        if cursor:
            c = cursor
        else:
            c = self.conn.cursor()

        if ss.minhash.num:
            raise ValueError("cannot store 'num' signatures in SqliteIndex")
        if ss.minhash.track_abundance:
            raise ValueError("cannot store signatures with abundance in SqliteIndex")

        if self.scaled is not None and self.scaled != ss.minhash.scaled:
            raise ValueError(f"this database can only store scaled values={self.scaled}")
        elif self.scaled is None:
            self.scaled = ss.minhash.scaled

        # ok, first create and insert a manifest row
        row = BaseCollectionManifest.make_manifest_row(ss, None,
                                                       include_signature=False)
        self.manifest._insert_row(c, row, call_is_from_index=True)

        # retrieve ID of row for retrieving hashes:
        c.execute("SELECT last_insert_rowid()")
        sketch_id, = c.fetchone()

        # insert all the hashes
        hashes_to_sketch = []
        for h in ss.minhash.hashes:
            hh = convert_hash_to(h)
            hashes_to_sketch.append((hh, sketch_id))

        c.executemany("INSERT INTO sourmash_hashes (hashval, sketch_id) VALUES (?, ?)",
                      hashes_to_sketch)

        if commit:
            self.conn.commit()

    @property
    def location(self):
        return self.dbfile

    def signatures(self):
        "Return an iterator over all signatures in the Index object."
        for ss, loc in self.signatures_with_location():
            yield ss

    def signatures_with_location(self):
        "Return an iterator over tuples (signature, location) in the Index."
        c = self.conn.cursor()

        for ss, loc, iloc in self._load_sketches(c):
            yield ss, loc

    def save(self, *args, **kwargs):
        raise NotImplementedError

    def find(self, search_fn, query, **kwargs):
        search_fn.check_is_compatible(query)

        # check compatibility, etc.
        query_mh = query.minhash
        if self.scaled > query_mh.scaled:
            query_mh = query_mh.downsample(scaled=self.scaled)

        picklist = None
        if self.manifest.selection_dict:
            picklist = self.manifest.selection_dict.get('picklist')

        c1 = self.conn.cursor()
        c2 = self.conn.cursor()

        debug_literal('running _get_matching_sketches...')
        t0 = time.time()
        xx = self._get_matching_sketches(c1, query_mh.hashes,
                                         query_mh._max_hash)
        for sketch_id, n_matching_hashes in xx:
            debug_literal(f"...got sketch {sketch_id}, with {n_matching_hashes} matching hashes in {time.time() - t0:.2f}")
            #
            # first, estimate sketch size using sql results.
            #
            query_size = len(query_mh)
            subj_size = self._load_sketch_size(c2, sketch_id,
                                               query_mh._max_hash)
            total_size = query_size + subj_size - n_matching_hashes
            shared_size = n_matching_hashes

            score = search_fn.score_fn(query_size, shared_size, subj_size,
                                       total_size)

            debug_literal(f"APPROX RESULT: score={score} qsize={query_size}, ssize={subj_size} total={total_size} overlap={shared_size}")

            # do we pass?
            if not search_fn.passes(score):
                debug_literal(f"FAIL score={score}")

            # CTB if we are doing containment only, we could break loop here.
            # but for Jaccard, we must continue.
            # see 'test_sqlite_jaccard_ordering'

            if search_fn.passes(score):
                subj = self._load_sketch(c2, sketch_id)
                if search_fn.collect(score, subj):
                    if picklist is None or subj in picklist:
                        yield IndexSearchResult(score, subj, self.location)

    def _select(self, *, num=0, track_abundance=False, **kwargs):
        "Run a select! This just modifies the manifest."
        # check SqliteIndex specific conditions on the 'select'
        if num:
            raise ValueError("cannot select on 'num' in SqliteIndex")
        if track_abundance:
            raise ValueError("cannot store or search signatures with abundance")
        # create manifest if needed
        manifest = self.manifest
        if manifest is None:
            manifest = SqliteCollectionManifest(self.conn,
                                                managed_by_index=True)

        # modify manifest
        manifest = manifest.select_to_manifest(**kwargs)

        return manifest

    def select(self, *args, **kwargs):
        sqlite_manifest = self._select(*args, **kwargs)

        # return a new SqliteIndex with a new manifest, but same old conn.
        return SqliteIndex(self.dbfile,
                           sqlite_manifest=sqlite_manifest,
                           conn=self.conn)

    #
    # Actual SQL queries, etc.
    #

    def _load_sketch_size(self, c1, sketch_id, max_hash):
        "Get sketch size for given sketch, downsampled by max_hash."
        if max_hash <= MAX_SQLITE_INT:
            c1.execute("""
            SELECT COUNT(hashval) FROM sourmash_hashes
            WHERE sketch_id=? AND hashval >= 0 AND hashval <= ?""",
                       (sketch_id, max_hash))
        else:
            c1.execute('SELECT COUNT(hashval) FROM sourmash_hashes WHERE sketch_id=?',
                       (sketch_id,))

        n_hashes, = c1.fetchone()
        return n_hashes

    def _load_sketch(self, c, sketch_id, *, match_scaled=None):
        "Load an individual sketch. If match_scaled is set, downsample."

        start = time.time()
        c.execute("""
        SELECT id, name, scaled, ksize, filename, moltype, seed
        FROM sourmash_sketches WHERE id=?""", (sketch_id,))
        debug_literal(f"load sketch {sketch_id}: got sketch info in {time.time() - start:.2f}")

        sketch_id, name, scaled, ksize, filename, moltype, seed = c.fetchone()
        if match_scaled is not None:
            scaled = max(scaled, match_scaled)

        is_protein = 1 if moltype=='protein' else 0
        is_dayhoff = 1 if moltype=='dayhoff' else 0
        is_hp = 1 if moltype=='hp' else 0

        mh = MinHash(n=0, ksize=ksize, scaled=scaled, seed=seed,
                     is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp)


        template_values = [sketch_id]

        hash_constraint_str = ""
        max_hash = mh._max_hash
        if max_hash <= MAX_SQLITE_INT:
            hash_constraint_str = "sourmash_hashes.hashval >= 0 AND sourmash_hashes.hashval <= ? AND"
            template_values.insert(0, max_hash)
        else:
            debug_literal('NOT EMPLOYING hash_constraint_str')

        debug_literal(f"finding hashes for sketch {sketch_id} in {time.time() - start:.2f}")
        c.execute(f"SELECT hashval FROM sourmash_hashes WHERE {hash_constraint_str} sourmash_hashes.sketch_id=?", template_values)

        debug_literal(f"loading hashes for sketch {sketch_id} in {time.time() - start:.2f}")
        for hashval, in c:
            hh = convert_hash_from(hashval)
            mh.add_hash(hh)

        debug_literal(f"done loading sketch {sketch_id} {time.time() - start:.2f})")

        return SourmashSignature(mh, name=name, filename=filename)

    def _load_sketches(self, c):
        "Load sketches based on manifest _id column."
        for row in self.manifest.rows:
            sketch_id = row['_id']
            assert row['num'] == 0

            moltype = row['moltype']
            is_protein = 1 if moltype=='protein' else 0
            is_dayhoff = 1 if moltype=='dayhoff' else 0
            is_hp = 1 if moltype=='hp' else 0

            ksize = row['ksize']
            scaled = row['scaled']
            seed = row['seed']

            mh = MinHash(n=0, ksize=ksize, scaled=scaled, seed=seed,
                         is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp)

            c.execute("SELECT hashval FROM sourmash_hashes WHERE sketch_id=?",
                       (sketch_id,))

            for hashval, in c:
                mh.add_hash(convert_hash_from(hashval))

            ss = SourmashSignature(mh, name=row['name'],
                                   filename=row['filename'])
            yield ss, self.dbfile, sketch_id

    def _get_matching_sketches(self, c, hashes, max_hash):
        """
        For hashvals in 'hashes', retrieve all matching sketches,
        together with the number of overlapping hashes for each sketch.

        CTB: we do not use sqlite manifest conditions on this select,
        because it slows things down in practice.
        """
        c.execute("DROP TABLE IF EXISTS sourmash_hash_query")
        c.execute("CREATE TEMPORARY TABLE sourmash_hash_query (hashval INTEGER PRIMARY KEY)")

        hashvals = [ (convert_hash_to(h),) for h in hashes ]
        c.executemany("INSERT OR IGNORE INTO sourmash_hash_query (hashval) VALUES (?)",
                      hashvals)

        #
        # set up SELECT conditions
        #

        conditions = []
        template_values = []

        # downsample? => add to conditions
        max_hash = min(max_hash, max(hashes))
        if max_hash <= MAX_SQLITE_INT:
            select_str = "sourmash_hashes.hashval >= 0 AND sourmash_hashes.hashval <= ?"
            conditions.append(select_str)
            template_values.append(max_hash)

        # format conditions
        conditions.append('sourmash_hashes.hashval=sourmash_hash_query.hashval')
        conditions = " AND ".join(conditions)

        c.execute(f"""
        SELECT DISTINCT sourmash_hashes.sketch_id,COUNT(sourmash_hashes.hashval) as CNT
        FROM sourmash_hashes, sourmash_hash_query
        WHERE {conditions}
        GROUP BY sourmash_hashes.sketch_id ORDER BY CNT DESC
        """, template_values)

        return c


class SqliteCollectionManifest(BaseCollectionManifest):
    """
    A SQLite-based manifest, used both for SqliteIndex and as a standalone
    manifest class.

    This class serves two purposes:
    * first, it is a fast, on-disk manifest that can be used in place of
      CollectionManifest.
    * second, it can be included within a SqliteIndex (which stores hashes
      too). In this case, however, new entries must be inserted by SqliteIndex
      rather than directly in this class.

    In the latter case, the SqliteCollectionManifest is created with
    managed_by_index set to True.
    """
    def __init__(self, conn, *, selection_dict=None, managed_by_index=False):
        """
        Here, 'conn' should already be connected and configured.

        Use 'create(filename)' to create a new database.

        Use 'create_from_manifest(filename, manifest) to create a new db
        from an existing manifest object.

        Use 'load_from_filename' to load from file.
        """
        assert conn is not None
        self.conn = conn
        self.selection_dict = selection_dict
        self.managed_by_index = managed_by_index
        self._num_rows = None

    @classmethod
    def create(cls, filename):
        "Connect to 'filename' and create the tables as a standalone manifest."
        conn = sqlite3.connect(filename)
        cursor = conn.cursor()
        cls._create_tables(cursor)
        return cls(conn)

    @classmethod
    def create_or_open(cls, filename):
        "Connect to 'filename' and create tables if not exist."
        conn = sqlite3.connect(filename)
        cursor = conn.cursor()
        try:
            cls._create_tables(cursor)
        except sqlite3.OperationalError:
            pass
        return cls(conn)

    @classmethod
    def load_from_manifest(cls, manifest, *, dbfile=":memory:", append=False):
        "Create a new sqlite manifest from an existing manifest object."
        return cls._create_manifest_from_rows(manifest.rows, location=dbfile,
                                              append=append)

    @classmethod
    def create_manifest(cls, locations_iter, *, include_signature=False):
        """Create a manifest from an iterator that yields (ss, location)

        Stores signatures in manifest rows by default.

        Note: do NOT catch exceptions here, so this passes through load excs.
        Note: this method ignores 'include_signature'.
        """
        def rows_iter():
            for ss, location in locations_iter:
                row = cls.make_manifest_row(ss, location,
                                            include_signature=False)
                yield row

        return cls._create_manifest_from_rows(rows_iter())

    @classmethod
    def _create_tables(cls, cursor):
        "Create the manifest table."
        # this is a class method so that it can be used by SqliteIndex to
        # create manifest-compatible tables.

        sqlite_utils.add_sourmash_internal(cursor, 'SqliteManifest', '1.0')
        cursor.execute("""
        CREATE TABLE sourmash_sketches
          (id INTEGER PRIMARY KEY,
           name TEXT,
           num INTEGER NOT NULL,
           scaled INTEGER NOT NULL,
           ksize INTEGER NOT NULL,
           filename TEXT,
           moltype TEXT NOT NULL,
           with_abundance BOOLEAN NOT NULL,
           md5sum TEXT NOT NULL,
           seed INTEGER NOT NULL,
           n_hashes INTEGER NOT NULL,
           internal_location TEXT,
        UNIQUE(internal_location, md5sum)
        )
        """)

    def add_row(self, row):
        c = self.conn.cursor()
        self._insert_row(c, row)

    def _insert_row(self, cursor, row, *, call_is_from_index=False):
        "Insert a new manifest row."
        # check - is this manifest managed by SqliteIndex? If so, prevent
        # insertions unless SqliteIndex is the one calling it.
        if self.managed_by_index and not call_is_from_index:
            raise Exception("must use SqliteIndex.insert to add to this manifest")

        row = dict(row)
        if 'seed' not in row:
            row['seed'] = 42

        cursor.execute("""
        INSERT OR IGNORE INTO sourmash_sketches
          (name, num, scaled, ksize, filename, md5sum, moltype,
           seed, n_hashes, with_abundance, internal_location)
        VALUES (:name, :num, :scaled, :ksize, :filename, :md5,
                :moltype, :seed, :n_hashes, :with_abundance,
                :internal_location)""", row)

        self._num_rows = None   # reset cache

    def __bool__(self):
        "Is this manifest empty?"
        if self._num_rows is not None:
            return bool(self._num_rows)

        try:
            next(iter(self.rows))
            return True
        except StopIteration:
            return False

    def __eq__(self, other):
        "Check equality on a row-by-row basis. May fail on out-of-order rows."
        for (a, b) in itertools.zip_longest(self.rows, other.rows):
            # ignore non-required keys.
            for k in self.required_keys:
                if a[k] != b[k]:
                    return False

        return True

    def __len__(self):
        "Number of rows."

        # can we use cached value?
        if self._num_rows is not None:
            return self._num_rows

        # self.rows is a generator, so can't use 'len'
        self._num_rows = sum(1 for _ in self.rows)
        return self._num_rows

    def __iadd__(self, other):
        c = self.conn.cursor()
        for row in other.rows:
            self._insert_row(c, row)
        return self

    def __add__(self, other):
        new_mf = self.create(":memory:")
        new_mf += self
        new_mf += other
        return new_mf

    def close(self):
        self.conn.commit()

    def _make_select(self):
        """Build a set of SQL SELECT conditions and matching value tuple
        that can be used to select the right sketches from the
        database.

        Returns a triple 'conditions', 'values', and 'picklist'.
        'conditions' is a list that should be joined with 'AND'.

        The picklist is simply retrieved from the selection dictionary.
        """
        conditions = []
        values = []
        picklist = None
        if self.selection_dict:
            select_d = self.selection_dict
            if 'ksize' in select_d and select_d['ksize']:
                conditions.append("sourmash_sketches.ksize = ?")
                values.append(select_d['ksize'])
            if 'num' in select_d and select_d['num'] > 0:
                conditions.append("sourmash_sketches.num > 0")
            if 'scaled' in select_d and select_d['scaled'] > 0:
                conditions.append("sourmash_sketches.scaled > 0")
            if 'containment' in select_d and select_d['containment']:
                conditions.append("sourmash_sketches.scaled > 0")
            if 'moltype' in select_d and select_d['moltype'] is not None:
                moltype = select_d['moltype']
                assert moltype in ('DNA', 'protein', 'dayhoff', 'hp'), moltype
                conditions.append(f"sourmash_sketches.moltype = '{moltype}'")

            picklist = select_d.get('picklist')

        return conditions, values, picklist

    def select_to_manifest(self, **kwargs):
        "Create a new SqliteCollectionManifest with the given select args."
        # Pass along all the selection kwargs to a new instance
        if self.selection_dict:
            debug_literal("sqlite manifest: merging selection dicts")
            # combine selects...
            d = dict(self.selection_dict)
            for k, v in kwargs.items():
                if k in d:
                    if d[k] is not None and d[k] != v:
                        raise ValueError(f"incompatible select on '{k}'")
                d[k] = v
            kwargs = d

        new_mf = SqliteCollectionManifest(self.conn, selection_dict=kwargs)

        # if picklist, make sure we fill in 'found'.
        picklist = kwargs.get('picklist')
        if picklist is not None:
            debug_literal("sqlite manifest: iterating through picklist")
            _ = len(self)       # this forces iteration through rows.

        return new_mf

    @property
    def rows(self):
        "Return rows that match the selection."
        c1 = self.conn.cursor()

        conditions, values, picklist = self._make_select()
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        debug_literal(f"sqlite manifest rows: executing select with '{conditions}'")
        c1.execute(f"""
        SELECT id, name, md5sum, num, scaled, ksize, filename, moltype,
        seed, n_hashes, internal_location FROM sourmash_sketches {conditions}
        """, values)

        debug_literal("sqlite manifest: entering row yield loop")
        for (_id, name, md5sum, num, scaled, ksize, filename, moltype,
             seed, n_hashes, iloc) in c1:
            row = dict(num=num, scaled=scaled, name=name, filename=filename,
                       n_hashes=n_hashes, with_abundance=0, ksize=ksize,
                       md5=md5sum, internal_location=iloc,
                       moltype=moltype, md5short=md5sum[:8],
                       seed=seed, _id=_id)
            if picklist is None or picklist.matches_manifest_row(row):
                yield row

    def filter_rows(self, row_filter_fn):
        """Create a new manifest filtered through row_filter_fn.

        This is done in memory, inserting each row one at a time.
        """
        def rows_iter():
            for row in self.rows:
                if row_filter_fn(row):
                    yield row

        return self._create_manifest_from_rows(rows_iter())

    def filter_on_columns(self, col_filter_fn, col_names):
        "Create a new manifest based on column matches."
        def row_filter_fn(row):
            x = [ row[col] for col in col_names if row[col] is not None ]
            return col_filter_fn(x)
        return self.filter_rows(row_filter_fn)

    def locations(self):
        """Return all possible locations for signatures.

        CTB: this may be a (big) superset of locations, if picklists are used.
        See test_sqlite_manifest_locations.

        Use set(row['internal_locations'] for row in self.rows)
        if you want an exact set of locations; will be slow for big manifests
        tho.
        """
        c1 = self.conn.cursor()

        conditions, values, picklist = self._make_select()
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        c1.execute(f"""
        SELECT DISTINCT internal_location FROM sourmash_sketches {conditions}
        """, values)

        return ( iloc for iloc, in c1 )

    def __contains__(self, ss):
        "Check to see if signature 'ss' is in this manifest."
        md5 = ss.md5sum()

        c = self.conn.cursor()
        c.execute('SELECT COUNT(*) FROM sourmash_sketches WHERE md5sum=?',
                  (md5,))
        val, = c.fetchone()

        if bool(val):
            picklist = self.picklist
            return picklist is None or ss in self.picklist
        return False

    @property
    def picklist(self):
        "Return the picklist, if any."
        if self.selection_dict:
            return self.selection_dict.get('picklist')
        return None

    def to_picklist(self):
        "Convert this manifest to a picklist."
        pickset = set()
        for row in self.rows:
            pickset.add(row['md5'])

        picklist = SignaturePicklist('md5')
        picklist.pickset = pickset
        return picklist

    @classmethod
    def _create_manifest_from_rows(cls, rows_iter, *, location=":memory:",
                                   append=False):
        """Create a SqliteCollectionManifest from a rows iterator.

        Internal utility function.

        CTB: should enable converting in-memory sqlite db to on-disk,
        probably with sqlite3 'conn.backup(...)' function.
        """
        try:
            mf = cls.create(location)
        except (sqlite3.OperationalError, sqlite3.DatabaseError) as exc:
            if not append:
                raise Exception(f"cannot create sqlite3 db at '{location}'; exception: {str(exc)}")
            db = load_sqlite_index(location, request_manifest=True)
            mf = db.manifest

        cursor = mf.conn.cursor()

        for row in rows_iter:
            mf._insert_row(cursor, row)

        mf.conn.commit()
        return mf


class LCA_SqliteDatabase(SqliteIndex):
    """
    A wrapper class for SqliteIndex + lineage db => LCA_Database functionality.
    """
    is_database = True

    def __init__(self, dbfile, *, lineage_db=None, sqlite_manifest=None):
        # CTB note: we need to let SqliteIndex open dbfile here, so can't
        # just pass in a conn.
        super().__init__(dbfile, sqlite_manifest=sqlite_manifest)

        c = self.conn.cursor()

        c.execute('SELECT DISTINCT ksize, moltype FROM sourmash_sketches')
        res = list(c)
        if len(res) > 1:
            raise TypeError("can only have one ksize & moltype in an LCA_SqliteDatabase")
        if len(res) == 0:
            raise ValueError("cannot load an LCA_SqliteDatabase")

        self.ksize, self.moltype = res[0]
        debug_literal(f"setting ksize and moltype to {self.ksize}, {self.moltype}")

        if lineage_db is not None:
            self.lineage_db = lineage_db

            ## the below is done once, but could be implemented as something
            ## ~dynamic.
            self._build_index()

    @classmethod
    def load(cls, filename):
        "Load LCA_SqliteDatabase from a single file."
        from sourmash.tax.tax_utils import LineageDB_Sqlite

        # first, load the SqliteIndex:
        try:
            debug_literal("sqlite_index: loading LCA_SqliteDatabase as SqliteIndex.")
            obj = cls(filename)
        except sqlite3.OperationalError:
            raise ValueError(f"cannot open '{filename}' as a SQLite index.")

        # now, toss in the lineage DB.
        lineage_db = LineageDB_Sqlite(obj.conn)
        obj.lineage_db = lineage_db
        obj._build_index()

        return obj

    @classmethod
    def create(cls, filename, idx, lineage_db):
        "Create a LCA_SqliteDatabase in a single file from existing idx/ldb."
        from sourmash.tax.tax_utils import MultiLineageDB

        # first, save/create signatures...
        sqlidx = SqliteIndex.create(filename)

        for ss in idx.signatures():
            sqlidx.insert(ss)

        # now, save the lineage_db into the same database
        out_lineage_db = MultiLineageDB()
        out_lineage_db.add(lineage_db)
        out_lineage_db._save_sqlite(None, conn=sqlidx.conn)

        # and voila! return, I guess?
        return cls.load(filename)

    def _build_index(self):
        "Rebuild the mappings that support identifier <-> lineage."
        mf = self.manifest
        lineage_db = self.lineage_db

        ident_to_idx = {}
        next_lid = 0
        idx_to_lid = {}
        lineage_to_lid = {}
        lid_to_lineage = {}

        for row in mf.rows:
            name = row['name']
            if name:
                # this is a bit of a hack. we try identifiers _with_ and
                # _without_ versions, and take whichever works. There is
                # definitely a better way to do this, but I can't think
                # of one right now.
                ident = name.split(' ')[0]

                lineage = lineage_db.get(ident) # try with identifier version
                if lineage is None:             # nope - remove version.x
                    ident = name.split('.')[0]
                    lineage = lineage_db.get(ident)

                idx = row['_id'] # this is only present in sqlite manifests.
                ident_to_idx[ident] = idx

                if lineage:
                    lid = lineage_to_lid.get(lineage)

                    # manufacture new lid?
                    if lid is None:
                        lid = next_lid
                        next_lid += 1

                        lineage_to_lid[lineage] = lid
                        lid_to_lineage[lid] = lineage

                    # assign idx <-> lid
                    idx_to_lid[idx] = lid

        self.ident_to_idx = ident_to_idx
        self.idx_to_lid = idx_to_lid
        self.lid_to_lineage = lid_to_lineage

    # prevent insertions
    def insert(self, *args, **kwargs):
        raise NotImplementedError

    # return correct type on select
    def select(self, *args, **kwargs):
        sqlite_manifest = self._select(*args, **kwargs)

        return LCA_SqliteDatabase(self.dbfile,
                                  sqlite_manifest=sqlite_manifest,
                                  lineage_db=self.lineage_db)

    ### LCA_Database API/protocol.

    def downsample_scaled(self, scaled):
        "Downsample the scaled for querying."
        if scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        # CTB: maybe return a new LCA_Database? Right now this isn't how
        # the lca_db protocol works tho.
        self.scaled = scaled

    def get_lineage_assignments(self, hashval, *, min_num=None):
        """
        Get a list of lineages for this hashval.
        """
        x = []

        idx_list = self.hashval_to_idx.get(hashval, [])
        if min_num is None or len(idx_list) >= min_num:
            for idx in idx_list:
                lid = self.idx_to_lid.get(idx, None)
                if lid is not None:
                    lineage = self.lid_to_lineage[lid]
                    x.append(lineage)

        return x

    @cached_property
    def idx_to_ident(self):
        "Map individual idx to ident."
        d = defaultdict(set)
        for ident, idx in self.ident_to_idx.items():
            assert idx not in d
            d[idx] = ident
        return d

    @property
    def hashval_to_idx(self):
        "Dynamically interpret the SQL 'hashes' table like it's a dict."
        return _SqliteIndexHashvalToIndex(self)

    @property
    def hashvals(self):
        "Return all hashvals"
        return iter(_SqliteIndexHashvalToIndex(self))

    def get_identifiers_for_hashval(self, hashval):
        "Return identifiers associated with this hashval"
        idxlist = self.hashval_to_idx[hashval]
        for idx in idxlist:
            yield self.idx_to_ident[idx]


class _SqliteIndexHashvalToIndex:
    """
    Internal wrapper class to retrieve keys and key/value pairs for 
    hashval -> [ list of idx ].
    """
    def __init__(self, sqlidx):
        self.sqlidx = sqlidx

    def __iter__(self):
        "Get all hashvals."
        c = self.sqlidx.conn.cursor()
        c.execute('SELECT DISTINCT hashval FROM sourmash_hashes')
        for hashval, in c:
            yield hashval

    def get(self, key, dv=None):
        "Retrieve idxlist for a given hash."
        sqlidx = self.sqlidx
        c = sqlidx.cursor()

        hh = convert_hash_to(key)

        c.execute('SELECT sketch_id FROM sourmash_hashes WHERE hashval=?',
                  (hh,))

        x = [ convert_hash_from(h) for h, in c ]
        return x or dv

    def __getitem__(self, key):
        "Retrieve idxlist for a given hash; raise KeyError if not present."
        v = self.get(key)
        if v is None:
            raise KeyError(key)
        return v
