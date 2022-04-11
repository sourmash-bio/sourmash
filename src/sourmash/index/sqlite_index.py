"""Provide SqliteIndex, a sqlite3-based Index class for storing and searching
sourmash signatures.

Note that SqliteIndex supports both storage and fast _search_ of scaled
signatures, via a reverse index.

Features and limitations:

* Currently we try to maintain only one database connection.  It's not 100%
  clear what happens if the same database is opened by multiple independent
  processes and one or more of them write to it. It should all work, because
  SQL...

* Unlike LCA_Database, SqliteIndex supports multiple ksizes and moltypes.

* SqliteIndex does not support 'num' signatures. It could store them
  easily, but since it cannot search them properly with 'find', we've
  omitted them.

* Likewise, SqliteIndex does not support 'abund' signatures because it cannot
  search them (just like SBTs cannot).

* @CTB document manifest stuff.

CTB consider:
* a SqliteIndex sqldb can store taxonomy table just fine. Is there any
    extra support that might be worthwhile?

* do we want to prevent storage of scaled=1 sketches and then
  dispense with the MAX_SQLITE_INT stuff? It's kind of a nice hack :laugh:

TODO:
@CTB don't do constraints if scaleds are equal?
"""
import time
import os
import sqlite3
from collections import Counter, defaultdict
from collections.abc import Mapping

from bitstring import BitArray

from sourmash.index import Index
import sourmash
from sourmash import MinHash, SourmashSignature
from sourmash.index import IndexSearchResult, StandaloneManifestIndex
from sourmash.picklist import PickStyle, SignaturePicklist
from sourmash.logging import debug_literal, notify
from sourmash import sqlite_utils

from sourmash.lca.lca_db import cached_property
from sourmash.manifest import BaseCollectionManifest, CollectionManifest
from sourmash.logging import debug_literal

# converters for unsigned 64-bit ints: if over MAX_SQLITE_INT,
# convert to signed int.

MAX_SQLITE_INT = 2 ** 63 - 1
convert_hash_to = lambda x: BitArray(uint=x, length=64).int if x > MAX_SQLITE_INT else x
convert_hash_from = lambda x: BitArray(int=x, length=64).uint if x < 0 else x

picklist_transforms = dict(
    name=lambda x: x,
    ident=lambda x: x + ' %',
    identprefix=lambda x: x + '%',
    md5short=lambda x: x[:8] + '%',
    md5prefix8=lambda x: x[:8] + '%',
    md5=lambda x: x,
    )

picklist_selects = dict(
    name='INSERT INTO pickset SELECT id FROM sketches WHERE name=?',
    ident='INSERT INTO pickset SELECT id FROM sketches WHERE name LIKE ?',
    identprefix='INSERT INTO pickset SELECT id FROM sketches WHERE name LIKE ?',
    md5short='INSERT INTO pickset SELECT id FROM sketches WHERE md5sum LIKE ?',
    md5prefix8='INSERT INTO pickset SELECT id FROM sketches WHERE md5sum LIKE ?',
    md5='INSERT INTO pickset SELECT id FROM sketches WHERE md5sum=?',
    )


# @CTB write tests that cross-product the various types.
def load_sqlite_file(filename, *, request_manifest=False):
    "Load a SqliteIndex or a SqliteCollectionManifest from a sqlite file."
    conn = sqlite_utils.open_sqlite_db(filename)

    if conn is None:
        debug_literal("load_sqlite_file: conn is None.")
        return

    c = conn.cursor()

    # now, use sourmash_internal table to figure out what it can do.
    try:
        c.execute('SELECT DISTINCT key, value FROM sourmash_internal')
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        debug_literal("load_sqlite_file: no sourmash_internal table")
        return

    results = c.fetchall()

    is_index = False
    is_manifest = False
    for k, v in results:
        if k == 'SqliteIndex':
            # @CTB: check how we do version errors on sbt
            if v != '1.0':
                raise Exception(f"unknown SqliteIndex version '{v}'")
            is_index = True
            debug_literal("load_sqlite_file: it's an index!")
        elif k == 'SqliteManifest':
            if v != '1.0':
                raise Exception(f"unknown SqliteManifest version '{v}'")
            assert v == '1.0'
            is_manifest = True
            debug_literal("load_sqlite_file: it's a manifest!")
        # it's ok if there's no match, that just means we added keys
        # for some other type of sourmash SQLite database. #futureproofing.

    # every Index is a Manifest!
    if is_index:
        assert is_manifest

    idx = None
    if is_index and not request_manifest:
        conn.close()
        idx = SqliteIndex(filename)
        debug_literal("load_sqlite_file: returning SqliteIndex")
    elif is_manifest:
        managed_by_index=False
        if is_index:
            assert request_manifest
            managed_by_index=True

        prefix = os.path.dirname(filename)
        mf = SqliteCollectionManifest(conn, managed_by_index=managed_by_index)
        idx = StandaloneManifestIndex(mf, filename, prefix=prefix)
        debug_literal("load_sqlite_file: returning StandaloneManifestIndex")

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
        c.execute("SELECT DISTINCT scaled FROM sketches")
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
            # note: here we will want to re-open database.
            conn = sqlite3.connect(dbfile,
                                   detect_types=sqlite3.PARSE_DECLTYPES)

            c = conn.cursor()

            c.execute("PRAGMA cache_size=10000000")
            c.execute("PRAGMA synchronous = OFF")
            c.execute("PRAGMA journal_mode = MEMORY")
            c.execute("PRAGMA temp_store = MEMORY")

            if not empty_ok:
                c.execute("SELECT * FROM hashes LIMIT 1")
                c.fetchone()
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            raise ValueError(f"cannot open '{dbfile}' as SqliteIndex database")

        return conn

    @classmethod
    def create(cls, dbfile):
        conn = cls._open(dbfile, empty_ok=True)
        cls._create_tables(conn.cursor())
        conn.commit()

        return cls(dbfile, conn=conn)

    @classmethod
    def _create_tables(cls, c):
        "Create sqlite tables for SqliteIndex"
        # @CTB check what happens when you try to append to existing.
        try:
            sqlite_utils.add_sourmash_internal(c, 'SqliteIndex', '1.0')
            SqliteCollectionManifest._create_table(c)

            c.execute("""
            CREATE TABLE IF NOT EXISTS hashes (
               hashval INTEGER NOT NULL,
               sketch_id INTEGER NOT NULL,
               FOREIGN KEY (sketch_id) REFERENCES sketches (id)
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS hashval_idx ON hashes (
               hashval,
               sketch_id
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS hashval_idx2 ON hashes (
               hashval
            )
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS sketch_idx ON hashes (
               sketch_id
            )
            """
            )
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            raise ValueError(f"cannot create SqliteIndex tables")

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
        hashes = []
        hashes_to_sketch = []
        for h in ss.minhash.hashes:
            hh = convert_hash_to(h)
            hashes_to_sketch.append((hh, sketch_id))

        c.executemany("INSERT INTO hashes (hashval, sketch_id) VALUES (?, ?)",
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
        c2 = self.conn.cursor()

        # have the manifest run a select...
        picklist = self.manifest._run_select(c)

        #... and then operate on the results of that in 'c'
        for ss, loc, iloc in self._load_sketches(c, c2):
            if picklist is None or ss in picklist:
                yield ss, loc

    def save(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def load(self, dbfile):
        return SqliteIndex(dbfile)

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
            # CTB if we are doing containment only, we could break here.
            # but for Jaccard, we must continue.
            # see 'test_sqlite_jaccard_ordering'

            if search_fn.passes(score):
                subj = self._load_sketch(c2, sketch_id)

                # for testing only, I guess? remove this after validation :) @CTB
                subj_mh = subj.minhash.downsample(scaled=query_mh.scaled)
                int_size, un_size = query_mh.intersection_and_union_size(subj_mh)

                query_size = len(query_mh)
                subj_size = len(subj_mh)
                score2 = search_fn.score_fn(query_size, int_size, subj_size,
                                            un_size)
                assert score == score2

                if search_fn.collect(score, subj):
                    if picklist is None or subj in picklist:
                        yield IndexSearchResult(score, subj, self.location)

    def select(self, *, num=0, track_abundance=False, **kwargs):
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

        # return a new SqliteIndex with a new manifest, but same old conn.
        return SqliteIndex(self.dbfile,
                           sqlite_manifest=manifest,
                           conn=self.conn)

    #
    # Actual SQL queries, etc.
    #

    def _load_sketch_size(self, c1, sketch_id, max_hash):
        "Get sketch size for given sketch, downsampled by max_hash."
        if max_hash <= MAX_SQLITE_INT:
            c1.execute("""
            SELECT COUNT(hashval) FROM hashes
            WHERE sketch_id=? AND hashval >= 0 AND hashval <= ?""",
                       (sketch_id, max_hash))
        else:
            c1.execute('SELECT COUNT(hashval) FROM hashes WHERE sketch_id=?',
                       (sketch_id,))

        n_hashes, = c1.fetchone()
        return n_hashes

    def _load_sketch(self, c, sketch_id, *, match_scaled=None):
        "Load an individual sketch. If match_scaled is set, downsample."

        start = time.time()
        c.execute("""
        SELECT id, name, scaled, ksize, filename, moltype, seed
        FROM sketches WHERE id=?""",
                   (sketch_id,))
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
            hash_constraint_str = "hashes.hashval >= 0 AND hashes.hashval <= ? AND"
            template_values.insert(0, max_hash)
        else:
            debug_literal('NOT EMPLOYING hash_constraint_str')

        debug_literal(f"finding hashes for sketch {sketch_id} in {time.time() - start:.2f}")
        c.execute(f"SELECT hashval FROM hashes WHERE {hash_constraint_str} hashes.sketch_id=?", template_values)

        debug_literal(f"loading hashes for sketch {sketch_id} in {time.time() - start:.2f}")
        xy = c.fetchall()
        debug_literal(f"adding hashes for sketch {sketch_id} in {time.time() - start:.2f}")
        for hashval, in xy:
            hh = convert_hash_from(hashval)
            mh.add_hash(hh)

        debug_literal(f"done loading sketch {sketch_id} {time.time() - start:.2f})")

        ss = SourmashSignature(mh, name=name, filename=filename)
        return ss

    def _load_sketches(self, c1, c2):
        """Load sketches based on results from 'c1', using 'c2'.

        Here, 'c1' should already have run an appropriate 'select' on
        'sketches'. 'c2' will be used to load the hash values.
        """
        for sketch_id, name, num, scaled, ksize, filename, moltype, seed in c1:
            assert num == 0

            is_protein = 1 if moltype=='protein' else 0
            is_dayhoff = 1 if moltype=='dayhoff' else 0
            is_hp = 1 if moltype=='hp' else 0

            mh = MinHash(n=0, ksize=ksize, scaled=scaled, seed=seed,
                         is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp)
            c2.execute("SELECT hashval FROM hashes WHERE sketch_id=?",
                       (sketch_id,))

            hashvals = c2.fetchall()
            for hashval, in hashvals:
                mh.add_hash(convert_hash_from(hashval))

            ss = SourmashSignature(mh, name=name, filename=filename)
            yield ss, self.dbfile, sketch_id

    def _get_matching_sketches(self, c, hashes, max_hash):
        """
        For hashvals in 'hashes', retrieve all matching sketches,
        together with the number of overlapping hashes for each sketh.

        CTB: we do not use sqlite manifest conditions on this select,
        because it slows things down in practice.
        """
        c.execute("DROP TABLE IF EXISTS hash_query")
        c.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER PRIMARY KEY)")

        hashvals = [ (convert_hash_to(h),) for h in hashes ]
        c.executemany("INSERT OR IGNORE INTO hash_query (hashval) VALUES (?)", hashvals)

        #
        # set up SELECT conditions
        #

        conditions = []
        template_values = []

        # downsample? => add to conditions
        max_hash = min(max_hash, max(hashes))
        if max_hash <= MAX_SQLITE_INT:
            select_str = "hashes.hashval >= 0 AND hashes.hashval <= ?"
            conditions.append(select_str)
            template_values.append(max_hash)

        # format conditions
        conditions.append('hashes.hashval=hash_query.hashval')
        conditions = " AND ".join(conditions)

        c.execute(f"""
        SELECT DISTINCT hashes.sketch_id,COUNT(hashes.hashval) as CNT
        FROM hashes,hash_query
        WHERE {conditions}
        GROUP BY hashes.sketch_id ORDER BY CNT DESC
        """, template_values)

        return c


class SqliteCollectionManifest(BaseCollectionManifest):
    def __init__(self, conn, *, selection_dict=None, managed_by_index=False):
        """
        Here, 'conn' should already be connected and configured.
        """
        assert conn is not None
        self.conn = conn
        self.selection_dict = selection_dict
        self.managed_by_index = managed_by_index

    @classmethod
    def create(cls, filename):
        conn = sqlite3.connect(filename)
        cursor = conn.cursor()
        cls._create_table(cursor)
        return cls(conn)

    @classmethod
    def _create_table(cls, cursor):
        "Create the manifest table."

        cursor.execute("""
        CREATE TABLE IF NOT EXISTS sourmash_internal (
            key TEXT,
            value TEXT
        )
        """)

        cursor.execute("""
        INSERT INTO sourmash_internal (key, value)
        VALUES ('SqliteManifest', '1.0')
        """)

        cursor.execute("""
        CREATE TABLE sketches
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

    def _insert_row(self, cursor, row, *, call_is_from_index=False):
        # check - is this manifest managed by SqliteIndex?
        if self.managed_by_index and not call_is_from_index:
            raise Exception("must use SqliteIndex.insert to add to this manifest")

        cursor.execute("""
        INSERT OR IGNORE INTO sketches
          (name, num, scaled, ksize, filename, md5sum, moltype,
           seed, n_hashes, with_abundance, internal_location)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (row['name'],
         row['num'],
         row['scaled'],
         row['ksize'],
         row['filename'],
         row['md5'],
         row['moltype'],
         row.get('seed', 42),
         row['n_hashes'],
         row['with_abundance'],
         row['internal_location']))

    @classmethod
    def create_from_manifest(cls, dbfile, manifest, *, append=False):
        return cls._create_manifest_from_rows(manifest.rows, location=dbfile,
                                              append=append)

    def __bool__(self):
        return bool(len(self))

    def __eq__(self, other):
        return list(self.rows) == list(other.rows)

    def __len__(self):
        c = self.conn.cursor()
        conditions, values, picklist = self._select_signatures(c)
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        c.execute(f"SELECT COUNT(*) FROM sketches {conditions}", values)
        count, = c.fetchone()

        # @CTB do we need to pay attention to picklist here?
        # @CTB yes - we don't do prefix checking, we do 'like' in SQL.
        # e.g. check ident.
        #
        # we can generate manifest and use 'picklist.matches_manifest_row'
        # on rows...? basically is there a place where this will be
        # different / can we find it and test it :grin:
        # count = 0
        # for row in self.rows:
        #    if picklist.matches_manifest_row(row):
        #       count += 1
        return count

    def _select_signatures(self, c):
        """
        Given cursor 'c', build a set of SQL SELECT conditions
        and matching value tuple that can be used to select the
        right sketches from the database.

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
                conditions.append("sketches.ksize = ?")
                values.append(select_d['ksize'])
            if 'num' in select_d and select_d['num'] > 0:
                # @CTB check num
                assert 0
            if 'scaled' in select_d and select_d['scaled'] > 0:
                conditions.append("sketches.scaled > 0")
            if 'containment' in select_d and select_d['containment']:
                conditions.append("sketches.scaled > 0")
            if 'moltype' in select_d and select_d['moltype'] is not None:
                moltype = select_d['moltype']
                assert moltype in ('DNA', 'protein', 'dayhoff', 'hp'), moltype
                conditions.append(f"moltype = '{moltype}'")

            picklist = select_d.get('picklist')

            # support picklists!
            if picklist is not None:
                c.execute("DROP TABLE IF EXISTS pickset")
                c.execute("CREATE TABLE pickset (sketch_id INTEGER)")

                transform = picklist_transforms[picklist.coltype]
                sql_stmt = picklist_selects[picklist.coltype]

                vals = [ (transform(v),) for v in picklist.pickset ]
                c.executemany(sql_stmt, vals)

                if picklist.pickstyle == PickStyle.INCLUDE:
                    conditions.append("""
                    sketches.id IN (SELECT sketch_id FROM pickset)
                    """)
                elif picklist.pickstyle == PickStyle.EXCLUDE:
                    conditions.append("""
                    sketches.id NOT IN (SELECT sketch_id FROM pickset)
                    """)

        return conditions, values, picklist

    def select_to_manifest(self, **kwargs):
        # Pass along all the selection kwargs to a new instance
        if self.selection_dict:
            # combine selects...
            d = dict(self.selection_dict)
            for k, v in kwargs.items():
                if k in d:
                    if d[k] is not None and d[k] != v:
                        raise ValueError(f"incompatible select on '{k}'")
                d[k] = v
            kwargs = d

        new_mf = SqliteCollectionManifest(self.conn, selection_dict=kwargs)
        if 'picklist' in kwargs:
            picklist = kwargs['picklist']
            for row in new_mf.rows:
                does_match = picklist.matches_manifest_row(row)
                if not does_match:
                    assert 0

        return new_mf

    def _run_select(self, c):
        conditions, values, picklist = self._select_signatures(c)
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        c.execute(f"""
        SELECT id, name, num, scaled, ksize, filename, moltype, seed
        FROM sketches {conditions}""",
                  values)

        return picklist

    def _extract_manifest(self):
        """
        Generate a CollectionManifest dynamically from the SQL database.
        """
        manifest_list = []
        for row in self.rows:
            if '_id' in row:
                del row['_id']
            manifest_list.append(row)

        return CollectionManifest(manifest_list)

    @property
    def rows(self):
        c1 = self.conn.cursor()

        conditions, values, picklist = self._select_signatures(c1)
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        c1.execute(f"""
        SELECT id, name, md5sum, num, scaled, ksize, filename, moltype,
        seed, n_hashes, internal_location FROM sketches {conditions}
        """, values)

        manifest_list = []
        for (_id, name, md5sum, num, scaled, ksize, filename, moltype,
             seed, n_hashes, iloc) in c1:
            row = dict(num=num, scaled=scaled, name=name, filename=filename,
                       n_hashes=n_hashes, with_abundance=0, ksize=ksize,
                       md5=md5sum, internal_location=iloc,
                       moltype=moltype, md5short=md5sum[:8],
                       _id=_id)
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
        c1 = self.conn.cursor()

        conditions, values, picklist = self._select_signatures(c1)
        if conditions:
            conditions = conditions = "WHERE " + " AND ".join(conditions)
        else:
            conditions = ""

        c1.execute(f"""
        SELECT DISTINCT internal_location FROM sketches {conditions}
        """, values)

        return ( iloc for iloc, in c1 )

    def __contains__(self, ss):
        md5 = ss.md5sum()

        c = self.conn.cursor()
        c.execute('SELECT COUNT(*) FROM sketches WHERE md5sum=?', (md5,))
        val, = c.fetchone()
        return bool(val)

    def to_picklist(self):
        "Convert this manifest to a picklist."
        picklist = SignaturePicklist('md5')

        c = self.conn.cursor()
        c.execute('SELECT DISTINCT md5sum FROM sketches')
        pickset = set()
        pickset.update(( val for val, in c ))
        picklist.pickset = pickset

        return picklist

    @classmethod
    def create_manifest(cls, locations_iter, *, include_signature=False):
        """Create a manifest from an iterator that yields (ss, location)

        Stores signatures in manifest rows by default.

        Note: do NOT catch exceptions here, so this passes through load excs.
        Note: ignores 'include_signature'.
        """
        def rows_iter():
            for ss, location in locations_iter:
                row = cls.make_manifest_row(ss, location,
                                            include_signature=False)
                yield row

        return cls._create_manifest_from_rows(rows_iter())

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
            db = load_sqlite_file(location, request_manifest=True)
            mf = db.manifest

        cursor = mf.conn.cursor()

        for row in rows_iter:
            mf._insert_row(cursor, row)

        mf.conn.commit()
        return mf


class LCA_Database_SqliteWrapper:
    # @CTB: test via roundtrip ;).

    def __init__(self, filename):
        from sourmash.tax.tax_utils import LineageDB_Sqlite

        try:
            sqlite_idx = SqliteIndex.load(filename)
            lineage_db = LineageDB_Sqlite(sqlite_idx.conn)

            conn = sqlite_idx.conn
            c = conn.cursor()
            c.execute('SELECT DISTINCT key, value FROM sourmash_internal')
            d = dict(c)
            #print(d)
            # @CTB
        except sqlite3.OperationalError:
            raise ValueError(f"cannot open '{filename}' as sqlite database.")

        c.execute('SELECT DISTINCT ksize FROM sketches')
        ksizes = set(( ksize for ksize, in c ))
        assert len(ksizes) == 1
        self.ksize = next(iter(ksizes))
        notify(f"setting ksize to {self.ksize}")

        c.execute('SELECT DISTINCT moltype FROM sketches')
        moltypes = set(( moltype for moltype, in c ))
        assert len(moltypes) == 1
        self.moltype = next(iter(moltypes))
        notify(f"setting moltype to {self.moltype}")

        self.scaled = sqlite_idx.scaled

        assert isinstance(sqlite_idx, SqliteIndex)
        assert sqlite_idx.scaled
        self.sqlidx = sqlite_idx
        self.lineage_db = lineage_db

        ## the below is done once, but could be implemented as something
        ## ~dynamic.
        self._build_index()

    def _build_index(self):
        mf = self.sqlidx.manifest
        lineage_db = self.lineage_db

        ident_to_name = {}
        ident_to_idx = {}
        next_lid = 0
        idx_to_lid = {}
        lineage_to_lid = {}
        lid_to_lineage = {}

        for row in mf.rows:
            name = row['name']
            if name:
                ident = name.split(' ')[0].split('.')[0]

                assert ident not in ident_to_name
                ident_to_name[ident] = name
                idx = row['_id'] # this is only present in sqlite manifests.
                ident_to_idx[ident] = idx

                lineage = lineage_db.get(ident)
                if lineage:
                    lid = lineage_to_lid.get(lineage)
                    if lid is None:
                        lid = next_lid
                        next_lid += 1
                        lineage_to_lid[lineage] = lid
                        lid_to_lineage[lid] = lineage
                        idx_to_lid[idx] = lid

        self.ident_to_name = ident_to_name
        self.ident_to_idx = ident_to_idx
        self._next_lid = next_lid
        self.idx_to_lid = idx_to_lid
        self.lineage_to_lid = lineage_to_lid
        self.lid_to_lineage = lid_to_lineage

    def __len__(self):
        assert 0
        return len(self.sqlidx)

    def signatures(self):
        return self.sqlidx.signatures()

    def search(self, *args, **kwargs):
        return self.sqlidx.search(*args, **kwargs)

    def gather(self, *args, **kwargs):
        return self.sqlidx.gather(*args, **kwargs)

    def _get_ident_index(self, ident, fail_on_duplicate=False):
        "Get (create if nec) a unique int id, idx, for each identifier."
        assert 0

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
        assert 0
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

    def insert(self, *args, **kwargs):
        raise NotImplementedError

    def __repr__(self):
        return "LCA_Database_SqliteWrapper('{}')".format(self.sqlidx.location)

    def load(self, *args, **kwargs):
        # this could do the appropriate MultiLineageDB stuff.
        raise NotImplementedError

    def downsample_scaled(self, scaled):
        # @CTB this is necessary for internal implementation reasons,
        # but is not required technically.
        if scaled < self.sqlidx.scaled:
            assert 0
        else:
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

    @property
    def hashval_to_idx(self):
        "Dynamically interpret the SQL 'hashes' table like it's a dict."
        return _SqliteIndexHashvalToIndex(self.sqlidx)

    @property
    def conn(self):
        return self.sqlidx.conn

    @property
    def hashvals(self):
        return iter(_SqliteIndexHashvalToIndex(self.sqlidx))

    def get_identifiers_for_hashval(self, hashval):
        idxlist = self.hashval_to_idx[hashval]
        for idx in idxlist:
            yield self.idx_to_ident[idx]


class _SqliteIndexHashvalToIndex:
    def __init__(self, sqlidx):
        self.sqlidx = sqlidx

    def __iter__(self):
        c = self.sqlidx.conn.cursor()
        c.execute('SELECT DISTINCT hashval FROM hashes')
        for hashval, in c:
            yield hashval

    def items(self):
        "Retrieve hashval, idxlist for all hashvals."
        sqlidx = self.sqlidx
        c = sqlidx.cursor()

        c.execute('SELECT hashval, sketch_id FROM hashes ORDER BY hashval')

        this_hashval = None
        idxlist = []
        for hashval, sketch_id in c:
            if hashval == this_hashval:
                idxlist.append(sketch_id)
            else:
                if idxlist:
                    hh = convert_hash_from(this_hashval)
                    yield hh, idxlist

                this_hashval = hashval
                idxlist = [sketch_id]

        if idxlist:
            hh = convert_hash_from(this_hashval)
            yield hh, idxlist

    def get(self, key, dv=None):
        "Retrieve idxlist for a given hash."
        sqlidx = self.sqlidx
        c = sqlidx.cursor()

        hh = convert_hash_to(key)

        c.execute('SELECT sketch_id FROM hashes WHERE hashval=?', (hh,))

        x = [ convert_hash_from(h) for h, in c ]
        return x or dv

    def __getitem__(self, key):
        v = self.get(key)
        if v is None:
            raise KeyError(key)
        return v
