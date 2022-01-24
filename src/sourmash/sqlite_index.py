import sqlite3

from .index import Index
import sourmash
from sourmash import MinHash, SourmashSignature
from sourmash.index import IndexSearchResult
from collections import Counter


MAX_SQLITE_INT = 2 ** 63 - 1
sqlite3.register_adapter(
    int, lambda x: hex(x) if x > MAX_SQLITE_INT else x)
sqlite3.register_converter(
    'integer', lambda b: int(b, 16 if b[:2] == b'0x' else 10))

def load_sketch(db, sketch_id):
    c2 = db.cursor()

    c2.execute("SELECT name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed FROM sketches WHERE id=?", (sketch_id,))

    name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed = c2.fetchone()

    mh = sourmash.MinHash(n=num, ksize=ksize, scaled=scaled, seed=seed, is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp, track_abundance=track_abundance)

    c2.execute("SELECT hashval FROM hashes WHERE sketch_id=?", (sketch_id,))

    for hashval, in c2:
        mh.add_hash(hashval)

    ss = sourmash.SourmashSignature(mh, name=name, filename=filename)
    return ss


def get_matching_sketches(db, unitig_mh):
    query_cursor = db.cursor()
    query_cursor.execute("DROP TABLE IF EXISTS hash_query")
    query_cursor.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER PRIMARY KEY)")
    for hashval in unitig_mh.hashes:
        query_cursor.execute("INSERT INTO hash_query (hashval) VALUES (?)", (hashval,))

    # do we have an overlap with any query at all??
    query_cursor.execute("SELECT DISTINCT sketches.id FROM sketches,hashes WHERE sketches.id=hashes.sketch_id AND hashes.hashval IN (SELECT hashval FROM hash_query)")

    for sketch_id, in query_cursor:
        yield load_sketch(db, sketch_id)


def get_matching_hashes(query_cursor, unitig_mh):
    query_cursor.execute("DROP TABLE IF EXISTS hash_query")
    query_cursor.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER PRIMARY KEY)")
    for hashval in unitig_mh.hashes:
        query_cursor.execute("INSERT INTO hash_query (hashval) VALUES (?)", (hashval,))

    query_cursor.execute("SELECT DISTINCT hashes.sketch_id,hashes.hashval FROM hashes,hash_query WHERE hashes.hashval=hash_query.hashval")

    for sketch_id, hashval in query_cursor:
        yield sketch_id, hashval


class SqliteIndex(Index):
    is_database = True
    
    def __init__(self, dbfile):
        self.dbfile = dbfile
        self.conn = sqlite3.connect(dbfile,
                                    detect_types=sqlite3.PARSE_DECLTYPES)

        c = self.conn.cursor()

        c.execute("PRAGMA cache_size=1000000")
        c.execute("PRAGMA synchronous = OFF")
        c.execute("PRAGMA journal_mode = MEMORY")

        c.execute("CREATE TABLE IF NOT EXISTS sketches (id INTEGER PRIMARY KEY, name TEXT, num INTEGER NOT NULL, scaled INTEGER NOT NULL, ksize INTEGER NOT NULL, filename TEXT, is_dna BOOLEAN, is_protein BOOLEAN, is_dayhoff BOOLEAN, is_hp BOOLEAN, track_abundance BOOLEAN, seed INTEGER NOT NULL)")
        c.execute("CREATE TABLE IF NOT EXISTS hashes (hashval INTEGER NOT NULL, sketch_id INTEGER NOT NULL, FOREIGN KEY (sketch_id) REFERENCES sketches (id))")
    def close(self):
        self.conn.close()

    def commit(self):
        self.conn.commit()

    def insert(self, ss, commit=True):
        c = self.conn.cursor()
        c.execute("INSERT INTO sketches (name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (ss.name, ss.minhash.num, ss.minhash.scaled, ss.minhash.ksize, ss.filename, ss.minhash.is_dna, ss.minhash.is_protein, ss.minhash.dayhoff, ss.minhash.hp, ss.minhash.track_abundance, ss.minhash.seed))
        c.execute("SELECT last_insert_rowid()")
        id, = c.fetchone()
        for h in ss.minhash.hashes:
            c.execute("INSERT INTO hashes (hashval, sketch_id) VALUES (?, ?)", (h, id))

        if commit:
            self.conn.commit()

    @property
    def location(self):
        return self.dbfile

    def signatures(self):
        "Return an iterator over all signatures in the Index object."
        for ss, loc, iloc in self._signatures_with_internal():
            yield ss

    def signatures_with_location(self):
        "Return an iterator over tuples (signature, location) in the Index."
        for ss, loc, iloc in self._signatures_with_internal():
            yield ss, loc

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, location, internal_location).

        This is an internal API for use in generating manifests, and may
        change without warning.

        This method should be implemented separately for each Index object.
        """
        c = self.conn.cursor()
        c2 = self.conn.cursor()

        c.execute("SELECT id, name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed FROM sketches")
        for (sketch_id, name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed) in c:
            mh = MinHash(n=num, ksize=ksize, scaled=scaled, seed=seed, is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp, track_abundance=track_abundance)
            c2.execute("SELECT hashval FROM hashes WHERE sketch_id=?", (sketch_id,))

            for hashval, in c2:
                mh.add_hash(hashval)

            ss = SourmashSignature(mh, name=name, filename=filename)
            yield ss, self.dbfile, sketch_id

    def save(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def load(self, dbfile):
        return SqliteIndex(dbfile)

    def find(self, search_fn, query, **kwargs):
        search_fn.check_is_compatible(query)

        # check compatibility, etc. @CTB
        query_mh = query.minhash

        cursor = self.conn.cursor()
        c = Counter()
        for sketch_id, hashval in get_matching_hashes(cursor, query_mh):
            c[sketch_id] += 1

        for sketch_id, count in c.most_common():
            subj = load_sketch(self.conn, sketch_id)

            # @CTB more goes here

            subj_mh = subj.minhash

            # all numbers calculated after downsampling --
            query_size = len(query_mh)
            subj_size = len(subj_mh)
            shared_size = query_mh.count_common(subj_mh)
            total_size = len(query_mh + subj_mh)

            score = search_fn.score_fn(query_size, shared_size, subj_size,
                                       total_size)

            if search_fn.passes(score):
                if search_fn.collect(score, subj):
                    if 1: #passes_all_picklists(subj, self.picklists):
                        yield IndexSearchResult(score, subj, self.location)

    def select(self, ksize=None, moltype=None, scaled=None, num=None,
               abund=None, containment=None, picklist=None):
        return self
