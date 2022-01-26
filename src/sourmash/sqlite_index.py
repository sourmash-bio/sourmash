import sqlite3

# add DISTINCT to sketch and hash select

from .index import Index
import sourmash
from sourmash import MinHash, SourmashSignature
from sourmash.index import IndexSearchResult
from collections import Counter

# register converters for unsigned 64-bit ints. @CTB stackoverflow link.
MAX_SQLITE_INT = 2 ** 63 - 1
sqlite3.register_adapter(
    int, lambda x: hex(x) if x > MAX_SQLITE_INT else x)
sqlite3.register_converter(
    'integer', lambda b: int(b, 16 if b[:2] == b'0x' else 10))


class SqliteIndex(Index):
    is_database = True
    
    def __init__(self, dbfile, selection_dict=None):
        self.dbfile = dbfile
        self.selection_dict = selection_dict
        try:
            self.conn = sqlite3.connect(dbfile,
                                        detect_types=sqlite3.PARSE_DECLTYPES)

            c = self.conn.cursor()

            c.execute("PRAGMA cache_size=1000000")
            c.execute("PRAGMA synchronous = OFF")
            c.execute("PRAGMA journal_mode = MEMORY")
            c.execute("PRAGMA temp_store = MEMORY")

            c.execute("""
            CREATE TABLE IF NOT EXISTS sketches
              (id INTEGER PRIMARY KEY,
               name TEXT, num INTEGER NOT NULL,
               scaled INTEGER NOT NULL,
               ksize INTEGER NOT NULL,
               filename TEXT,
               is_dna BOOLEAN NOT NULL,
               is_protein BOOLEAN NOT NULL,
               is_dayhoff BOOLEAN NOT NULL,
               is_hp BOOLEAN NOT NULL,
               track_abundance BOOLEAN NOT NULL,
               md5sum TEXT NOT NULL, 
               seed INTEGER NOT NULL)
            """)
            c.execute("""
            CREATE TABLE IF NOT EXISTS hashes
              (hashval INTEGER NOT NULL,
               sketch_id INTEGER NOT NULL,
               FOREIGN KEY (sketch_id) REFERENCES sketches (id))
            """)
            c.execute("""
            CREATE INDEX IF NOT EXISTS hashval_idx ON hashes (hashval)
            """)

        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            raise ValueError(f"cannot open '{dbfile}' as sqlite3 database")

    def cursor(self):
        return self.conn.cursor()

    def close(self):
        self.conn.close()

    def commit(self):
        self.conn.commit()

    def insert(self, ss, cursor=None, commit=True):
        if cursor:
            c = cursor
        else:
            c = self.conn.cursor()

        c.execute("""
        INSERT INTO sketches
          (name, num, scaled, ksize, filename, md5sum,
           is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (ss.name, ss.minhash.num, ss.minhash.scaled, ss.minhash.ksize,
         ss.filename, ss.md5sum(),
         ss.minhash.is_dna, ss.minhash.is_protein, ss.minhash.dayhoff,
         ss.minhash.hp, ss.minhash.track_abundance, ss.minhash.seed))

        c.execute("SELECT last_insert_rowid()")
        sketch_id, = c.fetchone()

        hashes = []
        for h in ss.minhash.hashes:
            hashes.append((h, sketch_id))

        c.executemany("INSERT INTO hashes (hashval, sketch_id) VALUES (?, ?)",
                      hashes)

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

        conditions = []
        values = []
        picklist = None
        if self.selection_dict:
            select_d = self.selection_dict
            if 'ksize' in select_d and select_d['ksize']:
                conditions.append("sketches.ksize = ?")
                values.append(select_d['ksize'])
            if 'scaled' in select_d and select_d['scaled'] > 0:
                conditions.append("sketches.scaled > 0")
            if 'containment' in select_d and select_d['containment']:
                conditions.append("sketches.scaled > 0")
            if 'moltype' in select_d:
                moltype = select_d['moltype']
                if moltype == 'DNA':
                    conditions.append("sketches.is_dna")
                elif moltype == 'protein':
                    conditions.append("sketches.is_protein")
                elif moltype == 'dayhoff':
                    conditions.append("sketches.is_dayhoff")
                elif moltype == 'hp':
                    conditions.append("sketches.is_hp")
            # TODO: num, abund
            picklist = select_d.get('picklist')

            # note: support md5, md5short, name, ident?, identprefix?
            c.execute("DROP TABLE IF EXISTS pickset")
            c.execute("CREATE TABLE pickset (sketch_id INTEGER)")

            if picklist.coltype == 'name':
                names = [ (name,) for name in picklist.pickset ]
                c.executemany('INSERT INTO pickset SELECT id FROM sketches WHERE name=?', names)
                conditions.append('sketches.id in (SELECT sketch_id FROM pickset)')
            elif picklist.coltype == 'ident':
                pidents = [ (p + ' %',) for p in picklist.pickset ]
                c.executemany('INSERT INTO pickset SELECT id FROM sketches WHERE name LIKE ?', pidents)
                conditions.append('sketches.id in (SELECT sketch_id FROM pickset)')
            elif picklist.coltype == 'identprefix':
                pidents = [ (p + '%',) for p in picklist.pickset ]
                c.executemany('INSERT INTO pickset SELECT id FROM sketches WHERE name LIKE ?', pidents)
                conditions.append('sketches.id in (SELECT sketch_id FROM pickset)')
            elif picklist.coltype == 'md5short': # @CTB no md5sum in our table!
                md5shorts = [ (p[:8] + '%',) for p in picklist.pickset ]
                c.executemany('INSERT INTO pickset SELECT id FROM sketches WHERE md5sum LIKE ?', md5shorts)
                conditions.append('sketches.id in (SELECT sketch_id FROM pickset)')
            elif picklist.coltype == 'md5': # @CTB no md5sum in our table!
                md5s = [ (p,) for p in picklist.pickset ]
                c.executemany('INSERT INTO pickset SELECT id FROM sketches WHERE md5sum=?', md5s)
                conditions.append('sketches.id in (SELECT sketch_id FROM pickset)')

        if conditions:
            conditions = "WHERE " + " AND ".join(conditions)

        c.execute(f"""
        SELECT id, name, num, scaled, ksize, filename, is_dna, is_protein,
        is_dayhoff, is_hp, track_abundance, seed FROM sketches {conditions}""",
                  values)

        for ss, loc, iloc in self._load_sketches(c, c2):
            if picklist is None or ss in picklist:
                yield ss, loc

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, location, internal_location).

        This is an internal API for use in generating manifests, and may
        change without warning.

        This method should be implemented separately for each Index object.
        """
        c = self.conn.cursor()
        c2 = self.conn.cursor()

        c.execute("""
        SELECT id, name, num, scaled, ksize, filename, is_dna, is_protein,
        is_dayhoff, is_hp, track_abundance, seed FROM sketches
        """)
        
        for ss, loc, iloc in self._load_sketches(c, c2):
            yield ss, loc, iloc

    def _load_sketch(self, c1, sketch_id):
        # here, c1 should already have run an appropriate 'select' on 'sketches'
        # c2 will be used to load the hash values.
        c1.execute("""
        SELECT id, name, num, scaled, ksize, filename, is_dna, is_protein,
        is_dayhoff, is_hp, track_abundance, seed FROM sketches WHERE id=?""",
                   (sketch_id,))

        (sketch_id, name, num, scaled, ksize, filename, is_dna,
         is_protein, is_dayhoff, is_hp, track_abundance, seed) = c1.fetchone()
        mh = MinHash(n=num, ksize=ksize, scaled=scaled, seed=seed,
                     is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp,
                     track_abundance=track_abundance)

        c1.execute("SELECT hashval FROM hashes WHERE sketch_id=?", (sketch_id,))

        for hashval, in c1:
            mh.add_hash(hashval)

        ss = SourmashSignature(mh, name=name, filename=filename)
        return ss

    def _load_sketches(self, c1, c2):
        # here, c1 should already have run an appropriate 'select' on 'sketches'
        # c2 will be used to load the hash values.
        for (sketch_id, name, num, scaled, ksize, filename, is_dna, is_protein,
             is_dayhoff, is_hp, track_abundance, seed) in c1:
            mh = MinHash(n=num, ksize=ksize, scaled=scaled, seed=seed,
                         is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp,
                         track_abundance=track_abundance)
            c2.execute("SELECT hashval FROM hashes WHERE sketch_id=?",
                       (sketch_id,))

            hashvals = c2.fetchall()
            for hashval, in hashvals:
                mh.add_hash(hashval)

            ss = SourmashSignature(mh, name=name, filename=filename)
            yield ss, self.dbfile, sketch_id

    def _get_matching_hashes(self, c, hashes):
        c.execute("DROP TABLE IF EXISTS hash_query")
        c.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER)")

        hashvals = [ (h,) for h in hashes ]
        c.executemany("INSERT INTO hash_query (hashval) VALUES (?)", hashvals)

        # @CTB do we want to add select stuff on here?
        c.execute("""
        SELECT DISTINCT hashes.sketch_id,hashes.hashval FROM
        hashes,hash_query WHERE hashes.hashval=hash_query.hashval""")

        return c.fetchall()

    def save(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def load(self, dbfile):
        return SqliteIndex(dbfile)

    def find(self, search_fn, query, **kwargs):
        search_fn.check_is_compatible(query)

        # check compatibility, etc. @CTB
        query_mh = query.minhash

        picklist = self.selection_dict.get('picklist')

        cursor = self.conn.cursor()
        c = Counter()
        # @CTB do select here
        for sketch_id, hashval in self._get_matching_hashes(cursor, query_mh.hashes):
            c[sketch_id] += 1

        for sketch_id, count in c.most_common():
            subj = self._load_sketch(cursor, sketch_id)

            # @CTB more goes here?

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
                    if picklist is None or subj in picklist:
                        yield IndexSearchResult(score, subj, self.location)

    def select(self, **kwargs):
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

        return SqliteIndex(self.dbfile, selection_dict=kwargs)
