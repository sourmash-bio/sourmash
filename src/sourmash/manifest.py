"""
Manifests for collections of signatures.
"""
import csv
import ast
import gzip
import os.path
from abc import abstractmethod
import itertools

from sourmash.picklist import SignaturePicklist


class BaseCollectionManifest:
    """
    Signature metadata for a collection of signatures.

    Manifests support selection and rapid lookup of signatures.

    * 'select_to_manifest(...)' matches the Index selector protocol
    * 'rows' is a public iterable that can be used to iterate over the manifest
       contents.
    * 'locations()' returns all distinct locations for e.g. lazy loading
    * supports container protocol for signatures, e.g. 'if ss in manifest: ...'
    """
    # each manifest row must have the following, although they may be empty.
    required_keys = ('internal_location',
                     'md5', 'md5short', 'ksize', 'moltype', 'num',
                     'scaled', 'n_hashes', 'with_abundance',
                     'name', 'filename')

    @classmethod
    @abstractmethod
    def load_from_manifest(cls, manifest, **kwargs):
        "Load this manifest from another manifest object."

    @classmethod
    def load_from_filename(cls, filename):
        # SQLite db?
        db = cls.load_from_sql(filename)
        if db is not None:
            return db

        # not a SQLite db? CTB: fix this to actually try loading this as .gz...
        if filename.endswith('.gz'):
            xopen = gzip.open
        else:
            xopen = open

        with xopen(filename, 'rt', newline="") as fp:
            return cls.load_from_csv(fp)

    @classmethod
    def load_from_csv(cls, fp):
        "load a manifest from a CSV file."
        manifest_list = []
        firstline = fp.readline().rstrip()
        if not firstline.startswith('# SOURMASH-MANIFEST-VERSION: '):
            raise ValueError("manifest is missing version header")

        version = firstline[len('# SOURMASH-MANIFEST-VERSION: '):]
        if float(version) != 1.0:
            raise ValueError(f"unknown manifest version number {version}")

        r = csv.DictReader(fp)
        if not r.fieldnames:
            raise ValueError("missing column headers in manifest")

        for k in cls.required_keys:
            if k not in r.fieldnames:
                raise ValueError(f"missing column '{k}' in manifest.")

        row = None

        # do row type conversion
        introws = ('num', 'scaled', 'ksize', 'n_hashes')
        boolrows = ('with_abundance',)

        for row in r:
            for k in introws:
                row[k] = int(row[k])
            for k in boolrows:
                row[k] = bool(ast.literal_eval(str(row[k])))
            row['signature'] = None
            manifest_list.append(row)

        return CollectionManifest(manifest_list)

    @classmethod
    def load_from_sql(cls, filename):
        from sourmash.index.sqlite_index import load_sqlite_index
        db = load_sqlite_index(filename, request_manifest=True)
        if db is not None:
            return db.manifest

        return None

    def write_to_filename(self, filename, *, database_format='csv',
                          ok_if_exists=False):
        if database_format == 'csv':
            from .sourmash_args import FileOutputCSV
            if ok_if_exists or not os.path.exists(filename):
                with FileOutputCSV(filename) as fp:
                    return self.write_to_csv(fp, write_header=True)
            elif os.path.exists(filename) and not ok_if_exists:
                raise Exception("output manifest already exists")

        elif database_format == 'sql':
            from sourmash.index.sqlite_index import SqliteCollectionManifest
            SqliteCollectionManifest.load_from_manifest(self, dbfile=filename,
                                                        append=ok_if_exists)

    @classmethod
    def write_csv_header(cls, fp):
        "write header for manifest CSV format"
        fp.write('# SOURMASH-MANIFEST-VERSION: 1.0\n')
        w = csv.DictWriter(fp, fieldnames=cls.required_keys)
        w.writeheader()

    def write_to_csv(self, fp, write_header=False):
        "write manifest CSV to specified file handle"
        w = csv.DictWriter(fp, fieldnames=self.required_keys,
                           extrasaction='ignore')

        if write_header:
            self.write_csv_header(fp)

        for row in self.rows:
            # don't write signature!
            if 'signature' in row:
                del row['signature']
            w.writerow(row)

    @classmethod
    def make_manifest_row(cls, ss, location, *, include_signature=True):
        "make a manifest row dictionary."
        row = {}
        row['md5'] = ss.md5sum()
        row['md5short'] = row['md5'][:8]
        row['ksize'] = ss.minhash.ksize
        row['moltype'] = ss.minhash.moltype
        row['num'] = ss.minhash.num
        row['scaled'] = ss.minhash.scaled
        row['n_hashes'] = len(ss.minhash)
        row['with_abundance'] = 1 if ss.minhash.track_abundance else 0
        row['name'] = ss.name
        row['filename'] = ss.filename
        row['internal_location'] = location

        assert set(row.keys()) == set(cls.required_keys)

        # if requested, include the signature in the manifest.
        if include_signature:
            row['signature'] = ss
        return row

    @classmethod
    def create_manifest(cls, locations_iter, *, include_signature=True):
        """Create a manifest from an iterator that yields (ss, location)

        Stores signatures in manifest rows by default.

        Note: do NOT catch exceptions here, so this passes through load excs.
        """
        manifest_list = []
        for ss, location in locations_iter:
            row = cls.make_manifest_row(ss, location,
                                        include_signature=include_signature)
            manifest_list.append(row)

        return cls(manifest_list)

    ## implement me
    @abstractmethod
    def __add__(self, other):
        "Add two manifests"

    @abstractmethod
    def __bool__(self):
        "Test if manifest is empty"

    @abstractmethod
    def __len__(self):
        "Get number of entries in manifest"

    @abstractmethod
    def __eq__(self, other):
        "Check for equality of manifest based on rows"

    @abstractmethod
    def select_to_manifest(self, **kwargs):
        "Select compatible signatures"

    @abstractmethod
    def filter_rows(self, row_filter_fn):
        "Filter rows based on a pattern matching function."

    @abstractmethod
    def filter_on_columns(self, col_filter_fn, col_names):
        "Filter on column values."

    @abstractmethod
    def locations(self):
        "Return a list of distinct locations"

    @abstractmethod
    def __contains__(self, ss):
        "Determine if a particular SourmashSignature is in this manifest."

    @abstractmethod
    def to_picklist(self):
        "Convert manifest to a picklist."


class CollectionManifest(BaseCollectionManifest):
    """
    An in-memory manifest that simply stores the rows in a list.
    """
    def __init__(self, rows=[]):
        "Initialize from an iterable of metadata dictionaries."
        self.rows = []
        self._md5_set = set()

        self._add_rows(rows)

    @classmethod
    def load_from_manifest(cls, manifest, **kwargs):
        "Load this manifest from another manifest object."
        return cls(manifest.rows)

    def add_row(self, row):
        self._add_rows([row])

    def _add_rows(self, rows):
        self.rows.extend(rows)

        # maintain a fast check for md5sums for __contains__ check.
        md5set = self._md5_set
        for row in self.rows:
            md5set.add(row['md5'])

    def __iadd__(self, other):
        self._add_rows(other.rows)
        return self

    def __add__(self, other):
        mf = CollectionManifest(self.rows)
        mf._add_rows(other.rows)
        return mf

    def __bool__(self):
        return bool(self.rows)

    def __len__(self):
        return len(self.rows)

    def __eq__(self, other):
        "Check equality on a row-by-row basis. May fail on out-of-order rows."
        for (a, b) in itertools.zip_longest(self.rows, other.rows):
            if a is None or b is None:
                return False

            # ignore non-required keys.
            for k in self.required_keys:
                if a[k] != b[k]:
                    return False

        return True

    def _select(self, *, ksize=None, moltype=None, scaled=0, num=0,
                containment=False, abund=None, picklist=None):
        """Yield manifest rows for sigs that match the specified requirements.

        Internal method; call `select_to_manifest` instead.
        """
        matching_rows = self.rows
        if ksize:
            matching_rows = ( row for row in matching_rows
                              if row['ksize'] == ksize )
        if moltype:
            matching_rows = ( row for row in matching_rows
                              if row['moltype'] == moltype )
        if scaled or containment:
            if containment and not scaled:
                raise ValueError("'containment' requires 'scaled' in Index.select'")

            matching_rows = ( row for row in matching_rows
                              if row['scaled'] and not row['num'] )
        if num:
            matching_rows = ( row for row in matching_rows
                              if row['num'] and not row['scaled'] )

        if abund:
            # only need to concern ourselves if abundance is _required_
            matching_rows = ( row for row in matching_rows
                              if row['with_abundance'] )

        if picklist:
            matching_rows = ( row for row in matching_rows
                              if picklist.matches_manifest_row(row) )

        # return only the internal filenames!
        for row in matching_rows:
            yield row

    def select_to_manifest(self, **kwargs):
        "Do a 'select' and return a new CollectionManifest object."
        new_rows = self._select(**kwargs)
        return CollectionManifest(new_rows)

    def filter_rows(self, row_filter_fn):
        "Create a new manifest filtered through row_filter_fn."
        new_rows = [ row for row in self.rows if row_filter_fn(row) ]

        return CollectionManifest(new_rows)

    def filter_on_columns(self, col_filter_fn, col_names):
        "Create a new manifest based on column matches."
        def row_filter_fn(row):
            x = [ row[col] for col in col_names if row[col] is not None ]
            return col_filter_fn(x)
        return self.filter_rows(row_filter_fn)

    def locations(self):
        "Return all distinct locations."
        seen = set()
        for row in self.rows:
            loc = row['internal_location']

            # track/remove duplicates
            if loc not in seen:
                seen.add(loc)
                yield loc

    def __contains__(self, ss):
        "Does this manifest contain this signature?"
        md5 = ss.md5sum()
        return md5 in self._md5_set

    def to_picklist(self):
        "Convert this manifest to a picklist."
        picklist = SignaturePicklist('md5')
        picklist.pickset = set(self._md5_set)

        return picklist
