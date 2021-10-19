"""
Manifests for collections of signatures.
"""
import csv

from sourmash.picklist import SignaturePicklist


class CollectionManifest:
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

    def __init__(self, rows):
        "Initialize from an iterable of metadata dictionaries."
        self.rows = tuple(rows)

        # build a fast lookup table for md5sums in particular
        md5set = set()
        for row in self.rows:
            md5set.add(row['md5'])
        self._md5_set = md5set

    def __bool__(self):
        return bool(self.rows)

    def __len__(self):
        return len(self.rows)

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
                row[k] = bool(row[k])
            row['signature'] = None
            manifest_list.append(row)

        return cls(manifest_list)

    @classmethod
    def write_csv_header(cls, fp):
        "write header for manifest CSV format"
        fp.write('# SOURMASH-MANIFEST-VERSION: 1.0\n')
        w = csv.DictWriter(fp, fieldnames=cls.required_keys)
        w.writeheader()

    def write_to_csv(self, fp, write_header=False):
        "write manifest CSV to specified file handle"
        w = csv.DictWriter(fp, fieldnames=self.required_keys)

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
            row = cls.make_manifest_row(ss, location, include_signature=True)
            manifest_list.append(row)

        return cls(manifest_list)

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
