"""
Index object/sigfile loading and signature saving code.

This is the middleware code responsible for loading and saving signatures
in a variety of ways.

---

Command-line functionality goes in sourmash_args.py.

Low-level JSON reading/writing is in signature.py.

Index objects are implemented in the index submodule.

Public API:

* load_file_as_index(filename, ...) -- load a sourmash.Index class
* SaveSignaturesToLocation(filename) - bulk signature output

APIs for plugins to use:

* class Base_SaveSignaturesToLocation - to implement a new output method.

CTB TODO:
* consider replacing ValueError with IndexNotLoaded in the future.
"""
import sys
import os
import gzip
from io import StringIO
import zipfile
import itertools
import traceback

import screed
import sourmash

from . import plugins as sourmash_plugins
from .logging import notify, debug_literal
from .exceptions import IndexNotLoaded

from .index.sqlite_index import load_sqlite_index, SqliteIndex
from .sbtmh import load_sbt_index
from .lca.lca_db import load_single_database
from . import signature as sigmod
from .index import (LinearIndex, ZipFileLinearIndex, MultiIndex)
from .manifest import CollectionManifest


def load_file_as_index(filename, *, yield_all_files=False):
    """Load 'filename' as a database; generic database loader.

    If 'filename' contains an SBT or LCA indexed database, or a regular
    Zip file, will return the appropriate objects. If a Zip file and
    yield_all_files=True, will try to load all files within zip, not just
    .sig files.

    If 'filename' is a JSON file containing one or more signatures, will
    return an Index object containing those signatures.

    If 'filename' is a directory, will load *.sig underneath
    this directory into an Index object. If yield_all_files=True, will
    attempt to load all files.
    """
    return _load_database(filename, yield_all_files)


def SaveSignaturesToLocation(location):
    """
    Provides a context manager that saves signatures in various output formats.

    Usage:

    with SaveSignaturesToLocation(filename_or_location) as save_sigs:
       save_sigs.add(sig_obj)
    """
    save_list = itertools.chain(_save_classes,
                                sourmash_plugins.get_save_to_functions())
    for priority, cls in sorted(save_list, key=lambda x:x[0]):
        debug_literal(f"trying to match save function {cls}, priority={priority}")

        if cls.matches(location):
            debug_literal(f"{cls} is a match!")
            return cls(location)

    raise Exception(f"cannot determine how to open location {location} for saving; this should never happen!?")

### Implementation machinery for _load_databases


def _load_database(filename, traverse_yield_all, *, cache_size=None):
    """Load file as a database - list of signatures, LCA, SBT, etc.

    Return Index object.

    This is an internal function used by other functions in sourmash_args.
    """
    loaded = False

    # load plugins
    plugin_fns = sourmash_plugins.get_load_from_functions()

    # aggregate with default load_from functions & sort by priority
    load_from_functions = sorted(itertools.chain(_loader_functions,
                                                 plugin_fns))
                                                 
    # iterate through loader functions, sorted by priority; try them all.
    # Catch ValueError & IndexNotLoaded but nothing else.
    for (priority, desc, load_fn) in load_from_functions:
        db = None
        try:
            debug_literal(f"_load_databases: trying loader fn - priority {priority} - '{desc}'")
            db = load_fn(filename,
                         traverse_yield_all=traverse_yield_all,
                         cache_size=cache_size)
        except (ValueError, IndexNotLoaded):
            debug_literal(f"_load_databases: FAIL with ValueError: on fn {desc}.")
            debug_literal(traceback.format_exc())
            debug_literal("(continuing past exception)")

        if db is not None:
            loaded = True
            debug_literal("_load_databases: success!")
            break

    if loaded:
        assert db is not None
        return db
    
    raise ValueError(f"Error while reading signatures from '{filename}'.")


_loader_functions = []
def add_loader(name, priority):
    "decorator to add name/priority to _loader_functions"
    def dec_priority(func):
        _loader_functions.append((priority, name, func))
        return func
    return dec_priority


@add_loader("load from stdin", 10)
def _load_stdin(filename, **kwargs):
    "Load collection from .sig file streamed in via stdin"
    db = None
    if filename == '-':
        # load as LinearIndex, then pass into MultiIndex to generate a
        # manifest.
        lidx = LinearIndex.load(sys.stdin, filename='-')
        db = MultiIndex.load((lidx,), (None,), parent="-")

    return db


@add_loader("load from standalone manifest", 30)
def _load_standalone_manifest(filename, **kwargs):
    from sourmash.index import StandaloneManifestIndex

    try:
        idx = StandaloneManifestIndex.load(filename)
    except gzip.BadGzipFile as exc:
        raise IndexNotLoaded(exc)

    return idx


@add_loader("load from list of paths", 50)
def _multiindex_load_from_pathlist(filename, **kwargs):
    "Load collection from a list of signature/database files"
    db = MultiIndex.load_from_pathlist(filename)

    return db


@add_loader("load from path (file or directory)", 40)
def _multiindex_load_from_path(filename, **kwargs):
    "Load collection from a directory."
    traverse_yield_all = kwargs['traverse_yield_all']
    db = MultiIndex.load_from_path(filename, traverse_yield_all)

    return db


@add_loader("load SBT", 60)
def _load_sbt(filename, **kwargs):
    "Load collection from an SBT."
    cache_size = kwargs.get('cache_size')

    try:
        db = load_sbt_index(filename, cache_size=cache_size)
    except (FileNotFoundError, TypeError) as exc:
        raise IndexNotLoaded(exc)

    return db


@add_loader("load revindex", 70)
def _load_revindex(filename, **kwargs):
    "Load collection from an LCA database/reverse index."
    db, _, _ = load_single_database(filename)
    return db


@add_loader("load collection from sqlitedb", 20)
def _load_sqlite_db(filename, **kwargs):
    return load_sqlite_index(filename)


@add_loader("load collection from zipfile", 80)
def _load_zipfile(filename, **kwargs):
    "Load collection from a .zip file."
    db = None
    if filename.endswith('.zip'):
        traverse_yield_all = kwargs['traverse_yield_all']
        try:
            db = ZipFileLinearIndex.load(filename,
                                         traverse_yield_all=traverse_yield_all)
        except FileNotFoundError as exc:
            # turn this into an IndexNotLoaded => proper exception handling by
            # _load_database.
            raise IndexNotLoaded(exc)

    return db


@add_loader("catch FASTA/FASTQ files and error", 1000)
def _error_on_fastaq(filename, **kwargs):
    "This is a tail-end loader that checks for FASTA/FASTQ sequences => err."
    success = False
    try:
        with screed.open(filename) as it:
            _ = next(iter(it))

            success = True
    except:
        pass

    if success:
        raise Exception(f"Error while reading signatures from '{filename}' - got sequences instead! Is this a FASTA/FASTQ file?")


### Implementation machinery for SaveSignaturesToLocation

class Base_SaveSignaturesToLocation:
    "Base signature saving class. Track location (if any) and count."
    def __init__(self, location):
        self.location = location
        self.count = 0

    @classmethod
    def matches(cls, location):
        "returns True when this class should handle a specific location"
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

    def __len__(self):
        return self.count

    def open(self):
        pass

    def close(self):
        pass

    def __enter__(self):
        "provide context manager functionality"
        self.open()
        return self

    def __exit__(self, type, value, traceback):
        "provide context manager functionality"
        self.close()

    def add(self, ss):
        self.count += 1

    def add_many(self, sslist):
        for ss in sslist:
            self.add(ss)


def _get_signatures_from_rust(siglist):
    # this function deals with a disconnect between the way Rust
    # and Python handle signatures; Python expects one
    # minhash (and hence one md5sum) per signature, while
    # Rust supports multiple. For now, go through serializing
    # and deserializing the signature! See issue #1167 for more.
    json_str = sourmash.save_signatures(siglist)
    for ss in sourmash.load_signatures(json_str):
        yield ss


class SaveSignatures_NoOutput(Base_SaveSignaturesToLocation):
    "Do not save signatures."
    def __repr__(self):
        return 'SaveSignatures_NoOutput()'

    @classmethod
    def matches(cls, location):
        return location is None

    def open(self):
        pass

    def close(self):
        pass


class SaveSignatures_Directory(Base_SaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)

    def __repr__(self):
        return f"SaveSignatures_Directory('{self.location}')"

    @classmethod
    def matches(cls, location):
        "anything ending in /"
        if location:
            return location.endswith('/')

    def close(self):
        pass

    def open(self):
        try:
            os.mkdir(self.location)
        except FileExistsError:
            pass
        except:
            notify(f"ERROR: cannot create signature output directory '{self.location}'")
            sys.exit(-1)

    def add(self, ss):
        super().add(ss)
        md5 = ss.md5sum()

        # don't overwrite even if duplicate md5sum
        outname = os.path.join(self.location, f"{md5}.sig.gz")
        if os.path.exists(outname):
            i = 0
            while 1:
                outname = os.path.join(self.location, f"{md5}_{i}.sig.gz")
                if not os.path.exists(outname):
                    break
                i += 1

        with gzip.open(outname, "wb") as fp:
            sigmod.save_signatures([ss], fp, compression=1)


class SaveSignatures_SqliteIndex(Base_SaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)
        self.location = location
        self.idx = None
        self.cursor = None

    @classmethod
    def matches(cls, location):
        "anything ending in .sqldb"
        if location:
            return location.endswith('.sqldb')

    def __repr__(self):
        return f"SaveSignatures_SqliteIndex('{self.location}')"

    def close(self):
        self.idx.commit()
        self.cursor.execute('VACUUM')
        self.idx.close()

    def open(self):
        self.idx = SqliteIndex.create(self.location, append=True)
        self.cursor = self.idx.cursor()

    def add(self, add_sig):
        for ss in _get_signatures_from_rust([add_sig]):
            super().add(ss)
            self.idx.insert(ss, cursor=self.cursor, commit=False)

            # commit every 1000 signatures.
            if self.count % 1000 == 0:
                self.idx.commit()


class SaveSignatures_SigFile(Base_SaveSignaturesToLocation):
    "Save signatures to a .sig JSON file."
    def __init__(self, location):
        super().__init__(location)
        self.keep = []
        self.compress = 0
        if self.location.endswith('.gz'):
            self.compress = 1

    @classmethod
    def matches(cls, location):
        # match anything that is not None or ""
        return bool(location)

    def __repr__(self):
        return f"SaveSignatures_SigFile('{self.location}')"

    def open(self):
        pass

    def close(self):
        if self.location == '-':
            sourmash.save_signatures(self.keep, sys.stdout)
        else:
            # text mode? encode in utf-8
            mode = "w"
            encoding = 'utf-8'

            # compressed? bytes & binary.
            if self.compress:
                encoding = None
                mode = "wb"

            with open(self.location, mode, encoding=encoding) as fp:
                sourmash.save_signatures(self.keep, fp,
                                         compression=self.compress)

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)


class SaveSignatures_ZipFile(Base_SaveSignaturesToLocation):
    "Save compressed signatures in an uncompressed Zip file."
    def __init__(self, location):
        super().__init__(location)
        self.storage = None

    @classmethod
    def matches(cls, location):
        "anything ending in .zip"
        if location:
            return location.endswith('.zip')

    def __repr__(self):
        return f"SaveSignatures_ZipFile('{self.location}')"

    def close(self):
        # finish constructing manifest object & save
        manifest = CollectionManifest(self.manifest_rows)
        manifest_name = "SOURMASH-MANIFEST.csv"

        manifest_fp = StringIO()
        manifest.write_to_csv(manifest_fp, write_header=True)
        manifest_data = manifest_fp.getvalue().encode("utf-8")

        self.storage.save(manifest_name, manifest_data, overwrite=True,
                          compress=True)
        self.storage.flush()
        self.storage.close()

    def open(self):
        from .sbt_storage import ZipStorage

        do_create = True
        if os.path.exists(self.location):
            do_create = False

        storage = None
        try:
            storage = ZipStorage(self.location, mode="w")
        except zipfile.BadZipFile:
            pass

        if storage is None:
            raise ValueError(f"File '{self.location}' cannot be opened as a zip file.")

        if not storage.subdir:
            storage.subdir = 'signatures'

        # now, try to load manifest
        try:
            manifest_data = storage.load('SOURMASH-MANIFEST.csv')
        except (FileNotFoundError, KeyError):
            # if file already exists must have manifest...
            if not do_create:
                raise ValueError(f"Cannot add to existing zipfile '{self.location}' without a manifest")
            self.manifest_rows = []
        else:
            # success! decode manifest_data, create manifest rows => append.
            manifest_data = manifest_data.decode('utf-8')
            manifest_fp = StringIO(manifest_data)
            manifest = CollectionManifest.load_from_csv(manifest_fp)
            self.manifest_rows = list(manifest._select())

        self.storage = storage

    def _exists(self, name):
        try:
            self.storage.load(name)
            return True
        except KeyError:
            return False

    def add(self, add_sig):
        if not self.storage:
            raise ValueError("this output is not open")

        for ss in _get_signatures_from_rust([add_sig]):
            buf = sigmod.save_signatures([ss], compression=1)
            md5 = ss.md5sum()

            storage = self.storage
            path = f'{storage.subdir}/{md5}.sig.gz'
            location = storage.save(path, buf)

            # update manifest
            row = CollectionManifest.make_manifest_row(ss, location,
                                                       include_signature=False)
            self.manifest_rows.append(row)
            super().add(ss)


_save_classes = [
    (10, SaveSignatures_NoOutput),
    (20, SaveSignatures_Directory),
    (30, SaveSignatures_ZipFile),
    (40, SaveSignatures_SqliteIndex),
    (1000, SaveSignatures_SigFile),
]
