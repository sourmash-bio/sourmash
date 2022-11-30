"""
Utility functions for sourmash CLI commands.

The sourmash_args submodule contains functions that help with various
command-line functions. Library functions in this module often directly
send output to stdout/stderr in support of the CLI, and/or call
sys.exit to exit.

argparse functionality:

* check_scaled_bounds(args) -- check that --scaled is reasonable
* check_num_bounds(args) -- check that --num is reasonable
* get_moltype(args) -- verify that moltype selected is legit
* calculate_moltype(args) -- confirm that only one moltype was selected
* load_picklist(args) -- create a SignaturePicklist from --picklist args
* report_picklist(args, picklist) -- report on picklist value usage/matches
* load_include_exclude_db_patterns(args) -- load --include-db-pattern / --exclude-db-pattern
* apply_picklist_and_pattern(db, ...) -- subselect db on picklist and pattern

signature/database loading functionality:

* load_query_signature(filename, ...) -- load a single signature for query
* traverse_find_sigs(filenames, ...) -- find all .sig and .sig.gz files
* load_dbs_and_sigs(filenames, query, ...) -- load databases & signatures
* load_file_as_index(filename, ...) -- load a sourmash.Index class
* load_file_as_signatures(filename, ...) -- load a list of signatures
* load_pathlist_from_file(filename) -- load a list of paths from a file
* load_many_signatures(locations) -- load many signatures from many files
* get_manifest(idx) -- retrieve or build a manifest from an Index
* class SignatureLoadingProgress - signature loading progress bar

signature and file output functionality:

* SaveSignaturesToLocation(filename) - bulk signature output
* class FileOutput - file output context manager that deals w/stdout well
* class FileOutputCSV - file output context manager for CSV files
"""
import sys
import os
import csv
from enum import Enum
import traceback
import gzip
from io import StringIO, TextIOWrapper
import re
import zipfile
import contextlib

import screed
import sourmash

from sourmash.sbtmh import load_sbt_index
from sourmash.lca.lca_db import load_single_database
import sourmash.exceptions

from .logging import notify, error, debug_literal

from .index import (LinearIndex, ZipFileLinearIndex, MultiIndex)
from .index.sqlite_index import load_sqlite_index, SqliteIndex
from . import signature as sigmod
from .picklist import SignaturePicklist, PickStyle
from .manifest import CollectionManifest
import argparse


DEFAULT_LOAD_K = 31


def check_scaled_bounds(arg):
    f = float(arg)

    if f < 0:
        raise argparse.ArgumentTypeError(f"ERROR: scaled value must be positive")
    if f < 100:
        notify('WARNING: scaled value should be >= 100. Continuing anyway.')
    if f > 1e6:
        notify('WARNING: scaled value should be <= 1e6. Continuing anyway.')
    return f


def check_num_bounds(arg):
    f = int(arg)

    if f < 0:
        raise argparse.ArgumentTypeError(f"ERROR: num value must be positive")
    if f < 50:
        notify('WARNING: num value should be >= 50. Continuing anyway.')
    if f > 50000:
        notify('WARNING: num value should be <= 50000. Continuing anyway.')
    return f


def get_moltype(sig, require=False):
    mh = sig.minhash
    if mh.moltype in ('DNA', 'dayhoff', 'hp', 'protein'):
        moltype = mh.moltype
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig))

    return moltype


def calculate_moltype(args, default=None):
    moltype = default

    n = 0
    if args.dna:
        moltype = 'DNA'
        n += 1
    if args.dayhoff:
        moltype = 'dayhoff'
        n += 1
    if args.hp:
        moltype = 'hp'
        n += 1
    if args.protein:
        moltype = 'protein'
        n += 1

    if n > 1:
        error("cannot specify more than one of --dna/--rna/--nucleotide/--protein/--hp/--dayhoff")
        sys.exit(-1)

    return moltype


def load_picklist(args):
    "Load a SignaturePicklist from --picklist arguments."
    picklist = None
    if args.picklist:
        try:
            picklist = SignaturePicklist.from_picklist_args(args.picklist)

            notify(f"picking column '{picklist.column_name}' of type '{picklist.coltype}' from '{picklist.pickfile}'")

            n_empty_val, dup_vals = picklist.load(picklist.pickfile, picklist.column_name)
        except ValueError as exc:
            error("ERROR: could not load picklist.")
            error(str(exc))
            sys.exit(-1)

        notify(f"loaded {len(picklist.pickset)} distinct values into picklist.")
        if n_empty_val:
            notify(f"WARNING: {n_empty_val} empty values in column '{picklist.column_name}' in picklist file")
        if dup_vals:
            notify(f"WARNING: {len(dup_vals)} values in picklist column '{picklist.column_name}' were not distinct")

    return picklist


def report_picklist(args, picklist):
    if picklist.pickstyle == PickStyle.INCLUDE:
        notify(f"for given picklist, found {len(picklist.found)} matches to {len(picklist.pickset)} distinct values")
        n_missing = len(picklist.pickset - picklist.found)
    elif picklist.pickstyle == PickStyle.EXCLUDE:
        notify(f"for given picklist, found {len(picklist.found)} matches by excluding {len(picklist.pickset)} distinct values")
        n_missing = 0
    if n_missing:
        notify(f"WARNING: {n_missing} missing picklist values.")
        # Note - picklist_require_all is currently only relevant for PickStyle.INCLUDE
        if args.picklist_require_all:
            error("ERROR: failing because --picklist-require-all was set")
            sys.exit(-1)


def load_include_exclude_db_patterns(args):
    if args.picklist and (args.include_db_pattern or args.exclude_db_pattern):
        error("ERROR: --picklist and --include-db-pattern/--exclude cannot be used together.")
        sys.exit(-1)

    if args.include_db_pattern and args.exclude_db_pattern:
        error("ERROR: --include-db-pattern and --exclude-db-pattern cannot be used together.")
        sys.exit(-1)

    if args.include_db_pattern:
        pattern = re.compile(args.include_db_pattern, re.IGNORECASE)
        search_pattern = lambda vals: any(pattern.search(val) for val in vals)
    elif args.exclude_db_pattern:
        pattern = re.compile(args.exclude_db_pattern, re.IGNORECASE)
        search_pattern = lambda vals: all(not pattern.search(val) for val in vals)
    else:
        search_pattern = None

    return search_pattern


def apply_picklist_and_pattern(db, picklist, pattern):
    assert not (picklist and pattern)
    if picklist:
        db = db.select(picklist=picklist)
    elif pattern:
        manifest = db.manifest
        if manifest is None:
            error(f"ERROR on filename '{db.location}'.")
            error("--include-db-pattern/--exclude-db-pattern require a manifest.")
            sys.exit(-1)

        manifest = manifest.filter_on_columns(pattern,
                                              ["name", "filename", "md5"])
        pattern_picklist = manifest.to_picklist()
        db = db.select(picklist=pattern_picklist)

    return db


def load_query_signature(filename, ksize, select_moltype, select_md5=None):
    """Load a single signature to use as a query.

    Uses load_file_as_signatures underneath, so can load from collections
    and indexed databases.
    """
    try:
        sl = load_file_as_signatures(filename, ksize=ksize,
                                     select_moltype=select_moltype)
        sl = list(sl)
    except (OSError, ValueError):
        error(f"Cannot open query file '{filename}'")
        sys.exit(-1)

    if len(sl) and select_md5:
        found_sig = None
        for sig in sl:
            sig_md5 = sig.md5sum()
            if sig_md5.startswith(select_md5.lower()):
                # make sure we pick only one --
                if found_sig is not None:
                    error(f"Error! Multiple signatures start with md5 '{select_md5}'")
                    error("Please use a longer --md5 selector.")
                    sys.exit(-1)
                else:
                    found_sig = sig

            sl = [found_sig]

    if len(sl) and ksize is None:
        ksizes = set([ ss.minhash.ksize for ss in sl ])
        if len(ksizes) == 1:
            ksize = ksizes.pop()
            sl = [ ss for ss in sl if ss.minhash.ksize == ksize ]
            notify(f'select query k={ksize} automatically.')
        elif DEFAULT_LOAD_K in ksizes:
            sl = [ ss for ss in sl if ss.minhash.ksize == DEFAULT_LOAD_K ]
            notify(f'selecting default query k={DEFAULT_LOAD_K}.')
    elif ksize:
        notify(f'selecting specified query k={ksize}')

    if len(sl) != 1:
        error(f"When loading query from '{filename}'", filename)
        error(f'{len(sl)} signatures matching ksize and molecule type;')
        error('need exactly one. Specify --ksize or --dna, --rna, or --protein.')
        sys.exit(-1)

    return sl[0]


def _check_suffix(filename, endings):
    for ending in endings:
        if filename.endswith(ending):
            return True
    return False


def traverse_find_sigs(filenames, yield_all_files=False):
    """Find all .sig and .sig.gz files in & beneath 'filenames'.

    By default, this function returns files with .sig and .sig.gz extensions.
    If 'yield_all_files' is True, this will return _all_ files
    (but not directories).
    """
    endings = ('.sig', '.sig.gz')
    for filename in filenames:
        # check for files in filenames:
        if os.path.isfile(filename):
            if yield_all_files or _check_suffix(filename, endings):
                yield filename

        # filename is a directory -- traverse beneath!
        elif os.path.isdir(filename):
            for root, dirs, files in os.walk(filename):
                for name in sorted(files):
                    fullname = os.path.join(root, name)
                    if yield_all_files or _check_suffix(fullname, endings):
                        yield fullname


def load_dbs_and_sigs(filenames, query, is_similarity_query, *,
                      cache_size=None, picklist=None, pattern=None,
                      fail_on_empty_database=False):
    """
    Load one or more Index objects to search - databases, etc.

    'select' on compatibility with query, and apply picklists & patterns.
    """
    query_mh = query.minhash

    # set selection parameter for containment
    containment = True
    if is_similarity_query:
        containment = False

    databases = []
    total_signatures_loaded = 0
    sum_signatures_after_select = 0
    for filename in filenames:
        notify(f"loading from '{filename}'...", end='\r')

        try:
            db = _load_database(filename, False, cache_size=cache_size)
        except ValueError as e:
            # cannot load database!
            notify(f"ERROR on loading from '{filename}':")
            notify(str(e))
            sys.exit(-1)

        total_signatures_loaded += len(db)

        # get compatible signatures - moltype/ksize/num/scaled
        try:
            db = db.select(moltype=query_mh.moltype,
                           ksize=query_mh.ksize,
                           num=query_mh.num,
                           scaled=query_mh.scaled,
                           containment=containment)
        except ValueError as exc:
            # incompatible collection specified!
            notify(f"ERROR: cannot use '{filename}' for this query.")
            notify(str(exc))
            if fail_on_empty_database:
                sys.exit(-1)
            else:
                db = LinearIndex([])

        # 'select' returns nothing => all signatures filtered out. fail!
        if not db:
            notify(f"no compatible signatures found in '{filename}'")
            if fail_on_empty_database:
                sys.exit(-1)

        sum_signatures_after_select += len(db)

        # last but not least, apply picklist!
        db = apply_picklist_and_pattern(db, picklist, pattern)

        databases.append(db)

    # display num loaded/num selected
    notify("--")
    notify(f"loaded {total_signatures_loaded} total signatures from {len(databases)} locations.")
    notify(f"after selecting signatures compatible with search, {sum_signatures_after_select} remain.")
    print('')

    return databases


def _load_stdin(filename, **kwargs):
    "Load collection from .sig file streamed in via stdin"
    db = None
    if filename == '-':
        # load as LinearIndex, then pass into MultiIndex to generate a
        # manifest.
        lidx = LinearIndex.load(sys.stdin, filename='-')
        db = MultiIndex.load((lidx,), (None,), parent="-")

    return db


def _load_standalone_manifest(filename, **kwargs):
    from sourmash.index import StandaloneManifestIndex

    try:
        idx = StandaloneManifestIndex.load(filename)
    except gzip.BadGzipFile as exc:
        raise ValueError(exc)

    return idx


def _multiindex_load_from_pathlist(filename, **kwargs):
    "Load collection from a list of signature/database files"
    db = MultiIndex.load_from_pathlist(filename)

    return db


def _multiindex_load_from_path(filename, **kwargs):
    "Load collection from a directory."
    traverse_yield_all = kwargs['traverse_yield_all']
    db = MultiIndex.load_from_path(filename, traverse_yield_all)

    return db


def _load_sbt(filename, **kwargs):
    "Load collection from an SBT."
    cache_size = kwargs.get('cache_size')

    try:
        db = load_sbt_index(filename, cache_size=cache_size)
    except (FileNotFoundError, TypeError) as exc:
        raise ValueError(exc)

    return db


def _load_revindex(filename, **kwargs):
    "Load collection from an LCA database/reverse index."
    db, _, _ = load_single_database(filename)
    return db


def _load_sqlite_db(filename, **kwargs):
    return load_sqlite_index(filename)


def _load_zipfile(filename, **kwargs):
    "Load collection from a .zip file."
    db = None
    if filename.endswith('.zip'):
        traverse_yield_all = kwargs['traverse_yield_all']
        try:
            db = ZipFileLinearIndex.load(filename,
                                         traverse_yield_all=traverse_yield_all)
        except FileNotFoundError as exc:
            # turn this into a ValueError => proper exception handling by
            # _load_database.
            raise ValueError(exc)

    return db


# all loader functions, in order.
_loader_functions = [
    ("load from stdin", _load_stdin),
    ("load collection from sqlitedb", _load_sqlite_db),
    ("load from standalone manifest", _load_standalone_manifest),
    ("load from path (file or directory)", _multiindex_load_from_path),
    ("load from file list", _multiindex_load_from_pathlist),
    ("load SBT", _load_sbt),
    ("load revindex", _load_revindex),
    ("load collection from zipfile", _load_zipfile),
    ]


def _load_database(filename, traverse_yield_all, *, cache_size=None):
    """Load file as a database - list of signatures, LCA, SBT, etc.

    Return Index object.

    This is an internal function used by other functions in sourmash_args.
    """
    loaded = False

    # iterate through loader functions, trying them all. Catch ValueError
    # but nothing else.
    for n, (desc, load_fn) in enumerate(_loader_functions):
        try:
            debug_literal(f"_load_databases: trying loader fn {n} '{desc}'")
            db = load_fn(filename,
                         traverse_yield_all=traverse_yield_all,
                         cache_size=cache_size)
        except ValueError:
            debug_literal(f"_load_databases: FAIL on fn {n} {desc}.")
            debug_literal(traceback.format_exc())

        if db is not None:
            loaded = True
            debug_literal("_load_databases: success!")
            break

    # check to see if it's a FASTA/FASTQ record (i.e. screed loadable)
    # so we can provide a better error message to users.
    if not loaded:
        successful_screed_load = False
        it = None
        try:
            # CTB: could be kind of time consuming for a big record, but at the
            # moment screed doesn't expose format detection cleanly.
            with screed.open(filename) as it:
                _ = next(iter(it))
            successful_screed_load = True
        except:
            pass

        if successful_screed_load:
            raise ValueError(f"Error while reading signatures from '{filename}' - got sequences instead! Is this a FASTA/FASTQ file?")

    if not loaded:
        raise ValueError(f"Error while reading signatures from '{filename}'.")

    if loaded:                  # this is a bit redundant but safe > sorry
        assert db is not None

    return db


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


def load_file_as_signatures(filename, *, select_moltype=None, ksize=None,
                            picklist=None,
                            yield_all_files=False,
                            progress=None,
                            pattern=None,
                            _use_manifest=True):
    """Load 'filename' as a collection of signatures. Return an iterable.

    If 'filename' contains an SBT or LCA indexed database, or a regular
    Zip file, will return a signatures() generator. If a Zip file and
    yield_all_files=True, will try to load all files within zip, not just
    .sig files.

    If 'filename' is a JSON file containing one or more signatures, will
    return a list of those signatures.

    If 'filename' is a directory, will load *.sig
    underneath this directory into a list of signatures. If
    yield_all_files=True, will attempt to load all files.

    Applies selector function if select_moltype, ksize or picklist are given.

    'pattern' is a function that returns True on matching values.
    """
    if progress:
        progress.notify(filename)

    db = _load_database(filename, yield_all_files)

    # test fixture ;)
    if not _use_manifest and db.manifest:
        db.manifest = None

    db = db.select(moltype=select_moltype, ksize=ksize)

    # apply pattern search & picklist
    db = apply_picklist_and_pattern(db, picklist, pattern)

    loader = db.signatures()

    if progress is not None:
        return progress.start_file(filename, loader)
    else:
        return loader


def load_pathlist_from_file(filename):
    "Load a list-of-files text file."
    try:
        with open(filename, 'rt') as fp:
            file_list = [ x.rstrip('\r\n') for x in fp ]
        file_list = set(file_list)
        if not file_list:
            raise ValueError("pathlist is empty")
        for checkfile in file_list:
            if not os.path.exists(checkfile):
                raise ValueError(f"file '{checkfile}' inside the pathlist does not exist")
    except IOError:
        raise ValueError(f"pathlist file '{filename}' does not exist")
    except OSError:
        raise ValueError(f"cannot open file '{filename}'")
    except UnicodeDecodeError:
        raise ValueError(f"cannot parse file '{filename}' as list of filenames")
    return file_list


class FileOutput:
    """A context manager for file outputs that handles sys.stdout gracefully.

    Usage:

       with FileOutput(filename, mode) as fp:
          ...

    does what you'd expect, but it handles the situation where 'filename'
    is '-' or None. This makes it nicely compatible with argparse usage,
    e.g.

    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    ...
    with FileOutput(args.output, 'wt') as fp:
       ...

    will properly handle no argument or '-' as sys.stdout.
    """
    def __init__(self, filename, mode='wt', *, newline=None, encoding='utf-8'):
        self.filename = filename
        self.mode = mode
        self.fp = None
        self.newline = newline
        self.encoding = encoding

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        self.fp = open(self.filename, self.mode, newline=self.newline,
                       encoding=self.encoding)
        return self.fp

    def close(self):
        if self.fp is not None: # in case of stdout
            self.fp.close()

    def __enter__(self):
        return self.open()

    def __exit__(self, type, value, traceback):
        # do we need to handle exceptions here?
        if self.fp:
            self.fp.close()

        return False


class FileOutputCSV(FileOutput):
    """A context manager for CSV file outputs.

    Usage:

       with FileOutputCSV(filename) as fp:
          ...

    does what you'd expect, but it handles the situation where 'filename'
    is '-' or None. This makes it nicely compatible with argparse usage,
    e.g.

    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    ...
    with FileOutputCSV(args.output) as w:
       ...

    will properly handle no argument or '-' as sys.stdout.
    """
    def __init__(self, filename):
        self.filename = filename
        self.fp = None

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        if self.filename.endswith('.gz'):
            self.fp = gzip.open(self.filename, 'wt', newline='')
        else:
            self.fp = open(self.filename, 'w', newline='')
        return self.fp


class _DictReader_with_version:
    """A version of csv.DictReader that allows a comment line with a version,
    e.g.

    # SOURMASH-MANIFEST-VERSION: 1.0

    The version is stored as a 2-tuple in the 'version_info' attribute.
    """
    def __init__(self, textfp, *, delimiter=','):
        self.version_info = []

        # is there a '#' in the raw buffer pos 0?
        ch = textfp.buffer.peek(1)

        try:
            ch = ch.decode('utf-8')
        except UnicodeDecodeError:
            raise csv.Error("unable to read CSV file")

        # yes - read a line from the text buffer => parse
        if ch.startswith('#'):
            line = textfp.readline()
            assert line.startswith('# '), line

            # note, this can set version_info to lots of different things.
            # revisit later, I guess. CTB.
            self.version_info = line[2:].strip().split(': ', 2)

        # build a DictReader from the remaining stream
        self.reader = csv.DictReader(textfp, delimiter=delimiter)
        self.fieldnames = self.reader.fieldnames

    def __iter__(self):
        for row in self.reader:
            yield row


@contextlib.contextmanager
def FileInputCSV(filename, *, encoding='utf-8', default_csv_name=None,
                 zipfile_obj=None, delimiter=','):
    """A context manager for reading in CSV files in gzip, zip or text format.

    Assumes comma delimiter, and uses csv.DictReader.

    Note: does not support stdin.

    Note: it seems surprisingly hard to write code that generically handles
    any file handle being passed in; the manifest loading code, in particular,
    uses ZipStorage.load => StringIO obj, which doesn't support peek etc.
    So for now, this context manager is focused on situations where it owns
    the file handle (opens/closes the file).
    """
    fp = None

    if zipfile_obj and not default_csv_name:
        raise ValueError("must provide default_csv_name with a zipfile_obj")

    # first, try to load 'default_csv_name' from a zipfile:
    if default_csv_name:
        # were we given a zipfile obj?
        if zipfile_obj:
            try:
                zi = zipfile_obj.getinfo(default_csv_name)
                with zipfile_obj.open(zi) as fp:
                    textfp = TextIOWrapper(fp,
                                           encoding=encoding,
                                           newline="")
                    r = _DictReader_with_version(textfp, delimiter=delimiter)
                    yield r
            except (zipfile.BadZipFile, KeyError):
                pass # uh oh, we were given a zipfile_obj and it FAILED.

            # no matter what, if given zipfile_obj don't try .gz or regular csv
            return
        else:
            try:
                with zipfile.ZipFile(filename, 'r') as zip_fp:
                    zi = zip_fp.getinfo(default_csv_name)
                    with zip_fp.open(zi) as fp:
                        textfp = TextIOWrapper(fp,
                                               encoding=encoding,
                                               newline="")
                        r = _DictReader_with_version(textfp, delimiter=delimiter)
                        yield r

                # if we got this far with no exceptions, we found
                # the CSV in the zip file. exit generator!
                return
            except (zipfile.BadZipFile, KeyError):
                # no zipfile_obj => it's ok to continue onwards to .gz
                # and regular CSV.
                pass

    # ok, not a zip file - try .gz:
    try:
        with gzip.open(filename, "rt", newline="", encoding=encoding) as fp:
            fp.buffer.peek(1)          # force exception if not a gzip file
            r = _DictReader_with_version(fp, delimiter=delimiter)
            yield r
        return
    except gzip.BadGzipFile:
        pass

    # neither zip nor gz; regular file!
    with open(filename, 'rt', newline="", encoding=encoding) as fp:
        r = _DictReader_with_version(fp, delimiter=delimiter)
        yield r


class SignatureLoadingProgress:
    """A wrapper for signature loading progress reporting.

    Instantiate this class once, and then pass it to load_file_as_signatures
    with progress=<obj>.

    Alternatively, call obj.start_file(location, iter) each time you
    start loading signatures from a new file via iter.

    You can optionally notify of reading a file with `.notify(location)`.
    """
    def __init__(self, reporting_interval=10):
        self.n_sig = 0
        self.interval = reporting_interval
        self.screen_width = 79

    def __len__(self):
        return self.n_sig

    def short_notify(self, msg_template, *args, **kwargs):
        """Shorten the notification message so that it fits on one line.

        Good for repeating notifications with end='\r' especially...
        """

        msg = msg_template.format(*args, **kwargs)
        end = kwargs.get('end', '\n')
        w = self.screen_width

        if len(msg) > w:
            truncate_len = len(msg) - w + 3
            msg = '<<<' + msg[truncate_len:]

        notify(msg, end=end)

    def notify(self, location):
        self.short_notify(f"...{self.n_sig} sigs so far. Now reading from file '{location}'", end='\r')

    def start_file(self, location, loader):
        n_this = 0
        n_before = self.n_sig

        try:
            for result in loader:
                # track n from this file, as well as total n
                n_this += 1
                n_total = n_before + n_this
                if n_this and n_total % self.interval == 0:
                    self.short_notify("...loading from '{}' / {} sigs total",
                                      location, n_total, end='\r')

                yield result
        except KeyboardInterrupt:
            # might as well nicely handle CTRL-C while we're at it!
            notify('\n(CTRL-C received! quitting.)')
            sys.exit(-1)
        finally:
            self.n_sig += n_this

        self.short_notify(f"Loaded {n_this} sigs from '{location}'",
                          end='\r')


def load_many_signatures(locations, progress, *, yield_all_files=False,
                         ksize=None, moltype=None, picklist=None, force=False,
                         pattern=None):
    """
    Load many signatures from multiple files, with progress indicators.

    Takes ksize, moltype, and picklist selectors.

    If 'yield_all_files=True' then tries to load all files in specified
    directories.

    If 'force=True' then continues past survivable errors.

    Yields (sig, location) tuples.
    """
    for loc in locations:
        try:
            # open index,
            idx = load_file_as_index(loc, yield_all_files=yield_all_files)
            idx = idx.select(ksize=ksize, moltype=moltype)

            idx = apply_picklist_and_pattern(idx, picklist, pattern)

            # start up iterator,
            loader = idx.signatures_with_location()

            # go!
            n = 0               # count signatures loaded
            for sig, sigloc in progress.start_file(loc, loader):
                yield sig, sigloc
                n += 1
            notify(f"loaded {n} signatures from '{loc}'", end='\r')
        except ValueError as exc:
            # trap expected errors, and either power through or display + exit.
            if force:
                notify(f"ERROR: {str(exc)}")
                notify("(continuing)")
                continue
            else:
                notify(f"ERROR: {str(exc)}")
                sys.exit(-1)
        except KeyboardInterrupt:
            notify("Received CTRL-C - exiting.")
            sys.exit(-1)

    n_files = len(locations)
    notify(f"loaded {len(progress)} signatures total, from {n_files} files")


def get_manifest(idx, *, require=True, rebuild=False):
    """
    Retrieve a manifest for this idx, loaded with `load_file_as_index`.

    Even if a manifest exists and `rebuild` is True, rebuild the manifest.
    If a manifest does not exist or `rebuild` is True, try to build one.
    If a manifest cannot be built and `require` is True, error exit.

    In the case where `require=False` and a manifest cannot be built,
    may return None. Otherwise always returns a manifest.
    """
    from sourmash.index import CollectionManifest

    m = idx.manifest

    # has one, and don't want to rebuild? easy! return!
    if m is not None and not rebuild:
        debug_literal("get_manifest: found manifest")
        return m

    debug_literal(f"get_manifest: no manifest found / rebuild={rebuild}")

    # need to build one...
    try:
        notify("Generating a manifest...")
        m = CollectionManifest.create_manifest(idx._signatures_with_internal(),
                                               include_signature=False)
        debug_literal("get_manifest: rebuilt manifest.")
    except NotImplementedError:
        if require:
            error(f"ERROR: manifests cannot be generated for {idx.location}")
            sys.exit(-1)
        else:
            debug_literal("get_manifest: cannot build manifest, not req'd")
            return None

    return m

#
# enum and classes for saving signatures progressively
#

def _get_signatures_from_rust(siglist):
    for ss in siglist:
        try:
            ss.md5sum()
            yield ss
        except sourmash.exceptions.Panic:
            # this deals with a disconnect between the way Rust
            # and Python handle signatures; Python expects one
            # minhash (and hence one md5sum) per signature, while
            # Rust supports multiple. For now, go through serializing
            # and deserializing the signature! See issue #1167 for more.
            json_str = sourmash.save_signatures([ss])
            for ss in sourmash.load_signatures(json_str):
                yield ss


class _BaseSaveSignaturesToLocation:
    "Base signature saving class. Track location (if any) and count."
    def __init__(self, location):
        self.location = location
        self.count = 0

    def __repr__(self):
        raise NotImplementedError

    def __len__(self):
        return self.count

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


class SaveSignatures_NoOutput(_BaseSaveSignaturesToLocation):
    "Do not save signatures."
    def __repr__(self):
        return 'SaveSignatures_NoOutput()'

    def open(self):
        pass

    def close(self):
        pass


class SaveSignatures_Directory(_BaseSaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)

    def __repr__(self):
        return f"SaveSignatures_Directory('{self.location}')"

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


class SaveSignatures_SqliteIndex(_BaseSaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)
        self.location = location
        self.idx = None
        self.cursor = None

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


class SaveSignatures_SigFile(_BaseSaveSignaturesToLocation):
    "Save signatures to a .sig JSON file."
    def __init__(self, location):
        super().__init__(location)
        self.keep = []
        self.compress = 0
        if self.location.endswith('.gz'):
            self.compress = 1

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


class SaveSignatures_ZipFile(_BaseSaveSignaturesToLocation):
    "Save compressed signatures in an uncompressed Zip file."
    def __init__(self, location):
        super().__init__(location)
        self.storage = None

    def __repr__(self):
        return f"SaveSignatures_ZipFile('{self.location}')"

    def close(self):
        # finish constructing manifest object & save
        manifest = CollectionManifest(self.manifest_rows)
        manifest_name = f"SOURMASH-MANIFEST.csv"

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


class SigFileSaveType(Enum):
    NO_OUTPUT = 0
    SIGFILE = 1
    SIGFILE_GZ = 2
    DIRECTORY = 3
    ZIPFILE = 4
    SQLITEDB = 5

_save_classes = {
    SigFileSaveType.NO_OUTPUT: SaveSignatures_NoOutput,
    SigFileSaveType.SIGFILE: SaveSignatures_SigFile,
    SigFileSaveType.SIGFILE_GZ: SaveSignatures_SigFile,
    SigFileSaveType.DIRECTORY: SaveSignatures_Directory,
    SigFileSaveType.ZIPFILE: SaveSignatures_ZipFile,
    SigFileSaveType.SQLITEDB: SaveSignatures_SqliteIndex,
}


def SaveSignaturesToLocation(filename, *, force_type=None):
    """Create and return an appropriate object for progressive saving of
    signatures."""
    save_type = None
    if not force_type:
        if filename is None:
            save_type = SigFileSaveType.NO_OUTPUT
        elif filename.endswith('/'):
            save_type = SigFileSaveType.DIRECTORY
        elif filename.endswith('.gz'):
            save_type = SigFileSaveType.SIGFILE_GZ
        elif filename.endswith('.zip'):
            save_type = SigFileSaveType.ZIPFILE
        elif filename.endswith('.sqldb'):
            save_type = SigFileSaveType.SQLITEDB
        else:
            # default to SIGFILE intentionally!
            save_type = SigFileSaveType.SIGFILE
    else:
        save_type = force_type

    cls = _save_classes.get(save_type)
    if cls is None:
        raise Exception("invalid save type; this should never happen!?")

    return cls(filename)
