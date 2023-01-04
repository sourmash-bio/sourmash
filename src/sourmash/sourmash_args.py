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
* load_pathlist_from_file(filename) -- load a list of paths from a file
* load_many_signatures(locations) -- load many signatures from many files
* get_manifest(idx) -- retrieve or build a manifest from an Index
* class SignatureLoadingProgress - signature loading progress bar
* load_file_as_signatures(filename, ...) -- load a list of signatures

signature and file output functionality:

* class FileOutput - file output context manager that deals w/stdout well
* class FileOutputCSV - file output context manager for CSV files
"""
import sys
import os
import csv
import gzip
from io import TextIOWrapper
import re
import zipfile
import contextlib
import argparse

from .logging import notify, error, debug_literal

from .index import LinearIndex
from .picklist import SignaturePicklist, PickStyle
from .manifest import CollectionManifest
from .save_load import (SaveSignaturesToLocation, load_file_as_index,
                        _load_database)


DEFAULT_LOAD_K = 31


def check_scaled_bounds(arg):
    f = float(arg)

    if f < 0:
        raise argparse.ArgumentTypeError("ERROR: scaled value must be positive")
    if f < 100:
        notify('WARNING: scaled value should be >= 100. Continuing anyway.')
    if f > 1e6:
        notify('WARNING: scaled value should be <= 1e6. Continuing anyway.')
    return f


def check_num_bounds(arg):
    f = int(arg)

    if f < 0:
        raise argparse.ArgumentTypeError("ERROR: num value must be positive")
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
