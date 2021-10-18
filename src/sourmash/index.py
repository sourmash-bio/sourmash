"""An Abstract Base Class for collections of signatures, plus implementations.

APIs and functionality
----------------------

Index classes support three sets of API functionality -

'select(...)', which selects subsets of signatures based on ksize, moltype,
and other criteria, including picklists.

'find(...)', and the 'search', 'gather', and 'counter_gather' implementations
built on top of 'find', which search for signatures that match a query.

'signatures()', which yields all signatures in the Index subject to the
selection criteria.

Classes defined in this file
----------------------------

Index - abstract base class for all Index objects.

LinearIndex - simple in-memory storage of signatures.

LazyLinearIndex - lazy selection and linear search of signatures.

ZipFileLinearIndex - simple on-disk storage of signatures.

class MultiIndex - in-memory storage and selection of signatures from multiple
index objects, using manifests.

LazyLoadedIndex - selection on manifests with loading of index on demand.

CounterGather - an ancillary class returned by the 'counter_gather()' method.
"""

import os
import sourmash
from abc import abstractmethod, ABC
from collections import namedtuple, Counter
import csv
from io import TextIOWrapper

from .search import make_jaccard_search_query, make_gather_query
from .manifest import CollectionManifest
from .logging import debug_literal
from .signature import load_signatures, save_signatures

# generic return tuple for Index.search and Index.gather
IndexSearchResult = namedtuple('Result', 'score, signature, location')

class Index(ABC):
    is_database = False
    manifest = None

    @property
    def location(self):
        "Return a resolvable location for this index, if possible."
        return None

    @abstractmethod
    def signatures(self):
        "Return an iterator over all signatures in the Index object."

    def signatures_with_location(self):
        "Return an iterator over tuples (signature, location) in the Index."
        for ss in self.signatures():
            yield ss, self.location

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, location, internal_location).

        This is an internal API for use in generating manifests, and may
        change without warning.

        This method should be implemented separately for each Index object.
        """
        raise NotImplementedError

    @abstractmethod
    def insert(self, signature):
        """ """

    @abstractmethod
    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        """ """

    @classmethod
    @abstractmethod
    def load(cls, location, leaf_loader=None, storage=None,
             print_version_warning=True):
        """ """

    def find(self, search_fn, query, **kwargs):
        """Use search_fn to find matching signatures in the index.

        search_fn follows the protocol in JaccardSearch objects.

        Returns a list.
        """
        # first: is this query compatible with this search?
        search_fn.check_is_compatible(query)

        # ok! continue!

        # this set of signatures may be heterogenous in scaled/num values;
        # define some processing functions to downsample appropriately.
        query_mh = query.minhash
        assert not query_mh.track_abundance
        if query_mh.scaled:
            # make query and subject compatible w/scaled.
            query_scaled = query_mh.scaled

            def prepare_subject(subj_mh):
                assert subj_mh.scaled
                if subj_mh.track_abundance:
                    subj_mh = subj_mh.flatten()

                # downsample subject to highest scaled
                subj_scaled = subj_mh.scaled
                if subj_scaled < query_scaled:
                    return subj_mh.downsample(scaled=query_scaled)
                else:
                    return subj_mh

            def prepare_query(query_mh, subj_mh):
                assert subj_mh.scaled

                # downsample query to highest scaled
                subj_scaled = subj_mh.scaled
                if subj_scaled > query_scaled:
                    return query_mh.downsample(scaled=subj_scaled)
                else:
                    return query_mh

        else:                   # num
            query_num = query_mh.num

            def prepare_subject(subj_mh):
                assert subj_mh.num
                if subj_mh.track_abundance:
                    subj_mh = subj_mh.flatten()

                # downsample subject to smallest num
                subj_num = subj_mh.num
                if subj_num > query_num:
                    return subj_mh.downsample(num=query_num)
                else:
                    return subj_mh

            def prepare_query(query_mh, subj_mh):
                assert subj_mh.num
                # downsample query to smallest num
                subj_num = subj_mh.num
                if subj_num < query_num:
                    return query_mh.downsample(num=subj_num)
                else:
                    return query_mh

        # now, do the search!
        for subj, location in self.signatures_with_location():
            subj_mh = prepare_subject(subj.minhash)
            # note: we run prepare_query here on the original query minhash.
            query_mh = prepare_query(query.minhash, subj_mh)

            assert not query_mh.track_abundance
            assert not subj_mh.track_abundance

            shared_size, total_size = query_mh.intersection_and_union_size(subj_mh)

            query_size = len(query_mh)
            subj_size = len(subj_mh)

            score = search_fn.score_fn(query_size,
                                       shared_size,
                                       subj_size,
                                       total_size)

            if search_fn.passes(score):
                # note: here we yield the original signature, not the
                # downsampled minhash.
                if search_fn.collect(score, subj):
                    yield IndexSearchResult(score, subj, location)

    def search_abund(self, query, *, threshold=None, **kwargs):
        """Return set of matches with angular similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.
        """
        if not query.minhash.track_abundance:
            raise TypeError("'search_abund' requires query signature with abundance information")

        # check arguments
        if threshold is None:
            raise TypeError("'search_abund' requires 'threshold'")
        threshold = float(threshold)

        # do the actual search:
        matches = []
        for subj, loc in self.signatures_with_location():
            if not subj.minhash.track_abundance:
                raise TypeError("'search_abund' requires subject signatures with abundance information")
            score = query.similarity(subj)
            if score >= threshold:
                matches.append(IndexSearchResult(score, subj, self.location))

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def search(self, query, *, threshold=None,
               do_containment=False, do_max_containment=False,
               best_only=False, **kwargs):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments accepted by all Index subclasses:
          * do_containment: default False. If True, use Jaccard containment.
          * best_only: default False. If True, allow optimizations that
            may. May discard matches better than threshold, but first match
            is guaranteed to be best.
        """
        # check arguments
        if threshold is None:
            raise TypeError("'search' requires 'threshold'")
        threshold = float(threshold)

        search_obj = make_jaccard_search_query(do_containment=do_containment,
                                               do_max_containment=do_max_containment,
                                               best_only=best_only,
                                               threshold=threshold)

        # do the actual search:
        matches = []

        for sr in self.find(search_obj, query, **kwargs):
            matches.append(sr)

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def prefetch(self, query, threshold_bp, **kwargs):
        "Return all matches with minimum overlap."
        if not self:            # empty database? quit.
            raise ValueError("no signatures to search")

        search_fn = make_gather_query(query.minhash, threshold_bp,
                                      best_only=False)

        for sr in self.find(search_fn, query, **kwargs):
            yield sr

    def gather(self, query, threshold_bp=None, **kwargs):
        "Return the match with the best Jaccard containment in the Index."

        results = []
        for result in self.prefetch(query, threshold_bp, **kwargs):
            results.append(result)

        # sort results by best score.
        results.sort(reverse=True,
                     key=lambda x: (x.score, x.signature.md5sum()))

        return results[:1]

    def peek(self, query_mh, threshold_bp=0):
        "Mimic CounterGather.peek() on top of Index. Yes, this is backwards."
        from sourmash import SourmashSignature

        # build a signature to use with self.gather...
        query_ss = SourmashSignature(query_mh)

        # run query!
        try:
            result = self.gather(query_ss, threshold_bp=threshold_bp)
        except ValueError:
            result = None

        if not result:
            return []

        # if matches, calculate intersection & return.
        sr = result[0]
        match_mh = sr.signature.minhash
        scaled = max(query_mh.scaled, match_mh.scaled)
        match_mh = match_mh.downsample(scaled=scaled).flatten()
        query_mh = query_mh.downsample(scaled=scaled)
        intersect_mh = match_mh & query_mh

        return [sr, intersect_mh]

    def consume(self, intersect_mh):
        "Mimic CounterGather.consume on top of Index. Yes, this is backwards."
        pass

    def counter_gather(self, query, threshold_bp, **kwargs):
        """Returns an object that permits 'gather' on top of the
        current contents of this Index.

        The default implementation uses `prefetch` underneath, and returns
        the results in a `CounterGather` object. However, alternate
        implementations need only return an object that meets the
        public `CounterGather` interface, of course.
        """
        # build a flat query
        prefetch_query = query.copy()
        prefetch_query.minhash = prefetch_query.minhash.flatten()

        # find all matches and construct a CounterGather object.
        counter = CounterGather(prefetch_query.minhash)
        for result in self.prefetch(prefetch_query, threshold_bp, **kwargs):
            counter.add(result.signature, result.location)

        # tada!
        return counter

    @abstractmethod
    def select(self, ksize=None, moltype=None, scaled=None, num=None,
               abund=None, containment=None):
        """Return Index containing only signatures that match requirements.

        Current arguments can be any or all of:
        * ksize
        * moltype
        * scaled
        * num
        * containment

        'select' will raise ValueError if the requirements are incompatible
        with the Index subclass.

        'select' may return an empty object or None if no matches can be
        found.
        """


def select_signature(ss, *, ksize=None, moltype=None, scaled=0, num=0,
                     containment=False, abund=None, picklist=None):
    "Check that the given signature matches the specificed requirements."
    # ksize match?
    if ksize and ksize != ss.minhash.ksize:
        return False

    # moltype match?
    if moltype and moltype != ss.minhash.moltype:
        return False

    # containment requires scaled; similarity does not.
    if containment:
        if not scaled:
            raise ValueError("'containment' requires 'scaled' in Index.select'")
        if not ss.minhash.scaled:
            return False

    # 'scaled' and 'num' are incompatible
    if scaled:
        if ss.minhash.num:
            return False
    if num:
        # note, here we check if 'num' is identical; this can be
        # changed later.
        if ss.minhash.scaled or num != ss.minhash.num:
            return False

    if abund:
        # note: minhash w/abund can always be flattened
        if not ss.minhash.track_abundance:
            return False

    if picklist is not None and ss not in picklist:
        return False

    return True


class LinearIndex(Index):
    """An Index for a collection of signatures. Can load from a .sig file.

    Note: See MultiIndex for an in-memory class that uses manifests.

    Concrete class; signatures held in memory; does not use manifests.
    """
    def __init__(self, _signatures=None, filename=None):
        self._signatures = []
        if _signatures:
            self._signatures = list(_signatures)
        self.filename = filename

    @property
    def location(self):
        return self.filename

    def signatures(self):
        return iter(self._signatures)

    def __bool__(self):
        return bool(self._signatures)

    def __len__(self):
        return len(self._signatures)

    def insert(self, node):
        self._signatures.append(node)

    def save(self, path):
        with open(path, 'wt') as fp:
            save_signatures(self.signatures(), fp)

    @classmethod
    def load(cls, location):
        "Load signatures from a JSON signature file."
        si = load_signatures(location, do_raise=True)

        lidx = LinearIndex(si, filename=location)
        return lidx

    def select(self, **kwargs):
        """Return new LinearIndex containing only signatures that match req's.

        Does not raise ValueError, but may return an empty Index.
        """
        siglist = []
        for ss in self._signatures:
            if select_signature(ss, **kwargs):
                siglist.append(ss)

        return LinearIndex(siglist, self.location)


class LazyLinearIndex(Index):
    """An Index for lazy linear search of another database.

    Wrapper class; does not use manifests.

    One of the main purposes of this class is to _force_ linear 'find'
    on index objects. So if this class wraps an SBT, for example, the
    SBT find method will be overriden with the linear 'find' from the
    base class. There are very few situations where this is an improvement,
    so use this class wisely!

    A few notes:
    * selection criteria defined by 'select' are only executed when
      signatures are actually requested (hence, 'lazy').
    * this class stores the provided index 'db' in memory. If you need
      a class that does lazy loading of signatures from disk and does not
      store signatures in memory, see LazyLoadedIndex.
    * if you want efficient manifest-based selection, consider
      MultiIndex (signatures in memory).
    """

    def __init__(self, db, selection_dict={}):
        self.db = db
        self.selection_dict = dict(selection_dict)

    def signatures(self):
        "Return the selected signatures."
        db = self.db.select(**self.selection_dict)
        for ss in db.signatures():
            yield ss

    def signatures_with_location(self):
        "Return the selected signatures, with a location."
        db = self.db.select(**self.selection_dict)
        for tup in db.signatures_with_location():
            yield tup

    def __bool__(self):
        try:
            next(iter(self.signatures()))
            return True
        except StopIteration:
            return False

    def __len__(self):
        raise NotImplementedError

    def insert(self, node):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    @classmethod
    def load(cls, path):
        raise NotImplementedError

    def select(self, **kwargs):
        """Return new object yielding only signatures that match req's.

        Does not raise ValueError, but may return an empty Index.
        """
        selection_dict = dict(self.selection_dict)
        for k, v in kwargs.items():
            if k in selection_dict:
                if selection_dict[k] != v:
                    raise ValueError(f"cannot select on two different values for {k}")
            selection_dict[k] = v

        return LazyLinearIndex(self.db, selection_dict)


class ZipFileLinearIndex(Index):
    """\
    A read-only collection of signatures in a zip file.

    Does not support `insert` or `save`.

    Concrete class; signatures dynamically loaded from disk; uses manifests.
    """
    is_database = True

    def __init__(self, storage, *, selection_dict=None,
                 traverse_yield_all=False, manifest=None, use_manifest=True):
        self.storage = storage
        self.selection_dict = selection_dict
        self.traverse_yield_all = traverse_yield_all
        self.use_manifest = use_manifest

        # do we have a manifest already? if not, try loading.
        if use_manifest:
            if manifest is not None:
                debug_literal('ZipFileLinearIndex using passed-in manifest')
                self.manifest = manifest
            else:
                self._load_manifest()
        else:
            self.manifest = None

        if self.manifest is not None:
            assert not self.selection_dict, self.selection_dict
        if self.selection_dict:
            assert self.manifest is None

    def _load_manifest(self):
        "Load a manifest if one exists"
        try:
            manifest_data = self.storage.load('SOURMASH-MANIFEST.csv')
        except KeyError:
            self.manifest = None
        else:
            debug_literal(f'found manifest on load for {self.storage.path}')

            # load manifest!
            from io import StringIO
            manifest_data = manifest_data.decode('utf-8')
            manifest_fp = StringIO(manifest_data)
            self.manifest = CollectionManifest.load_from_csv(manifest_fp)

    def __bool__(self):
        "Are there any matching signatures in this zipfile? Avoid calling len."
        try:
            next(iter(self.signatures()))
        except StopIteration:
            return False

        return True

    def __len__(self):
        n = 0
        for _ in self.signatures():
            n += 1
        return n

    @property
    def location(self):
        return self.storage.path

    def insert(self, signature):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    @classmethod
    def load(cls, location, traverse_yield_all=False, use_manifest=True):
        "Class method to load a zipfile."
        from .sbt_storage import ZipStorage
        storage = ZipStorage(location)
        return cls(storage, traverse_yield_all=traverse_yield_all,
                   use_manifest=use_manifest)

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, location, internal_location).

        Note: does not limit signatures to subsets.
        """
        zf = self.storage.zipfile

        # list all the files, without using the Storage interface; currently,
        # 'Storage' does not provide a way to list all the files, so :shrug:.
        for zipinfo in zf.infolist():
            # should we load this file? if it ends in .sig OR we are forcing:
            if zipinfo.filename.endswith('.sig') or \
               zipinfo.filename.endswith('.sig.gz') or \
               self.traverse_yield_all:
                fp = zf.open(zipinfo)
                for ss in load_signatures(fp):
                    yield ss, zf.filename, zipinfo.filename

    def signatures(self):
        "Load all signatures in the zip file."
        selection_dict = self.selection_dict
        manifest = None
        if self.manifest is not None:
            manifest = self.manifest
            assert not selection_dict

            # yield all signatures found in manifest
            for filename in manifest.locations():
                data = self.storage.load(filename)
                for ss in load_signatures(data):
                    # in case multiple signatures are in the file, check
                    # to make sure we want to return each one.
                    if ss in manifest:
                        yield ss

        # no manifest! iterate.
        else:
            storage = self.storage
            # if no manifest here, break Storage class encapsulation
            # and go for all the files. (This is necessary to support
            # ad-hoc zipfiles that have no manifests.)
            for zipinfo in storage.zipfile.infolist():
                # should we load this file? if it ends in .sig OR force:
                if zipinfo.filename.endswith('.sig') or \
                   zipinfo.filename.endswith('.sig.gz') or \
                   self.traverse_yield_all:
                    if selection_dict:
                        select = lambda x: select_signature(x,
                                                            **selection_dict)
                    else:
                        select = lambda x: True

                    data = self.storage.load(zipinfo.filename)
                    for ss in load_signatures(data):
                        if select(ss):
                            yield ss

    def select(self, **kwargs):
        "Select signatures in zip file based on ksize/moltype/etc."

        # if we have a manifest, run 'select' on the manifest.
        manifest = self.manifest
        traverse_yield_all = self.traverse_yield_all

        if manifest is not None:
            manifest = manifest.select_to_manifest(**kwargs)
            return ZipFileLinearIndex(self.storage,
                                      selection_dict=None,
                                      traverse_yield_all=traverse_yield_all,
                                      manifest=manifest,
                                      use_manifest=True)
        else:
            # no manifest? just pass along all the selection kwargs to
            # the new ZipFileLinearIndex.

            assert manifest is None
            if self.selection_dict:
                # combine selects...
                d = dict(self.selection_dict)
                for k, v in kwargs.items():
                    if k in d:
                        if d[k] is not None and d[k] != v:
                            raise ValueError(f"incompatible select on '{k}'")
                    d[k] = v
                kwargs = d

            return ZipFileLinearIndex(self.storage,
                                      selection_dict=kwargs,
                                      traverse_yield_all=traverse_yield_all,
                                      manifest=None,
                                      use_manifest=False)


class CounterGather:
    """
    Track and summarize matches for efficient 'gather' protocol.  This
    could be used downstream of prefetch (for example).

    The public interface is `peek(...)` and `consume(...)` only.
    """
    def __init__(self, query_mh):
        if not query_mh.scaled:
            raise ValueError('gather requires scaled signatures')

        # track query
        self.orig_query_mh = query_mh.copy().flatten()
        self.scaled = query_mh.scaled

        # track matching signatures & their locations
        self.siglist = []
        self.locations = []

        # ...and overlaps with query
        self.counter = Counter()

        # cannot add matches once query has started.
        self.query_started = 0

    def add(self, ss, location=None, require_overlap=True):
        "Add this signature in as a potential match."
        if self.query_started:
            raise ValueError("cannot add more signatures to counter after peek/consume")

        # upon insertion, count & track overlap with the specific query.
        overlap = self.orig_query_mh.count_common(ss.minhash, True)
        if overlap:
            i = len(self.siglist)

            self.counter[i] = overlap
            self.siglist.append(ss)
            self.locations.append(location)

            # note: scaled will be max of all matches.
            self.downsample(ss.minhash.scaled)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

    def downsample(self, scaled):
        "Track highest scaled across all possible matches."
        if scaled > self.scaled:
            self.scaled = scaled

    def calc_threshold(self, threshold_bp, scaled, query_size):
        # CTB: this code doesn't need to be in this class.
        threshold = 0.0
        n_threshold_hashes = 0

        if threshold_bp:
            # if we have a threshold_bp of N, then that amounts to N/scaled
            # hashes:
            n_threshold_hashes = float(threshold_bp) / scaled

            # that then requires the following containment:
            threshold = n_threshold_hashes / query_size

        return threshold, n_threshold_hashes

    def peek(self, cur_query_mh, threshold_bp=0):
        "Get next 'gather' result for this database, w/o changing counters."
        self.query_started = 1
        scaled = cur_query_mh.scaled

        # empty? nothing to search.
        counter = self.counter
        if not counter:
            return []

        siglist = self.siglist
        assert siglist

        self.downsample(scaled)
        scaled = self.scaled
        cur_query_mh = cur_query_mh.downsample(scaled=scaled)

        if not cur_query_mh:             # empty query? quit.
            return []

        if cur_query_mh.contained_by(self.orig_query_mh, downsample=True) < 1:
            raise ValueError("current query not a subset of original query")

        # are we setting a threshold?
        threshold, n_threshold_hashes = self.calc_threshold(threshold_bp,
                                                            scaled,
                                                            len(cur_query_mh))
        # is it too high to ever match? if so, exit.
        if threshold > 1.0:
            return []

        # Find the best match -
        most_common = counter.most_common()
        dataset_id, match_size = most_common[0]

        # below threshold? no match!
        if match_size < n_threshold_hashes:
            return []

        ## at this point, we must have a legitimate match above threshold!

        # pull match and location.
        match = siglist[dataset_id]

        # calculate containment
        cont = cur_query_mh.contained_by(match.minhash, downsample=True)
        assert cont
        assert cont >= threshold

        # calculate intersection of this "best match" with query.
        match_mh = match.minhash.downsample(scaled=scaled).flatten()
        intersect_mh = cur_query_mh & match_mh
        location = self.locations[dataset_id]

        # build result & return intersection
        return (IndexSearchResult(cont, match, location), intersect_mh)

    def consume(self, intersect_mh):
        "Remove the given hashes from this counter."
        self.query_started = 1

        if not intersect_mh:
            return

        siglist = self.siglist
        counter = self.counter

        most_common = counter.most_common()

        # Prepare counter for finding the next match by decrementing
        # all hashes found in the current match in other datasets;
        # remove empty datasets from counter, too.
        for (dataset_id, _) in most_common:
            # CTB: note, remaining_mh may not be at correct scaled here.
            # this means that counters that _should_ be empty might not
            # _be_ empty in some situations.  This does not
            # lead to incorrect results, merely potentially overfull
            # 'counter' objects. The tradeoffs to fixing this would
            # need to be examined! (This could be fixed in self.downsample().)
            remaining_mh = siglist[dataset_id].minhash
            intersect_count = intersect_mh.count_common(remaining_mh,
                                                        downsample=True)
            if intersect_count:
                counter[dataset_id] -= intersect_count
                if counter[dataset_id] == 0:
                    del counter[dataset_id]


class MultiIndex(Index):
    """
    Load a collection of signatures, and retain their original locations.

    One specific use for this is when loading signatures from a directory;
    MultiIndex will record which specific files provided which
    signatures.

    Creates a manifest on load.

    Note: this is an in-memory collection, and does not do lazy loading:
    all signatures are loaded upon instantiation and kept in memory.

    Concrete class; signatures held in memory; builds and uses manifests.
    """
    def __init__(self, manifest, parent=""):
        """Constructor; takes manifest containing signatures, together with
        optional top-level location to prepend to internal locations.
        """
        self.manifest = manifest
        self.parent = parent

    def signatures(self):
        for row in self.manifest.rows:
            yield row['signature']

    def signatures_with_location(self):
        for row in self.manifest.rows:
            loc = row['internal_location']
            # here, 'parent' may have been removed from internal_location
            # for directories; if so, add it back in.
            if self.parent:
                loc = os.path.join(self.parent, loc)
            yield row['signature'], loc

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, parent, location)

        CTB note: here, 'internal_location' is the source file for the
        index. This is a special feature of this (in memory) class.
        """
        parent = self.parent
        for row in self.manifest.rows:
            yield row['signature'], parent, row['internal_location']


    def __len__(self):
        return len(self.manifest)

    def insert(self, *args):
        raise NotImplementedError

    @classmethod
    def load(cls, index_list, source_list, parent=""):
        """Create a MultiIndex from already-loaded indices.

        Takes two arguments: a list of Index objects, and a matching list
        of source strings (filenames, etc.)  If the source is not None,
        then it will be used to override the location provided by the
        matching Index object.
        """
        assert len(index_list) == len(source_list)

        # yield all signatures + locations
        def sigloc_iter():
            for idx, iloc in zip(index_list, source_list):
                # override internal location if location is explicitly provided
                if iloc is None:
                    iloc = idx.location
                for ss in idx.signatures():
                    yield ss, iloc

        # build manifest; note, signatures are stored in memory.
        manifest = CollectionManifest.create_manifest(sigloc_iter())

        # create!
        return cls(manifest, parent=parent)

    @classmethod
    def load_from_directory(cls, pathname, *, force=False):
        """Create a MultiIndex from a directory.

        Takes directory path plus optional boolean 'force'. Attempts to
        load all files ending in .sig or .sig.gz, by default; if 'force' is
        True, will attempt to load _all_ files, ignoring errors.
        """
        from .sourmash_args import traverse_find_sigs

        if not os.path.isdir(pathname):
            raise ValueError(f"'{pathname}' must be a directory.")

        index_list = []
        source_list = []

        traversal = traverse_find_sigs([pathname], yield_all_files=force)
        for thisfile in traversal:
            try:
                idx = LinearIndex.load(thisfile)
                index_list.append(idx)

                rel = os.path.relpath(thisfile, pathname)
                source_list.append(rel)
            except (IOError, sourmash.exceptions.SourmashError):
                if force:
                    continue    # ignore error
                else:
                    raise       # stop loading!

        # did we load anything? if not, error
        if not index_list:
            raise ValueError(f"no signatures to load under directory '{pathname}'")

        return cls.load(index_list, source_list, parent=pathname)

    @classmethod
    def load_from_path(cls, pathname, force=False):
        """
        Create a MultiIndex from a path (filename or directory).

        Note: this only uses LinearIndex.load(...), so will only load
        signature JSON files.
        """
        if not os.path.exists(pathname):
            raise ValueError(f"'{pathname}' must exist.")

        if os.path.isdir(pathname): # traverse
            return cls.load_from_directory(pathname, force=force)
        else:                   # load as a .sig/JSON file
            index_list = []
            source_list = []
            try:
                idx = LinearIndex.load(pathname)
                index_list = [idx]
                source_list = [pathname]
            except (IOError, sourmash.exceptions.SourmashError):
                if not force:
                    raise ValueError(f"no signatures to load from '{pathname}'")
                return None

            return cls.load(index_list, source_list)

    @classmethod
    def load_from_pathlist(cls, filename):
        """Create a MultiIndex from all files listed in a text file.

        Note: this will load signatures from directories and databases, too,
        if they are listed in the text file; it uses 'load_file_as_index'
        underneath.
        """
        from .sourmash_args import (load_pathlist_from_file,
                                    load_file_as_index)
        idx_list = []
        src_list = []

        file_list = load_pathlist_from_file(filename)
        for fname in file_list:
            idx = load_file_as_index(fname)
            src = fname

            idx_list.append(idx)
            src_list.append(src)

        return cls.load(idx_list, src_list)

    def save(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on the manifest."
        new_manifest = self.manifest.select_to_manifest(**kwargs)
        return MultiIndex(new_manifest, parent=self.parent)


class LazyLoadedIndex(Index):
    """Given an index location and a manifest, do select only on the manifest
    until signatures are actually requested, and only then load the index.

    This class is useful when you have an index object that consume
    memory when it is loaded (e.g. JSON signature files, or LCA
    databases) and you want to avoid keeping them in memory.  The
    downside of using this class is that it will load the signatures
    from disk every time they are needed (e.g. 'find(...)', 'signatures()').

    Wrapper class; signatures dynamically loaded from disk; uses manifests.
    """
    def __init__(self, filename, manifest):
        "Create an Index with given filename and manifest."
        self.filename = filename
        self.manifest = manifest

    @property
    def location(self):
        "the 'location' attribute for this index will be the filename."
        return self.filename

    def signatures(self):
        "yield all signatures from the manifest."
        if not len(self):
            # nothing in manifest? done!
            return []

        # ok - something in manifest, let's go get those signatures!
        picklist = self.manifest.to_picklist()
        idx = sourmash.load_file_as_index(self.location)

        # convert remaining manifest into picklist
        # CTB: one optimization down the road is, for storage-backed
        # Index objects, to just reach in and get the signatures directly,
        # without going through 'select'. Still, this is nice for abstraction
        # because we don't need to care what the index is - it'll work on
        # anything. It just might be a bit slower.
        idx = idx.select(picklist=picklist)

        # extract signatures.
        for ss in idx.signatures():
            yield ss

    def find(self, *args, **kwargs):
        """Run find after loading and selecting; this provides 'search',
        "'gather', and 'prefetch' functionality, which are built on 'find'.
        """
        if not len(self):
            # nothing in manifest? done!
            return []

        # ok - something in manifest, let's go get those signatures!
        picklist = self.manifest.to_picklist()
        idx = sourmash.load_file_as_index(self.location)

        # convert remaining manifest into picklist
        # CTB: one optimization down the road is, for storage-backed
        # Index objects, to just reach in and get the signatures directly,
        # without going through 'select'. Still, this is nice for abstraction
        # because we don't need to care what the index is - it'll work on
        # anything. It just might be a bit slower.
        idx = idx.select(picklist=picklist)

        for x in idx.find(*args, **kwargs):
            yield x

    def __len__(self):
        "track index size based on the manifest."
        return len(self.manifest)

    def __bool__(self):
        return bool(self.manifest)

    @classmethod
    def load(cls, location):
        """Load index from given location, but retain only the manifest.

        Fail if no manifest.
        """
        idx = sourmash.load_file_as_index(location)
        manifest = idx.manifest

        if not idx.manifest:
            raise ValueError(f"no manifest on index at {location}")

        del idx
        # NOTE: index is not retained outside this scope, just location.

        return cls(location, manifest)

    def insert(self, *args):
        raise NotImplementedError

    def save(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on manifest, return new object with new manifest."
        manifest = self.manifest
        new_manifest = manifest.select_to_manifest(**kwargs)

        return LazyLoadedIndex(self.filename, new_manifest)
