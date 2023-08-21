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

MultiIndex - in-memory storage and selection of signatures from multiple
index objects, using manifests. All signatures are kept in memory.

StandaloneManifestIndex - load manifests directly, and do lazy loading of
signatures on demand. No signatures are kept in memory.

CounterGather - an ancillary class returned by the 'counter_gather()' method.
"""

import os
import sourmash
from abc import abstractmethod, ABC
from collections import namedtuple, Counter

from sourmash.search import (make_jaccard_search_query,
                             make_containment_query,
                             calc_threshold_from_bp)
from sourmash.manifest import CollectionManifest
from sourmash.logging import debug_literal
from sourmash.signature import load_signatures, save_signatures
from sourmash.minhash import (flatten_and_downsample_scaled,
                              flatten_and_downsample_num,
                              flatten_and_intersect_scaled)

# generic return tuple for Index.search and Index.gather
IndexSearchResult = namedtuple('Result', 'score, signature, location')

class Index(ABC):
    # this will be removed soon; see sourmash#1894.
    is_database = False

    # 'manifest', when set, implies efficient selection and direct
    # access to signatures. Signatures may be stored in the manifest
    # or loaded on demand from disk depending on the class, however.
    manifest = None

    @abstractmethod
    def __len__(self):
        "Return the number of signatures in this Index object."

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
        """Return an iterator of tuples (ss, internal_location).

        Unlike 'signatures_with_location()', this iterator should return
        _all_ signatures in the object, not just those that remain after
        selection/filtering.

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

        Generator. Returns 0 or more IndexSearchResult objects.
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
                return flatten_and_downsample_scaled(subj_mh, query_scaled)

            def prepare_query(query_mh, subj_mh):
                return flatten_and_downsample_scaled(query_mh, subj_mh.scaled)

        else:                   # num
            query_num = query_mh.num

            def prepare_subject(subj_mh):
                return flatten_and_downsample_num(subj_mh, query_num)

            def prepare_query(query_mh, subj_mh):
                return flatten_and_downsample_num(query_mh, subj_mh.num)

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
        """Return list of IndexSearchResult with angular similarity above 'threshold'.

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
            score = query.similarity(subj, downsample=True)
            if score >= threshold:
                matches.append(IndexSearchResult(score, subj, loc))

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def search(self, query, *, threshold=None,
               do_containment=False, do_max_containment=False,
               best_only=False, **kwargs):
        """Return list of IndexSearchResult with similarity above 'threshold'.

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
        matches = list(self.find(search_obj, query, **kwargs))

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def prefetch(self, query, threshold_bp, **kwargs):
        """Return all matches with minimum overlap.

        Generator. Returns 0 or more IndexSearchResult namedtuples.
        """
        if not self:            # empty database? quit.
            raise ValueError("no signatures to search")

        # default best_only to False
        best_only = kwargs.get('best_only', False)

        search_fn = make_containment_query(query.minhash, threshold_bp,
                                           best_only=best_only)

        for sr in self.find(search_fn, query, **kwargs):
            yield sr

    def best_containment(self, query, threshold_bp=None, **kwargs):
        """Return the match with the best Jaccard containment in the Index.

        Returns an IndexSearchResult namedtuple or None.
        """

        results = self.prefetch(query, threshold_bp, best_only=True, **kwargs)
        results = sorted(results,
                         key=lambda x: (-x.score, x.signature.md5sum()))

        try:
            return next(iter(results))
        except StopIteration:
            return None

    def peek(self, query_mh, *, threshold_bp=0):
        """Mimic CounterGather.peek() on top of Index.

        This is implemented for situations where we don't want to use
        'prefetch' functionality. It is a light wrapper around the
        'best_containment(...)' method.
        """
        from sourmash import SourmashSignature

        # build a signature to use with self.gather...
        query_ss = SourmashSignature(query_mh)

        # run query!
        try:
            result = self.best_containment(query_ss, threshold_bp=threshold_bp)
        except ValueError:
            result = None

        if not result:
            return []

        # if matches, calculate intersection & return.
        intersect_mh = flatten_and_intersect_scaled(result.signature.minhash,
                                                    query_mh)

        return [result, intersect_mh]

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
        with query.update() as prefetch_query:
            prefetch_query.minhash = prefetch_query.minhash.flatten()

        # find all matches and construct a CounterGather object.
        counter = CounterGather(prefetch_query)
        for result in self.prefetch(prefetch_query, threshold_bp, **kwargs):
            counter.add(result.signature, location=result.location)

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
    "Check that the given signature matches the specified requirements."
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
    def load(cls, location, filename=None):
        "Load signatures from a JSON signature file."
        si = load_signatures(location, do_raise=True)

        if filename is None:
            filename=location
        lidx = LinearIndex(si, filename=filename)
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
      store signatures in memory, see StandaloneManifestIndex.
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
        db = self.db.select(**self.selection_dict)
        return len(db)

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
        except (KeyError, FileNotFoundError):
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
        "calculate number of signatures."

        # use manifest, if available.
        m = self.manifest
        if self.manifest is not None:
            return len(m)

        # otherwise, iterate across all signatures.
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
        from ..sbt_storage import ZipStorage

        # we can only load from existing zipfiles in this method.
        if not os.path.exists(location):
            raise FileNotFoundError(location)

        storage = ZipStorage(location)
        return cls(storage, traverse_yield_all=traverse_yield_all,
                   use_manifest=use_manifest)

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, internal_location).

        Note: does not limit signatures to subsets.
        """
        # list all the files, without using the Storage interface; currently,
        # 'Storage' does not provide a way to list all the files, so :shrug:.
        for filename in self.storage._filenames():
            # should we load this file? if it ends in .sig OR we are forcing:
            if filename.endswith('.sig') or \
               filename.endswith('.sig.gz') or \
               self.traverse_yield_all:
                sig_data = self.storage.load(filename)
                for ss in load_signatures(sig_data):
                    yield ss, filename

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
            for filename in storage._filenames():
                # should we load this file? if it ends in .sig OR force:
                if filename.endswith('.sig') or \
                   filename.endswith('.sig.gz') or \
                   self.traverse_yield_all:
                    if selection_dict:
                        select = lambda x: select_signature(x,
                                                            **selection_dict)
                    else:
                        select = lambda x: True

                    data = self.storage.load(filename)
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
    """This is an ancillary class that is used to implement "fast
    gather", post-prefetch. It tracks and summarize matches for
    efficient min-set-cov/'gather'.

    The class constructor takes a query MinHash that must be scaled, and
    then takes signatures that have overlaps with the query (via 'add').

    After all overlapping signatures have been loaded, the 'peek'
    method is then used at each stage of the 'gather' procedure to
    find the best match, and the 'consume' method is used to remove
    a match from this counter.

    This particular implementation maintains a collections.Counter that
    is used to quickly find the best match when 'peek' is called, but
    other implementations are possible ;).

    Note that redundant matches (SourmashSignature objects) with
    duplicate md5s are collapsed inside the class, because we use the
    md5sum as a key into the dictionary used to store matches.
    """
    def __init__(self, query):
        "Constructor - takes a query SourmashSignature."
        query_mh = query.minhash
        if not query_mh.scaled:
            raise ValueError('gather requires scaled signatures')

        # track query
        self.orig_query_mh = query_mh.copy().flatten()
        self.scaled = query_mh.scaled

        # use these to track loaded matches & their locations
        self.siglist = {}
        self.locations = {}

        # ...and also track overlaps with the progressive query
        self.counter = Counter()

        # fence to make sure we do add matches once query has started.
        self.query_started = 0

    def add(self, ss, *, location=None, require_overlap=True):
        "Add this signature in as a potential match."
        if self.query_started:
            raise ValueError("cannot add more signatures to counter after peek/consume")

        # upon insertion, count & track overlap with the specific query.
        overlap = self.orig_query_mh.count_common(ss.minhash, True)
        if overlap:
            md5 = ss.md5sum()

            self.counter[md5] = overlap
            self.siglist[md5] = ss
            self.locations[md5] = location

            # note: scaled will be max of all matches.
            self.downsample(ss.minhash.scaled)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

    def downsample(self, scaled):
        "Track highest scaled across all possible matches."
        if scaled > self.scaled:
            self.scaled = scaled
        return self.scaled

    def signatures(self):
        "Return all signatures."
        for ss in self.siglist.values():
            yield ss

    @property
    def union_found(self):
        """Return a MinHash containing all found hashes in the query.

        This calculates the union of the found matches, intersected
        with the original query.
        """
        orig_query_mh = self.orig_query_mh

        # create empty MinHash from orig query
        found_mh = orig_query_mh.copy_and_clear()

        # for each match, intersect match with query & then add to found_mh.
        for ss in self.siglist.values():
            intersect_mh = flatten_and_intersect_scaled(ss.minhash,
                                                        orig_query_mh)
            found_mh.add_many(intersect_mh)

        return found_mh

    def peek(self, cur_query_mh, *, threshold_bp=0):
        "Get next 'gather' result for this database, w/o changing counters."
        self.query_started = 1

        # empty? nothing to search.
        counter = self.counter
        if not counter:
            return []

        siglist = self.siglist
        assert siglist

        scaled = self.downsample(cur_query_mh.scaled)
        cur_query_mh = cur_query_mh.downsample(scaled=scaled)

        if not cur_query_mh:             # empty query? quit.
            return []

        # CTB: could probably remove this check unless debug requested.
        if cur_query_mh.contained_by(self.orig_query_mh, downsample=True) < 1:
            raise ValueError("current query not a subset of original query")

        # are we setting a threshold?
        try:
            x = calc_threshold_from_bp(threshold_bp, scaled, len(cur_query_mh))
            threshold, n_threshold_hashes = x
        except ValueError:
            # too high to ever match => exit
            return []

        # Find the best match using the internal Counter.
        most_common = counter.most_common()
        dataset_id, match_size = most_common[0]

        # below threshold? no match!
        if match_size < n_threshold_hashes:
            return []

        ## at this point, we have a legitimate match above threshold!

        # pull match and location.
        match = siglist[dataset_id]

        # calculate containment
        # CTB: this check is probably redundant with intersect_mh calc, below.
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
        "Maintain the internal counter by removing the given hashes."
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

    There are a variety of loading functions:
    * `load` takes a list of already-loaded Index objects,
      together with a list of their locations.
    * `load_from_directory` traverses a directory to load files within.
    * `load_from_path` takes an arbitrary pathname and tries to load it
      as a directory, or as a .sig file.
    * `load_from_pathlist` takes a text file full of pathnames and tries
      to load them all.

    Concrete class; signatures held in memory; builds and uses manifests.
    """
    def __init__(self, manifest, parent, *, prepend_location=False):
        """Constructor; takes manifest containing signatures, together with
        the top-level location.
        """
        self.manifest = manifest
        self.parent = parent
        self.prepend_location = prepend_location

        if prepend_location and self.parent is None:
            raise ValueError("must set 'parent' if 'prepend_location' is set")

    @property
    def location(self):
        return self.parent

    def signatures(self):
        for row in self.manifest.rows:
            yield row['signature']

    def signatures_with_location(self):
        for row in self.manifest.rows:
            loc = row['internal_location']
            # here, 'parent' may have been removed from internal_location
            # for directories; if so, add it back in.
            if self.prepend_location:
                loc = os.path.join(self.parent, loc)
            yield row['signature'], loc

    def _signatures_with_internal(self):
        """Return an iterator of tuples (ss, location)

        CTB note: here, 'internal_location' is the source file for the
        index. This is a special feature of this (in memory) class.
        """
        for row in self.manifest.rows:
            yield row['signature'], row['internal_location']


    def __len__(self):
        if self.manifest is None:
            return 0

        return len(self.manifest)

    def insert(self, *args):
        raise NotImplementedError

    @classmethod
    def load(cls, index_list, source_list, parent, *, prepend_location=False):
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

        # build manifest; note, ALL signatures are stored in memory.
        # CTB: could do this on demand?
        # CTB: should we use get_manifest functionality?
        # CTB: note here that the manifest is created by iteration
        # *even if it already exists.* This could be changed to be more
        # efficient... but for now, use StandaloneManifestIndex if you
        # want to avoid this when loading from multiple files.
        manifest = CollectionManifest.create_manifest(sigloc_iter())

        # create!
        return cls(manifest, parent, prepend_location=prepend_location)

    @classmethod
    def load_from_directory(cls, pathname, *, force=False):
        """Create a MultiIndex from a directory.

        Takes directory path plus optional boolean 'force'. Attempts to
        load all files ending in .sig or .sig.gz, by default; if 'force' is
        True, will attempt to load _all_ files, ignoring errors.

        Will not load anything other than JSON signature files.
        """
        from ..sourmash_args import traverse_find_sigs

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
            except (IOError, sourmash.exceptions.SourmashError) as exc:
                if force:
                    continue    # ignore error
                else:
                    raise ValueError(exc)      # stop loading!

        # did we load anything? if not, error
        if not index_list:
            raise ValueError(f"no signatures to load under directory '{pathname}'")

        return cls.load(index_list, source_list, pathname,
                        prepend_location=True)

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

        # load as a .sig/JSON file
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

        return cls.load(index_list, source_list, pathname)

    @classmethod
    def load_from_pathlist(cls, filename):
        """Create a MultiIndex from all files listed in a text file.

        Note: this will attempt to load signatures from each file,
        including zip collections, etc; it uses 'load_file_as_index'
        underneath.
        """
        from ..sourmash_args import (load_pathlist_from_file,
                                    load_file_as_index)
        idx_list = []
        src_list = []

        file_list = load_pathlist_from_file(filename)
        for fname in file_list:
            idx = load_file_as_index(fname)
            src = fname

            idx_list.append(idx)
            src_list.append(src)

        return cls.load(idx_list, src_list, filename)

    def save(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on the manifest."
        new_manifest = self.manifest.select_to_manifest(**kwargs)
        return MultiIndex(new_manifest, self.parent,
                          prepend_location=self.prepend_location)


class StandaloneManifestIndex(Index):
    """Load a standalone manifest as an Index.

    This class is useful for the situation where you have a directory
    with many signature collections underneath it, and you don't want to load
    every collection each time you run sourmash.

    Instead, you can run 'sourmash sig manifest <directory> -o mf.csv' to
    output a manifest and then use this class to load 'mf.csv' directly.
    Sketch type selection, picklists, and pattern matching will all work
    directly on the manifest and will load signatures only upon demand.

    One feature of this class is that absolute paths to sketches in
    the 'internal_location' field of the manifests will be loaded properly.
    This permits manifests to be constructed for various collections of
    signatures that reside elsewhere, and not just below a single directory
    prefix.

    StandaloneManifestIndex does _not_ store signatures in memory.

    This class also overlaps in concept with MultiIndex when
    MultiIndex.load_from_pathlist is used to load other Index
    objects. However, this class does not store any signatures in
    memory, unlike MultiIndex.
    """
    is_database = True

    def __init__(self, manifest, location, *, prefix=None):
        """Create object. 'location' is path of manifest file, 'prefix' is
        prepended to signature paths when loading non-abspaths."""
        assert manifest is not None
        self.manifest = manifest
        self._location = location
        self.prefix = prefix

    @classmethod
    def load(cls, location, *, prefix=None):
        """Load manifest file from given location.

        If prefix is None (default), it is automatically set from dirname.
        Set prefix='' to avoid this, or provide an explicit prefix.
        """
        if not os.path.isfile(location):
            raise ValueError(f"provided manifest location '{location}' is not a file")

        m = CollectionManifest.load_from_filename(location)

        if prefix is None:
            prefix = os.path.dirname(location)

        return cls(m, location, prefix=prefix)

    @property
    def location(self):
        "Return the path to this manifest."
        return self._location

    def signatures_with_location(self):
        "Return an iterator over all signatures and their locations."
        for ss, loc in self._signatures_with_internal():
            yield ss, loc

    def signatures(self):
        "Return an iterator over all signatures."
        for ss, loc in self._signatures_with_internal():
            yield ss

    def _signatures_with_internal(self):
        """Return an iterator over all sigs of (sig, internal_location)

        Note that this is implemented differently from most Index
        objects in that it only lists subselected parts of the
        manifest, and not the original manifest. This was done out of
        convenience: we don't currently have access to the original
        manifest in this class.
        """
        # collect all internal locations
        picklist = self.manifest.to_picklist()
        for iloc in self.manifest.locations():
            # prepend location with prefix?
            if not iloc.startswith('/') and self.prefix:
                iloc = os.path.join(self.prefix, iloc)

            idx = sourmash.load_file_as_index(iloc)
            idx = idx.select(picklist=picklist)
            for ss in idx.signatures():
                yield ss, iloc

    def __len__(self):
        "Number of signatures in this manifest (after any select)."
        return len(self.manifest)

    def __bool__(self):
        "Is this manifest empty?"
        return bool(self.manifest)

    def save(self, *args):
        raise NotImplementedError

    def insert(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on the manifest."
        new_manifest = self.manifest.select_to_manifest(**kwargs)
        return StandaloneManifestIndex(new_manifest, self._location,
                                       prefix=self.prefix)
