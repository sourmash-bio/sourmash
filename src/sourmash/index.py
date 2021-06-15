"An Abstract Base Class for collections of signatures."

import os
import sourmash
from abc import abstractmethod, ABC
from collections import namedtuple, Counter
import zipfile
import csv
from io import TextIOWrapper

from .search import make_jaccard_search_query, make_gather_query

# generic return tuple for Index.search and Index.gather
IndexSearchResult = namedtuple('Result', 'score, signature, location')

class Index(ABC):
    is_database = False

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
        query_mh = query.minhash
        scaled = query_mh.scaled

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
                     containment=False, picklist=None):
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

    if picklist is not None and ss not in picklist:
        return False

    return True


class LinearIndex(Index):
    """An Index for a collection of signatures. Can load from a .sig file.

    Note: does not use manifests. See LoadedCollection for that functionality.
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
        from .signature import save_signatures
        with open(path, 'wt') as fp:
            save_signatures(self.signatures(), fp)

    @classmethod
    def load(cls, location):
        from .signature import load_signatures
        si = load_signatures(location, do_raise=True)

        lidx = LinearIndex(si, filename=location)
        return lidx

    def select(self, **kwargs):
        """Return new LinearIndex containing only signatures that match req's.

        Does not raise ValueError, but may return an empty Index.
        """
        # eliminate things from kwargs with None or zero value
        kw = { k : v for (k, v) in kwargs.items() if v }

        siglist = []
        for ss in self._signatures:
            if select_signature(ss, **kwargs):
                siglist.append(ss)

        return LinearIndex(siglist, self.location)


class LazyLinearIndex(Index):
    """An Index for lazy linear search of another database.

    The defining feature of this class is that 'find' is inherited
    from the base Index class, which does a linear search with
    signatures().
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
            first_sig = next(iter(self.signatures()))
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
    """
    is_database = True

    def __init__(self, zf, *, selection_dict=None,
                 traverse_yield_all=False, manifest=None):
        self.zf = zf
        self.selection_dict = selection_dict
        self.traverse_yield_all = traverse_yield_all

        # load manifest?
        if manifest is None:
            try:
                zi = self.zf.getinfo('SOURMASH-MANIFEST.csv')
            except KeyError:
                self.manifest = None
            else:
                # CTB: maybe support passing manifest in on constructor?
                # otherwise manifest is loaded on 'select'.
                print(f'found manifest when loading {self.zf.filename}')

                with self.zf.open(zi, 'r') as mfp:
                    # wrap as text, since ZipFile.open only supports 'r' mode.
                    mfp = TextIOWrapper(mfp, 'utf-8')
                    # load manifest!
                    self.manifest = CollectionManifest.load_from_csv(mfp)
        else:
            print(f'using passed-in manifest')
            self.manifest = manifest

    def __bool__(self):
        "Are there any matching signatures in this zipfile? Avoid calling len."
        try:
            first_sig = next(iter(self.signatures()))
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
        return self.zf.filename

    def insert(self, signature):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    @classmethod
    def load(cls, location, traverse_yield_all=False):
        "Class method to load a zipfile."
        zf = zipfile.ZipFile(location, 'r')
        return cls(zf, traverse_yield_all=traverse_yield_all)

    def signatures_with_internal(self):
        # @CTB used for creating manifests
        from .signature import load_signatures
        for zipinfo in self.zf.infolist():
            # should we load this file? if it ends in .sig OR we are forcing:
            if zipinfo.filename.endswith('.sig') or \
               zipinfo.filename.endswith('.sig.gz') or \
               self.traverse_yield_all:
                fp = self.zf.open(zipinfo)

                # now load all the signatures and select on ksize/moltype:
                selection_dict = self.selection_dict

                # note: if 'fp' doesn't contain a valid JSON signature,
                # load_signatures will silently fail & yield nothing.
                for ss in load_signatures(fp):
                    if selection_dict:
                        if select_signature(ss, **self.selection_dict):
                            yield ss, self.zf.filename, zipinfo.filename
                    else:
                        yield ss, self.zf.filename, zipinfo.filename

    def signatures(self):
        "Load all signatures in the zip file."
        from .signature import load_signatures

        manifest = None
        if self.manifest:
            print('.signatures() found manifest!')
            picklist = None
            if self.selection_dict:
                manifest = self.manifest
                def yield_fp():
                    for filename in manifest.select_filenames(**self.selection_dict):
                        zi = self.zf.getinfo(filename)
                        yield self.zf.open(zi)

        if not manifest:
            def yield_fp():
                for zipinfo in self.zf.infolist():
                    # should we load this file? if it ends in .sig OR we are forcing:
                    if zipinfo.filename.endswith('.sig') or \
                       zipinfo.filename.endswith('.sig.gz') or \
                       self.traverse_yield_all:
                        yield self.zf.open(zipinfo)

        for fp in yield_fp():
            # now load all the signatures and select on ksize/moltype:
            selection_dict = self.selection_dict

            # note: if 'fp' doesn't contain a valid JSON signature,
            # load_signatures will silently fail & yield nothing.
            for ss in load_signatures(fp):
                if selection_dict:
                    if select_signature(ss, **self.selection_dict):
                        yield ss
                else:
                    yield ss

    def select(self, **kwargs):
        "Select signatures in zip file based on ksize/moltype/etc."
        return ZipFileLinearIndex(self.zf,
                                  selection_dict=kwargs,
                                  traverse_yield_all=self.traverse_yield_all,
                                  manifest=self.manifest)


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


class LoadedCollection(Index):
    """
    Load a collection of signatures, and retain their original locations.

    One specific use for this is when loading signatures from a directory;
    LoadedCollection will record which specific files provided which
    signatures.

    Creates a manifest on load.

    Note: this is an in-memory collection, and does not do lazy loading:
    all signatures are loaded upon instantiation and kept in memory.
    """
    def __init__(self, manifest):
        self.manifest = manifest

    def signatures(self):
        for row in self.manifest.rows:
            yield row['signature']

    def signatures_with_location(self):
        for row in self.manifest.rows:
            yield row['signature'], row['internal_location']

    def __len__(self):
        return len(self.manifest)

    def insert(self, *args):
        raise NotImplementedError

    @classmethod
    def load(cls, index_list, source_list):
        """Create a LoadedCollection from already-loaded indices.

        Takes two arguments: a list of Index objects, and a matching list
        of source strings (filenames, etc.)  If the source is not None,
        then it will be used to override the location provided by the
        matching Index object.
        """

        # yield all signatures + locations
        def sigloc_iter():
            for idx, loc in zip(index_list, source_list):
                if loc is None:
                    loc = idx.location
                for ss in idx.signatures():
                    yield ss, loc

        # build manifest
        manifest = CollectionManifest.create_manifest(sigloc_iter())

        # create!
        return cls(manifest)

    @classmethod
    def load_from_path(cls, pathname, force=False):
        """
        Create a LoadedCollection from a path (filename or directory).

        Note: this only uses LinearIndex.load(...), so will only load
        signature JSON files.
        """
        from .sourmash_args import traverse_find_sigs
        if not os.path.exists(pathname): # CTB consider changing to isdir...
            raise ValueError(f"'{pathname}' must exist.")

        index_list = []
        source_list = []
        for thisfile in traverse_find_sigs([pathname], yield_all_files=force):
            try:
                idx = LinearIndex.load(thisfile)

                if idx:
                    index_list.append(idx)
                    source_list.append(thisfile)
            except (IOError, sourmash.exceptions.SourmashError):
                if force:
                    continue    # ignore error
                else:
                    raise       # stop loading!

        if not index_list:
            raise ValueError(f"no signatures to load under directory '{pathname}'")

        return cls.load(index_list, source_list)

    @classmethod
    def load_from_pathlist(cls, filename):
        """Create a LoadedCollection from all files listed in a text file.

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
        return LoadedCollection(new_manifest)


class CollectionManifest:
    "@CTB"
    def __init__(self, rows):
        self.rows = tuple(rows)

    def __bool__(self):
        return bool(self.rows)

    def __len__(self):
        return len(self.rows)

    @classmethod
    def load_from_csv(cls, fp):
        "load a manifest from a CSV file."
        manifest_list = []
        r = csv.DictReader(fp)
        for k in ('internal_location', 'md5', 'md5short', 'ksize',
                  'moltype', 'num', 'scaled', 'n_hashes', 'seed',
                  'with_abundance', 'name'):
            if not r.fieldnames:
                return None

            if k not in r.fieldnames:
                raise ValueError(f"missing column '{k}' in manifest.")

        row = None
        for row in r:
            # @CTB revisit
            row['signature'] = None
            manifest_list.append(row)

        return cls(manifest_list)

    @classmethod
    def create_manifest(cls, locations_iter):
        """create a manifest from an iterator that yields (ss, location)

        Stores signatures in manifest rows.

        Note: do NOT catch exceptions here, so this passes through load excs.
        """
        manifest_list = []
        for ss, location in locations_iter:
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
            # @CTB: do we want filename in manifests?
            row['internal_location'] = location
            # @CTB: change name, maybe just make it 'location'

            # CTB: track signature when creating manifest w/this info.
            row['signature'] = ss

            manifest_list.append(row)

        return cls(manifest_list)

    def _select(self, *, ksize=None, moltype=None, scaled=0, num=0,
               containment=False, picklist=None):
        """Yield manifest rows for sigs that match the specified requirements.

        Internal method; call `select_to_manifest` or `select_filenames`
        instead.
        """
        matching_rows = self.rows
        if picklist:
            matching_rows = ( row for row in matching_rows
                              if picklist.matches_siginfo(row) )

        if ksize:
            matching_rows = ( row for row in matching_rows
                              if int(row['ksize']) == ksize )
        if moltype:
            matching_rows = ( row for row in matching_rows
                              if row['moltype'] == moltype )
        if scaled or containment:
            # CTB: check scaled AND containment per select_signature
            # CTB: check num, per select_signature?
            matching_rows = ( row for row in matching_rows
                              if int(row['scaled']) )
        if num:
            # CTB: check scaled, per select_signature?
            matching_rows = ( row for row in matching_rows
                              if int(row['num']) )

        # return only the internal filenames!
        for row in matching_rows:
            yield row

    def select_to_manifest(self, **kwargs):
        "Do a 'select' and return a new CollectionManifest object."
        new_rows = self._select(**kwargs)
        return CollectionManifest(new_rows)

    def select_filenames(self, **kwargs):
        "Do a 'select' and return all of the locations"
        # @CTB return to 'select_to_locations' or something?
        # @CTB or, support lazy loading signatures? or ...?
        for row in self._select(**kwargs):
            yield row['internal_location']

    def __contains__(self, ss):
        "Does this manifest contain this signature?"
        # @CTB currently iterative, we should optimize!
        # @CTB probably should change to 'contains_md5' and
        # @CTB 'contains_filename' or something.
        md5 = ss.md5sum()
        for row in self.rows:
            if md5 == row['md5']:
                return True
        return False
