"An Abstract Base Class for collections of signatures."

import sourmash
from abc import abstractmethod, ABC
from enum import Enum
from collections import namedtuple
import os


class SearchType(Enum):
    JACCARD = 1
    CONTAINMENT = 2
    MAX_CONTAINMENT = 3
    #ANGULAR_SIMILARITY = 4


def get_search_obj(do_containment, do_max_containment, best_only, threshold):
    if do_containment and do_max_containment:
        raise TypeError("'do_containment' and 'do_max_containment' cannot both be True")

    # configure search - containment? ignore abundance? best only?
    search_cls = IndexSearch
    if best_only:
        search_cls = IndexSearchBestOnly

    if do_containment:
        search_obj = search_cls(SearchType.CONTAINMENT, threshold)
    elif do_max_containment:
        search_obj = search_cls(SearchType.MAX_CONTAINMENT, threshold)
    else:
        search_obj = search_cls(SearchType.JACCARD, threshold)

    return search_obj


def get_gather_obj(query_mh, threshold_bp):
    scaled = query_mh.scaled
    if not scaled: raise TypeError #  @CTB

    # are we setting a threshold?
    threshold=0
    if threshold_bp:
        # if we have a threshold_bp of N, then that amounts to N/scaled
        # hashes:
        n_threshold_hashes = threshold_bp / scaled

        # that then requires the following containment:
        threshold = n_threshold_hashes / len(query_mh)

        # is it too high to ever match? if so, exit.
        if threshold > 1.0:
            return None

    search_obj = IndexSearch(SearchType.CONTAINMENT, threshold=threshold)

    return search_obj

class IndexSearch:
    def __init__(self, search_type, threshold=None):
        score_fn = None
        require_scaled = False

        if search_type == SearchType.JACCARD:
            score_fn = self.score_jaccard
        elif search_type == SearchType.CONTAINMENT:
            score_fn = self.score_containment
            require_scaled = True
        elif search_type == SearchType.MAX_CONTAINMENT:
            score_fn = self.score_max_containment
            require_scaled = True
        self.score_fn = score_fn
        self.require_scaled = require_scaled # @CTB

        if threshold is None:
            threshold = 0
        self.threshold = float(threshold)

    def check_is_compatible(self, sig):
        if self.require_scaled:
            if not sig.minhash.scaled:
                raise TypeError("this search requires a scaled signature")
        if sig.minhash.track_abundance:
            raise TypeError("this search cannot be done with an abund signature")

    def passes(self, score):
        if score and score >= self.threshold:
            return True
        return False

    def collect(self, score):
        pass

    def score_jaccard(self, query_size, shared_size, subject_size, total_size):
        return shared_size / total_size

    def score_containment(self, query_size, shared_size, subject_size,
                          total_size):
        if query_size == 0:
            return 0
        return shared_size / query_size

    def score_max_containment(self, query_size, shared_size, subject_size,
                              total_size):
        min_denom = min(query_size, subject_size)
        if min_denom == 0:
            return 0
        return shared_size / min_denom


class IndexSearchBestOnly(IndexSearch):
    def collect(self, score):
        self.threshold = max(self.threshold, score)


class Index(ABC):
    is_database = False

    @property
    def location(self):
        "Return a resolvable location for this index, if possible."
        return None

    @abstractmethod
    def signatures(self):
        "Return an iterator over all signatures in the Index object."

    @abstractmethod
    def insert(self, signature):
        """ """

    @abstractmethod
    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        """ """

    @classmethod
    @abstractmethod
    def load(cls, location, leaf_loader=None, storage=None, print_version_warning=True):
        """ """

    def find(self, search_fn, query, *args, **kwargs):
        """Use search_fn to find matching signatures in the index.

        search_fn(other_sig, *args) should return a boolean that indicates
        whether other_sig is a match.

        Returns a list.
        """
        search_fn.check_is_compatible(query)
        query_mh = query.minhash

        if query_mh.scaled:
            def downsample(a, b):
                max_scaled = max(a.scaled, b.scaled)
                return a.downsample(scaled=max_scaled), \
                    b.downsample(scaled=max_scaled)
        else:                   # num
            def downsample(a, b):
                min_num = min(a.num, b.num)
                return a.downsample(num=min_num), b.downsample(num=min_num)

        for subj in self.signatures():
            subj_mh = subj.minhash
            if subj_mh.track_abundance:
                subj_mh = subj_mh.flatten()
            qmh, subj_mh = downsample(query_mh, subj_mh)
            query_size = len(qmh)
            subj_size = len(subj_mh)

            # respects num
            merged = qmh + subj_mh
            intersect = set(qmh.hashes) & set(subj_mh.hashes) & set(merged.hashes)
            shared_size = len(intersect)
            total_size = len(merged)

            score = search_fn.score_fn(query_size,
                                       shared_size,
                                       subj_size,
                                       total_size)
            if search_fn.passes(score):
                search_fn.collect(score)
                yield subj, score

    def search_abund(self, query, threshold=None, **kwargs):
        """Return set of matches with angular similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.
        """
        assert query.minhash.track_abundance

        # check arguments
        if threshold is None:
            raise TypeError("'search' requires 'threshold'")
        threshold = float(threshold)

        # do the actual search:
        matches = []
        for subj in self.signatures():
            assert subj.minhash.track_abundance
            score = query.similarity(subj)
            if score >= threshold:
                matches.append((score, subj, self.location))

        # sort!
        matches.sort(key=lambda x: -x[0])
        return matches

    def search(self, query, threshold=None,
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

        search_obj = get_search_obj(do_containment,
                                    do_max_containment,
                                    best_only,
                                    threshold)

        # do the actual search:
        matches = []

        for subj, score in self.find(search_obj, query):
            matches.append((score, subj, self.location))

        # sort!
        matches.sort(key=lambda x: -x[0])
        return matches

    def gather(self, query, *args, **kwargs):
        "Return the match with the best Jaccard containment in the Index."
        if not query.minhash:             # empty query? quit.
            return []

        scaled = query.minhash.scaled
        if not scaled:
            raise ValueError('gather requires scaled signatures')

        threshold_bp = kwargs.get('threshold_bp', 0.0)
        search_obj = get_gather_obj(query.minhash, threshold_bp)
        if not search_obj:
            return []

        # actually do search!
        results = []
        for subj, score in self.find(search_obj, query):
            results.append((score, subj, self.location))

        results.sort(reverse=True, key=lambda x: (x[0], x[1].md5sum()))

        return results[:1]

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


def select_signature(ss, ksize=None, moltype=None, scaled=0, num=0,
                     containment=False):
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

    return True


class LinearIndex(Index):
    "An Index for a collection of signatures. Can load from a .sig file."
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

        return LinearIndex(siglist, self.filename)


class MultiIndex(Index):
    """An Index class that wraps other Index classes.

    The MultiIndex constructor takes two arguments: a list of Index
    objects, and a matching list of sources (filenames, etc.)  If the
    source is not None, then it will be used to override the 'filename'
    in the triple that is returned by search and gather.

    One specific use for this is when loading signatures from a directory;
    MultiIndex will properly record which files provided which signatures.
    """
    def __init__(self, index_list, source_list):
        self.index_list = list(index_list)
        self.source_list = list(source_list)
        assert len(index_list) == len(source_list)

    def signatures(self):
        for idx in self.index_list:
            for ss in idx.signatures():
                yield ss

    def signatures_with_location(self):
        for idx, loc in zip(self.index_list, self.source_list):
            for ss in idx.signatures():
                yield ss, loc

    def __len__(self):
        return sum([ len(idx) for idx in self.index_list ])

    def insert(self, *args):
        raise NotImplementedError

    @classmethod
    def load(self, *args):
        raise NotImplementedError

    @classmethod
    def load_from_path(cls, pathname, force=False):
        "Create a MultiIndex from a path (filename or directory)."
        from .sourmash_args import traverse_find_sigs
        if not os.path.exists(pathname):
            raise ValueError(f"'{pathname}' must be a directory")

        index_list = []
        source_list = []
        for thisfile in traverse_find_sigs([pathname], yield_all_files=force):
            try:
                idx = LinearIndex.load(thisfile)
                index_list.append(idx)
                source_list.append(thisfile)
            except (IOError, sourmash.exceptions.SourmashError):
                if force:
                    continue    # ignore error
                else:
                    raise       # continue past error!

        db = None
        if index_list:
            db = cls(index_list, source_list)
        else:
            raise ValueError(f"no signatures to load under directory '{pathname}'")

        return db

    @classmethod
    def load_from_pathlist(cls, filename):
        "Create a MultiIndex from all files listed in a text file."
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

        db = MultiIndex(idx_list, src_list)
        return db

    def save(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on all indices within this MultiIndex."
        new_idx_list = []
        new_src_list = []
        for idx, src in zip(self.index_list, self.source_list):
            idx = idx.select(**kwargs)
            new_idx_list.append(idx)
            new_src_list.append(src)

        return MultiIndex(new_idx_list, new_src_list)

    def filter(self, filter_fn):
        new_idx_list = []
        new_src_list = []
        for idx, src in zip(self.index_list, self.source_list):
            idx = idx.filter(filter_fn)
            new_idx_list.append(idx)
            new_src_list.append(src)

        return MultiIndex(new_idx_list, new_src_list)

    def search(self, query, *args, **kwargs):
        # do the actual search:
        matches = []
        for idx, src in zip(self.index_list, self.source_list):
            for (score, ss, filename) in idx.search(query, *args, **kwargs):
                best_src = src or filename # override if src provided
                matches.append((score, ss, best_src))
                
        # sort!
        matches.sort(key=lambda x: -x[0])
        return matches

    def gather(self, query, *args, **kwargs):
        "Return the match with the best Jaccard containment in the Index."
        # actually do search!
        results = []
        for idx, src in zip(self.index_list, self.source_list):
            for (score, ss, filename) in idx.gather(query, *args, **kwargs):
                best_src = src or filename # override if src provided
                results.append((score, ss, best_src))
            
        results.sort(reverse=True, key=lambda x: (x[0], x[1].md5sum()))

        return results
