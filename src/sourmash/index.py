"An Abstract Base Class for collections of signatures."

from abc import abstractmethod, ABC
from enum import Enum
from collections import namedtuple


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
            qmh, subj_mh = downsample(query_mh, subj.minhash)
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

    def search(self, query, threshold=None,
               do_containment=False, do_max_containment=False,
               ignore_abundance=False, best_only=False, **kwargs):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments accepted by all Index subclasses:
          * do_containment: default False. If True, use Jaccard containment.
          * best_only: default False. If True, allow optimizations that
            may. May discard matches better than threshold, but first match
            is guaranteed to be best.
          * ignore_abundance: default False. If True, and query signature
            and database support k-mer abundances, ignore those abundances.
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
    def select(self, ksize=None, moltype=None):
        ""

class LinearIndex(Index):
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

    def select(self, ksize=None, moltype=None):
        def select_sigs(siglist, ksize, moltype):
            for ss in siglist:
                if (ksize is None or ss.minhash.ksize == ksize) and \
                   (moltype is None or ss.minhash.moltype == moltype):
                   yield ss

        siglist=select_sigs(self._signatures, ksize, moltype)
        return LinearIndex(siglist, self.filename)
