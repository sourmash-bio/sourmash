"""
Code for searching collections of signatures.
"""
import csv
import numpy as np
from enum import Enum
import numpy as np
from dataclasses import dataclass

from .signature import SourmashSignature, MinHash
from .sketchcomparison import FracMinHashComparison, NumMinHashComparison


def calc_threshold_from_bp(threshold_bp, scaled, query_size):
    """
    Convert threshold_bp (threshold in estimated bp) to
    fraction of query & minimum number of hashes needed.
    """
    threshold = 0.0
    n_threshold_hashes = 0

    if threshold_bp:
        if threshold_bp < 0:
            raise TypeError("threshold_bp must be non-negative")

        # if we have a threshold_bp of N, then that amounts to N/scaled
        # hashes:
        n_threshold_hashes = float(threshold_bp) / scaled

        # that then requires the following containment:
        threshold = n_threshold_hashes / query_size

        # is it too high to ever match?
        if threshold > 1.0:
            raise ValueError("requested threshold_bp is unattainable with this query")
    return threshold, n_threshold_hashes


class SearchType(Enum):
    JACCARD = 1
    CONTAINMENT = 2
    MAX_CONTAINMENT = 3


def make_jaccard_search_query(*,
                              do_containment=False,
                              do_max_containment=False,
                              best_only=False,
                              threshold=None):
    """\
    Make a "flat" search object for Jaccard search & containment.
    """
    if do_containment and do_max_containment:
        raise TypeError("'do_containment' and 'do_max_containment' cannot both be True")

    # configure search - containment? ignore abundance? best only?
    search_cls = JaccardSearch
    if best_only:
        search_cls = JaccardSearchBestOnly

    if do_containment:
        search_obj = search_cls(SearchType.CONTAINMENT, threshold)
    elif do_max_containment:
        search_obj = search_cls(SearchType.MAX_CONTAINMENT, threshold)
    else:
        search_obj = search_cls(SearchType.JACCARD, threshold)

    return search_obj


def make_containment_query(query_mh, threshold_bp, *, best_only=True):
    "Make a search object for containment, with threshold_bp."
    if not query_mh:
        raise ValueError("query is empty!?")

    scaled = query_mh.scaled
    if not scaled:
        raise TypeError("query signature must be calculated with scaled")

    # are we setting a threshold?
    threshold, _ = calc_threshold_from_bp(threshold_bp, scaled, len(query_mh))

    if best_only:
        search_obj = JaccardSearchBestOnly(SearchType.CONTAINMENT,
                                           threshold=threshold)
    else:
        search_obj = JaccardSearch(SearchType.CONTAINMENT,
                                   threshold=threshold)

    return search_obj


class JaccardSearch:
    """
    A class used by Index classes for searching/gathering.
    """
    def __init__(self, search_type, threshold=None):
        "Constructor. Takes type of search, and optional threshold."
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
        self.require_scaled = require_scaled

        if threshold is None:
            threshold = 0
        self.threshold = float(threshold)

    def check_is_compatible(self, sig):
        """
        Is this query compatible with this type of search? Raise TypeError
        if not.
        """
        if self.require_scaled:
            if not sig.minhash.scaled:
                raise TypeError("this search requires a scaled signature")

        if sig.minhash.track_abundance:
            raise TypeError("this search cannot be done with an abund signature")

    def passes(self, score):
        """Return True if this score meets or exceeds the threshold.

        Note: this can be used whenever a score or estimate is available
        (e.g. internal nodes on an SBT). `collect(...)`, below, decides
        whether a particular signature should be collected, and/or can
        update the threshold (used for BestOnly behavior).
        """
        if score and score >= self.threshold:
            return True
        return False

    def collect(self, score, match_sig):
        "Return True if this match should be collected."
        return True

    def score_jaccard(self, query_size, shared_size, subject_size, total_size):
        "Calculate Jaccard similarity."
        if total_size == 0:
            return 0
        return shared_size / total_size

    def score_containment(self, query_size, shared_size, subject_size,
                          total_size):
        "Calculate Jaccard containment."
        if query_size == 0:
            return 0
        return shared_size / query_size

    def score_max_containment(self, query_size, shared_size, subject_size,
                              total_size):
        "Calculate Jaccard max containment."
        min_denom = min(query_size, subject_size)
        if min_denom == 0:
            return 0
        return shared_size / min_denom


class JaccardSearchBestOnly(JaccardSearch):
    "A subclass of JaccardSearch that implements best-only."
    def collect(self, score, match):
        "Raise the threshold to the best match found so far."
        self.threshold = max(self.threshold, score)
        return True

@dataclass
class BaseResult:
    """
    Base class for sourmash search results.
    Since we need some additional info (scaled vs num minhashes) to
    properly initialize a SketchComparison, this class doesn't actually do
    anything other than define some functions needed by *Result classes.
    """
    query: SourmashSignature
    match: SourmashSignature
    filename: str = None
    ignore_abundance: bool = False # optionally ignore abundances
    # need these for scaled result comparisons
    estimate_ani_ci: bool = False
    ani_confidence: float = 0.95
    threshold_bp: int = None
    cmp_scaled: int = None
    write_cols: list = None
    potential_false_negative: bool = False

    def init_result(self):
        self.mh1 = self.query.minhash
        self.mh2 = self.match.minhash

    def build_fracminhashcomparison(self):
        self.cmp = FracMinHashComparison(self.mh1, self.mh2, cmp_scaled=self.cmp_scaled,
                                        threshold_bp=self.threshold_bp,
                                        ignore_abundance=self.ignore_abundance,
                                        estimate_ani_ci=self.estimate_ani_ci,
                                        ani_confidence=self.ani_confidence)
        self.cmp_scaled = self.cmp.cmp_scaled
        self.query_scaled = self.mh1.scaled
        self.match_scaled = self.mh2.scaled
        self.size_may_be_inaccurate = self.cmp.size_may_be_inaccurate

    def build_numminhashcomparison(self, cmp_num=None):
        self.cmp = NumMinHashComparison(self.mh1, self.mh2, cmp_num=cmp_num, ignore_abundance=self.ignore_abundance)
        self.cmp_num = self.cmp.cmp_num
        self.query_num = self.mh1.num
        self.match_num = self.mh2.num
        self.size_may_be_inaccurate = self.cmp.size_may_be_inaccurate

    def get_cmpinfo(self):
        # grab signature /minhash metadata
        # note, with so few columns written for current SearchResult, we could move these to PrefetchResult initialization
        # I've left here for now because I think at some point we want to output more info from search..
        self.ksize = self.mh1.ksize
        self.moltype = self.mh1.moltype
        self.query_name = self.query.name
        self.query_filename = self.query.filename
        self.query_md5 = self.query.md5sum()
        self.match_name = self.match.name
        self.match_filename = self.match.filename
        # sometimes filename is not set in sig (match_filename is None),
        # and `search` is able to pass in the filename.
        if self.filename is None and self.match_filename is not None:
            self.filename = self.match_filename
        self.match_md5 = self.match.md5sum()
        # set these from self.match_*
        self.md5= self.match_md5
        self.name = self.match_name
        # could define in PrefetchResult instead, same reasoning as above
        self.query_abundance = self.mh1.track_abundance
        self.match_abundance = self.mh2.track_abundance
        self.query_n_hashes = len(self.mh1)
        self.match_n_hashes = len(self.mh2)

    @property
    def pass_threshold(self):
        return self.cmp.pass_threshold

    def shorten_md5(self, md5):
        return md5[:8]

    def to_write(self, columns=[]):
        # convert comparison attrs into a dictionary
        # that can be used by csv dictwriter
        info = {k: v for k, v in self.__dict__.items()
                if k in columns and v is not None}
        return info

    def init_dictwriter(self, csv_handle):
        # initialize the csv, return writer (do once)
        w = csv.DictWriter(csv_handle, fieldnames=self.write_cols)
        w.writeheader()
        return w

    def prep_result(self):
        # shorten or modify any columns that need modifying
        self.query_md5 = self.shorten_md5(self.query_md5)

    def write(self, w):
        # write result dictionary using csv dictwriter
        self.prep_result()
        w.writerow(self.to_write(columns=w.fieldnames))

    @property
    def resultdict(self):
        # instead of writing, just return dictionary of what we want to write
        self.prep_result()
        return self.to_write(columns=self.write_cols)


@dataclass
class SearchResult(BaseResult):
    """
    SearchResult class supports 'sourmash search' operations.
    """
    similarity: float = None
    cmp_num: int = None
    searchtype: SearchType = None

    #columns for standard SearchResult output
    search_write_cols = ['similarity', 'md5', 'filename', 'name',  # here we use 'filename'
                         'query_filename', 'query_name', 'query_md5', 'ani']

    ci_cols = ["ani_low", "ani_high"]

    search_write_cols_ci = search_write_cols + ci_cols

    def init_sigcomparison(self):
        self.init_result()
        if any([self.mh1.scaled, self.mh2.scaled]):
            self.build_fracminhashcomparison()
        elif any([self.mh1.num, self.mh2.num]):
            self.build_numminhashcomparison(cmp_num=self.cmp_num)
        self.get_cmpinfo() # grab comparison metadata

    def __post_init__(self):
        self.init_sigcomparison() # build sketch comparison
        self.check_similarity()
        if self.cmp_scaled is not None and self.searchtype is not None:
            self.estimate_search_ani()
        # define columns we want to write
        self.write_cols = self.search_write_cols
        if self.estimate_ani_ci:
            self.write_cols = self.search_write_cols_ci

    def check_similarity(self):
        # for now, require similarity for SearchResult
        # future: consider returning SearchResult *during* search, and passing SearchType in.
        # then allow similarity to be calculated here according to SearchType.
        if self.similarity is None:
            raise ValueError("Error: Must provide 'similarity' for SearchResult.")

    def estimate_search_ani(self):
        #future: could estimate ANI from abund searches if we want (use query containment?)
        if self.cmp_scaled is None:
            raise TypeError("Error: ANI can only be estimated from scaled signatures.")
        if self.searchtype == SearchType.CONTAINMENT:
            self.cmp.estimate_ani_from_mh1_containment_in_mh2()
            self.ani = self.cmp.ani_from_mh1_containment_in_mh2
            if self.estimate_ani_ci:
                self.ani_low = self.cmp.ani_from_mh1_containment_in_mh2_low
                self.ani_high = self.cmp.ani_from_mh1_containment_in_mh2_high
        elif self.searchtype == SearchType.MAX_CONTAINMENT:
            self.cmp.estimate_max_containment_ani()
            self.ani = self.cmp.max_containment_ani
            if self.estimate_ani_ci:
                self.ani_low = self.cmp.max_containment_ani_low
                self.ani_high = self.cmp.max_containment_ani_high
        elif self.searchtype == SearchType.JACCARD:
            self.cmp.estimate_jaccard_ani(jaccard=self.similarity)
            self.jaccard_ani_untrustworthy = self.cmp.jaccard_ani_untrustworthy
            self.ani = self.cmp.jaccard_ani
        # this can be set from any of the above
        self.potential_false_negative = self.cmp.potential_false_negative


@dataclass
class PrefetchResult(BaseResult):
    """
    PrefetchResult class supports 'sourmash prefetch' operations.
    """

    # current prefetch columns
    prefetch_write_cols = ['intersect_bp', 'jaccard', 'max_containment', 'f_query_match',
                           'f_match_query', 'match_filename', 'match_name', # here we use 'match_filename'
                           'match_md5', 'match_bp', 'query_filename', 'query_name',
                           'query_md5', 'query_bp', 'ksize', 'moltype', 'scaled',
                           'query_n_hashes', 'query_abundance', 'query_containment_ani',
                           'match_containment_ani', 'average_containment_ani', 'max_containment_ani',
                           'potential_false_negative'] #'match_abundance'

    ci_cols = ["query_containment_ani_low", "query_containment_ani_high",
                   "match_containment_ani_low", "match_containment_ani_high"]

    prefetch_write_cols_ci = prefetch_write_cols + ci_cols

    def init_sigcomparison(self):
        # shared prefetch/gather initialization
        self.init_result()
        if all([self.mh1.scaled, self.mh2.scaled]):
            self.build_fracminhashcomparison()
        else:
            raise TypeError("Error: prefetch and gather results must be between scaled signatures.")
        self.get_cmpinfo() # grab comparison metadata
        self.intersect_bp = self.cmp.total_unique_intersect_hashes
        self.max_containment = self.cmp.max_containment
        self.query_bp = self.mh1.unique_dataset_hashes
        self.match_bp = self.mh2.unique_dataset_hashes
        self.threshold = self.threshold_bp
        self.estimate_containment_ani()

    def estimate_containment_ani(self):
        self.cmp.estimate_all_containment_ani()
        self.query_containment_ani = self.cmp.ani_from_mh1_containment_in_mh2
        self.match_containment_ani = self.cmp.ani_from_mh2_containment_in_mh1
        self.average_containment_ani = self.cmp.avg_containment_ani
        self.max_containment_ani = self.cmp.max_containment_ani
        self.potential_false_negative = self.cmp.potential_false_negative
        if self.estimate_ani_ci:
            self.handle_ani_ci()

    def handle_ani_ci(self):
        self.query_containment_ani_low = self.cmp.ani_from_mh1_containment_in_mh2_low
        self.query_containment_ani_high = self.cmp.ani_from_mh1_containment_in_mh2_high
        self.match_containment_ani_low = self.cmp.ani_from_mh2_containment_in_mh1_low
        self.match_containment_ani_high = self.cmp.ani_from_mh2_containment_in_mh1_high

    def build_prefetch_result(self):
        # unique prefetch values
        self.jaccard = self.cmp.jaccard
        self.f_query_match = self.cmp.mh2_containment_in_mh1 #db_mh.contained_by(query_mh)
        self.f_match_query = self.cmp.mh1_containment_in_mh2 #query_mh.contained_by(db_mh)
        # set write columns for prefetch result
        self.write_cols = self.prefetch_write_cols
        if self.estimate_ani_ci:
            self.write_cols = self.prefetch_write_cols_ci

    def __post_init__(self):
        self.init_sigcomparison()
        self.build_prefetch_result()

    def prep_prefetch_result(self):
        # explicitly name so we can use this within GatherResult too
        self.scaled = self.cmp_scaled
        # in prefetch, we shorten all md5's
        self.query_md5 = self.shorten_md5(self.query_md5)
        self.md5 = self.shorten_md5(self.md5)
        self.match_md5 = self.shorten_md5(self.match_md5)

    def prep_result(self):
        # overwrite base prep_result
        self.prep_prefetch_result()

    @property
    def prefetchresultdict(self):
        # just return dictionary of what we want to write
        self.prep_prefetch_result()
        return self.to_write(columns=self.write_cols)


@dataclass
class GatherResult(PrefetchResult):
    gather_querymh: MinHash = None
    gather_result_rank: int = None
    orig_query_len: int = None
    orig_query_abunds: list = None
    sum_weighted_found: int = None
    total_weighted_hashes: int = None

    gather_write_cols = ['intersect_bp', 'f_orig_query', 'f_match', 'f_unique_to_query',
                         'f_unique_weighted','average_abund', 'median_abund', 'std_abund', 'filename', # here we use 'filename'
                         'name', 'md5', 'f_match_orig', 'unique_intersect_bp', 'gather_result_rank',
                         'remaining_bp', 'query_filename', 'query_name', 'query_md5', 'query_bp', 'ksize',
                         'moltype', 'scaled', 'query_n_hashes', 'query_abundance', 'query_containment_ani',
                         'match_containment_ani', 'average_containment_ani', 'max_containment_ani',
                         'potential_false_negative',
                         'n_unique_weighted_found', 'sum_weighted_found',
                         'total_weighted_hashes']

    ci_cols = ["query_containment_ani_low", "query_containment_ani_high",
                   "match_containment_ani_low", "match_containment_ani_high"]

    gather_write_cols_ci = gather_write_cols + ci_cols

    def init_gathersketchcomparison(self):
        # compare remaining gather hashes with match. Force at cmp_scaled. Force match flatten(), bc we don't need abunds.
        self.gather_comparison = FracMinHashComparison(self.gather_querymh, self.match.minhash.flatten())

    def check_gatherresult_input(self):
        # check we have what we need:
        if self.cmp_scaled is None:
            raise ValueError("Error: must provide comparison scaled value ('cmp_scaled') for GatherResult")
        if self.gather_querymh is None:
            raise ValueError("Error: must provide current gather sketch (remaining hashes) for GatherResult")
        if self.gather_result_rank is None:
            raise ValueError("Error: must provide 'gather_result_rank' to GatherResult")
        if not self.total_weighted_hashes: # catch total_weighted_hashes = 0 as well
            raise ValueError("Error: must provide sum of all abundances ('total_weighted_hashes') to GatherResult")
        if not self.orig_query_abunds:
            raise ValueError("Error: must provide original query abundances ('orig_query_abunds') to GatherResult")

    def build_gather_result(self):
        # build gather-specific attributes
    
        # the 'query' that is passed into gather is all _matched_ hashes, after subtracting noident_mh
        # this affects estimation of original query information, and requires us to pass in orig_query_len and orig_query_abunds.
        # we also need to overwrite self.query_bp, self.query_n_hashes, and self.query_abundance
        # todo: find a better solution?
        self.query_bp = self.orig_query_len * self.query.minhash.scaled
        self.query_n_hashes = self.orig_query_len

        # calculate intersection with query hashes:
        self.unique_intersect_bp = self.gather_comparison.total_unique_intersect_hashes
    
        # calculate fraction of subject match with orig query
        self.f_match_orig = self.cmp.mh2_containment_in_mh1

        # calculate fractions wrt first denominator - genome size
        self.f_match = self.gather_comparison.mh2_containment_in_mh1 # unique match containment
        self.f_orig_query = len(self.cmp.intersect_mh) / self.orig_query_len
        assert self.gather_comparison.intersect_mh.contained_by(self.gather_comparison.mh1_cmp) == 1.0
    
        # calculate fractions wrt second denominator - metagenome size
        assert self.gather_comparison.intersect_mh.contained_by(self.gather_comparison.mh2_cmp) == 1.0
        self.f_unique_to_query = len(self.gather_comparison.intersect_mh)/self.orig_query_len

        # here, need to make sure to use the mh1_cmp (bc was downsampled to cmp_scaled)
        self.remaining_bp = (self.gather_comparison.mh1_cmp.unique_dataset_hashes - self.gather_comparison.total_unique_intersect_hashes)

        # calculate stats on abundances, if desired.
        self.average_abund, self.median_abund, self.std_abund = None, None, None
        if not self.ignore_abundance:
            self.query_weighted_unique_intersection = self.gather_comparison.weighted_intersection(from_abundD = self.orig_query_abunds)
            self.average_abund = self.query_weighted_unique_intersection.mean_abundance
            self.median_abund = self.query_weighted_unique_intersection.median_abundance
            self.std_abund = self.query_weighted_unique_intersection.std_abundance
            # 'query' will be flattened by default. reset track abundance if we have abunds
            self.query_abundance = self.query_weighted_unique_intersection.track_abundance
            # calculate scores weighted by abundances
            self.n_unique_weighted_found = self.query_weighted_unique_intersection.sum_abundances
            self.f_unique_weighted = self.n_unique_weighted_found / self.total_weighted_hashes
        else:
            self.f_unique_weighted = self.f_unique_to_query
            self.query_abundance = False

    def __post_init__(self):
        self.check_gatherresult_input()
        self.init_sigcomparison() # initialize original sketch vs match sketch comparison (inherited from PrefetchResult)
        self.init_gathersketchcomparison() # initialize remaining gather sketch vs match sketch comparison
        self.build_gather_result() # build gather-specific attributes
        # set write columns for prefetch result
        self.write_cols = self.gather_write_cols
        if self.estimate_ani_ci:
            self.write_cols = self.gather_write_cols_ci

    def prep_gather_result(self):
        # for gather, we only shorten the query_md5
        self.scaled = self.cmp_scaled
        self.query_md5 = self.shorten_md5(self.query_md5)

    def prep_result(self):
        # overwrite base prep_result
        self.prep_gather_result()

    @property
    def gatherresultdict(self):
        # just return dictionary of what we want to write
        self.prep_gather_result()
        return self.to_write(columns=self.write_cols)

    @property
    def prefetchresultdict(self):
        # get all prefetch cols from gatherresult
        prefetch_cols = self.prefetch_write_cols
        if self.estimate_ani_ci:
            prefetch_cols = self.prefetch_write_cols_ci
        self.jaccard = self.cmp.jaccard
        self.f_query_match = self.cmp.mh2_containment_in_mh1 #db_mh.contained_by(query_mh)
        self.f_match_query = self.cmp.mh1_containment_in_mh2 #query_mh.contained_by(db_mh)
        self.prep_prefetch_result()
        return self.to_write(columns=prefetch_cols)


def format_bp(bp):
    "Pretty-print bp information."
    bp = float(bp)
    if bp < 500:
        return '{:.0f} bp '.format(bp)
    elif bp <= 500e3:
        return '{:.1f} kbp'.format(round(bp / 1e3, 1))
    elif bp < 500e6:
        return '{:.1f} Mbp'.format(round(bp / 1e6, 1))
    elif bp < 500e9:
        return '{:.1f} Gbp'.format(round(bp / 1e9, 1))
    return '???'


def search_databases_with_flat_query(query, databases, **kwargs):
    results = []
    found_md5 = set()

    for db in databases:
        search_iter = db.search(query, **kwargs)
        for (score, match, filename) in search_iter:
            md5 = match.md5sum()
            if md5 not in found_md5:
                results.append((score, match, filename))
                found_md5.add(md5)

    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x[0])

    # redefine searchtype and pass in here
    # repetitive/not optimal - would it be better to produce SearchResult from db.search?
    estimate_ani_ci = False
    search_type = SearchType.JACCARD
    if kwargs.get('do_containment'):
        search_type = SearchType.CONTAINMENT
        if kwargs.get('estimate_ani_ci'):
            estimate_ani_ci = True
    elif kwargs.get('do_max_containment'):
        search_type = SearchType.MAX_CONTAINMENT
        if kwargs.get('estimate_ani_ci'):
            estimate_ani_ci = True

    x = []
    for (score, match, filename) in results:
        x.append(SearchResult(query, match,
                              similarity=score,
                              filename = filename,
                              searchtype=search_type,
                              estimate_ani_ci=estimate_ani_ci))
    return x


def search_databases_with_abund_query(query, databases, **kwargs):
    results = []
    found_md5 = set()

    if kwargs.get('do_containment') or kwargs.get('do_max_containment'):
        raise TypeError("containment searches cannot be done with abund sketches")

    for db in databases:
        search_iter = db.search_abund(query, **kwargs) # could return SearchResult here instead of tuple?
        for (score, match, filename) in search_iter:
            md5 = match.md5sum()
            if md5 not in found_md5:
                results.append((score, match, filename))
                found_md5.add(md5)

    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x[0])

    x = []
    for (score, match, filename) in results:
        x.append(SearchResult(query, match,
                               similarity=score,
                               filename = filename))
    return x

###
### gather code
###

def _find_best(counters, query, threshold_bp):
    """
    Search for the best containment, return precisely one match.
    """
    best_result = None
    best_intersect_mh = None

    # find the best score across multiple counters, without consuming
    for counter in counters:
        result = counter.peek(query.minhash, threshold_bp=threshold_bp)
        if result:
            (sr, intersect_mh) = result

            if best_result is None or sr.score > best_result.score:
                best_result = sr
                best_intersect_mh = intersect_mh

    if best_result:
        # remove the best result from each counter
        for counter in counters:
            counter.consume(best_intersect_mh)

        # and done!
        return best_result, best_intersect_mh
    return None, None


class GatherDatabases:
    "Iterator object for doing gather/min-set-cov."

    def __init__(self, query, counters, *,
                 threshold_bp=0, ignore_abundance=False, noident_mh=None, ident_mh=None, estimate_ani_ci=False):
        # track original query information for later usage?
        track_abundance = query.minhash.track_abundance and not ignore_abundance
        self.orig_query = query
        self.orig_query_bp = query.minhash.unique_dataset_hashes
        self.orig_query_filename = query.filename
        self.orig_query_name = query.name
        self.orig_query_md5 = query.md5sum()[:8]

        # do we pay attention to abundances?
        query_mh = query.minhash
        query_hashes = query_mh.hashes
        orig_query_abunds = { k: 1 for k in query_hashes }
        if track_abundance:
            orig_query_abunds = query_hashes

        # adjust for not found...
        if noident_mh is None:  # create empty
            noident_mh = query_mh.copy_and_clear()
        self.noident_mh = noident_mh.to_frozen()

        if ident_mh is None:
            query_mh = query_mh.to_mutable()
            query_mh.remove_many(noident_mh)
        else:
            query_mh = ident_mh.to_mutable()

        orig_query_mh = query_mh.flatten()

        # query.minhash will be assigned to repeatedly in gather; make mutable.
        query = query.to_mutable()
        query.minhash = orig_query_mh

        cmp_scaled = query.minhash.scaled    # initialize with resolution of query

        self.result_n = 0
        self.query = query
        self.counters = counters
        self.threshold_bp = threshold_bp

        self.track_abundance = track_abundance
        self.orig_query_mh = orig_query_mh
        self.orig_query_abunds = orig_query_abunds

        self.cmp_scaled = 0     # initialize with something very low!
        self._update_scaled(cmp_scaled)

        self.estimate_ani_ci = estimate_ani_ci # by default, do not report ANI confidence intervals

    def _update_scaled(self, scaled):
        max_scaled = max(self.cmp_scaled, scaled)
        if self.cmp_scaled != max_scaled:
            self.cmp_scaled = max_scaled

            # CTB note: this can be expensive
            self.orig_query_mh = self.orig_query_mh.downsample(scaled=scaled)
            self.noident_mh = self.noident_mh.downsample(scaled=scaled)

            # NOTE: orig_query_abunds can be used w/o downsampling
            orig_query_abunds = self.orig_query_abunds
            self.noident_query_sum_abunds = sum(( orig_query_abunds[k] \
                                                  for k in self.noident_mh.hashes ))
            self.sum_abunds = sum(( orig_query_abunds[k] \
                                    for k in self.orig_query_mh.hashes ))
            self.sum_abunds += self.noident_query_sum_abunds

        if max_scaled != scaled:
            return max_scaled
        return max_scaled

    @property
    def scaled(self):
        return self.cmp_scaled

    def __iter__(self):
        return self

    def __next__(self):
        query = self.query
        if not self.query.minhash:
            raise StopIteration

        # may be changed:
        counters = self.counters
        cmp_scaled = self.cmp_scaled

        # will not be changed::
        track_abundance = self.track_abundance
        threshold_bp = self.threshold_bp
        orig_query_abunds = self.orig_query_abunds

        # find the best match!
        best_result, intersect_mh = _find_best(counters, query, threshold_bp)

        if not best_result:          # no matches at all for this cutoff!
            raise StopIteration

        best_match = best_result.signature
        filename = best_result.location

        # Is the best match computed with scaled? Die if not.
        match_scaled = best_match.minhash.scaled
        assert match_scaled

        # pick the highest scaled / lowest resolution.
        scaled = self._update_scaled(match_scaled)
        # CTB note: this means that if a high scaled/low res signature is
        # found early on, resolution will be low from then on.

        # retrieve various saved things, after potential downsampling
        orig_query_mh = self.orig_query_mh
        sum_abunds = self.sum_abunds
        noident_mh = self.noident_mh
        orig_query_len = len(orig_query_mh) + len(noident_mh)

        # eliminate hashes under this new resolution.
        query_mh = query.minhash.downsample(scaled=scaled)
        found_mh = best_match.minhash.downsample(scaled=scaled).flatten()

        # construct a new query, subtracting hashes found in previous one.
        new_query_mh = query_mh.to_mutable()
        new_query_mh.remove_many(found_mh)
        new_query = SourmashSignature(new_query_mh)

        # compute weighted information for remaining query hashes
        query_hashes = set(query_mh.hashes) - set(found_mh.hashes)
        n_weighted_missed = sum((orig_query_abunds[k] for k in query_hashes))
        n_weighted_missed += self.noident_query_sum_abunds
        sum_weighted_found = sum_abunds - n_weighted_missed

        # build a GatherResult
        result = GatherResult(self.orig_query, best_match,
                              cmp_scaled=scaled,
                              filename=filename,
                              gather_result_rank=self.result_n,
                              gather_querymh=query.minhash,
                              ignore_abundance= not track_abundance,
                              threshold_bp=threshold_bp,
                              orig_query_len=orig_query_len,
                              orig_query_abunds = self.orig_query_abunds,
                              estimate_ani_ci=self.estimate_ani_ci,
                              sum_weighted_found=sum_weighted_found,
                              total_weighted_hashes=sum_abunds,
                              )

        self.result_n += 1
        self.query = new_query
        self.orig_query_mh = orig_query_mh

        return result


###
### prefetch code
###

def prefetch_database(query, database, threshold_bp, *, estimate_ani_ci=False):
    """
    Find all matches to `query_mh` >= `threshold_bp` in `database`.
    """
    scaled = query.minhash.scaled
    assert scaled
    # iterate over all signatures in database, find matches
    for result in database.prefetch(query, threshold_bp): # future: could return PrefetchResult directly here
        result = PrefetchResult(query, result.signature, threshold_bp=threshold_bp, estimate_ani_ci=estimate_ani_ci)
        assert result.pass_threshold
        yield result
