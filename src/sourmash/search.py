"""
Code for searching collections of signatures.
"""
from collections import namedtuple
from enum import Enum
from multiprocessing.sharedctypes import Value
from re import I
import numpy as np
from dataclasses import dataclass

from .signature import SourmashSignature, MinHash
from .sketchcomparison import FracMinHashComparison, NumMinHashComparison


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


def make_gather_query(query_mh, threshold_bp, *, best_only=True):
    "Make a search object for gather."
    if not query_mh:
        raise ValueError("query is empty!?")

    scaled = query_mh.scaled
    if not scaled:
        raise TypeError("query signature must be calculated with scaled")

    # are we setting a threshold?
    threshold = 0
    if threshold_bp:
        if threshold_bp < 0:
            raise TypeError("threshold_bp must be non-negative")

        # if we have a threshold_bp of N, then that amounts to N/scaled
        # hashes:
        n_threshold_hashes = threshold_bp / scaled

        # that then requires the following containment:
        threshold = n_threshold_hashes / len(query_mh)

        # is it too high to ever match? if so, exit.
        if threshold > 1.0:
            raise ValueError("requested threshold_bp is unattainable with this query")

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

    def init_result(self):
        self.mh1 = self.query.minhash
        self.mh2 = self.match.minhash

    def build_fracminhashcomparison(self, cmp_scaled=None, threshold_bp=None):
        self.cmp = FracMinHashComparison(self.mh1, self.mh2, cmp_scaled=cmp_scaled, threshold_bp=threshold_bp, ignore_abundance=self.ignore_abundance)
        self.scaled = self.cmp.cmp_scaled
        self.query_scaled = self.mh1.scaled
        self.match_scaled = self.mh2.scaled

    def build_numminhashcomparison(self, cmp_num=None):
        self.cmp = NumMinHashComparison(self.mh1, self.mh2, cmp_num=cmp_num, ignore_abundance=self.ignore_abundance)
        self.num = self.cmp.cmp_num
        self.query_num = self.mh1.num
        self.match_num = self.mh2.num

    def get_cmpinfo(self):
        # grab signature /minhash metadata
        self.ksize = self.mh1.ksize
        self.moltype = self.mh1.moltype
        self.query_name = self.query.name
        self.query_filename = self.query.filename
        self.query_md5 = self.query.md5sum()
        self.match_name = self.match.name
        self.match_filename = self.match.filename
        # sometimes filename is not set in sig (match_filename is None),
        # and `search` is able to pass in the filename...
        if self.filename is None and self.match_filename is not None:
            self.filename = self.match_filename
        self.match_md5 = self.match.md5sum()
        # set these from self.match_*
        self.md5= self.match_md5
        self.name = self.match_name
        # do we actually need these here? Def need for prefetch, but maybe define there?
        self.query_abundance = self.mh1.track_abundance
        self.match_abundance = self.mh2.track_abundance
        self.query_n_hashes = len(self.mh1.hashes)
        self.match_n_hashes = len(self.mh2.hashes)

    @property
    def pass_threshold(self):
        return self.cmp.pass_threshold

    def shorten_md5(self, md5):
        return md5[:8]

    def to_write(self, columns=[]):
        info = {k: v for k, v in self.__dict__.items()
                if k in columns and v is not None}
        return info


@dataclass
class SearchResult(BaseResult):
    """
    SearchResult class supports 'sourmash search' operations.
    """
    similarity: float = None
    cmp_scaled: int = None
    cmp_num: int = None
    threshold_bp: int = None

    #columns for standard SearchResult output
    search_write_cols = ['similarity', 'md5', 'filename', 'name',  # here we use 'filename'
                         'query_filename', 'query_name', 'query_md5']

    def init_sigcomparison(self):
        self.init_result()
        if any([self.mh1.scaled, self.mh2.scaled]):
            self.build_fracminhashcomparison(cmp_scaled = self.cmp_scaled, threshold_bp = self.threshold_bp)
        elif any([not self.mh1.num, not self.mh2.num]):
            self.build_numminhashcomparison(cmp_num=self.cmp_num)
        self.get_cmpinfo() # grab comparison metadata

    def __post_init__(self):
        self.init_sigcomparison() # build sketch comparison
        self.check_similarity() # set similarity (if not passed in)

    def check_similarity(self):
        # Require similarity for SearchResult
        if self.similarity is None:
            raise ValueError("Error: Must provide 'similarity' for SearchResult.")
            #OR, if don't pass in similarity, return jaccard?
            #if self.cmp.cmp_scaled is not None:
            #    self.similarity = self.cmp.mh1_containment
            #else:
            #    self.similarity = self.cmp.jaccard

    @property
    def writedict(self):
        self.query_md5 = self.shorten_md5(self.query_md5)
        return self.to_write(columns=self.search_write_cols)



@dataclass
class PrefetchResult(BaseResult):
    """
    PrefetchResult class supports 'sourmash prefetch' operations.
    """
    cmp_scaled: int = None
    threshold_bp: int = None

    # current prefetch columns
    prefetch_write_cols = ['intersect_bp', 'jaccard', 'max_containment', 'f_query_match',
                           'f_match_query', 'match_filename', 'match_name', # here we use 'match_filename'
                           'match_md5', 'match_bp', 'query_filename', 'query_name',
                           'query_md5', 'query_bp', 'ksize', 'moltype', 'scaled',
                           'query_n_hashes', 'query_abundance'] #, 'match_abundance'

    def init_sigcomparison(self):
        # shared prefetch/gather initialization
        self.init_result()
        if all([self.mh1.scaled, self.mh2.scaled]):
            self.build_fracminhashcomparison(cmp_scaled = self.cmp_scaled, threshold_bp = self.threshold_bp)
        else:
            raise TypeError("Error: prefetch and gather results must be between scaled signatures.")
        self.get_cmpinfo() # grab comparison metadata
        self.intersect_bp = self.cmp.intersect_bp
        self.max_containment = self.cmp.max_containment
        self.query_bp = self.mh1.covered_bp
        self.match_bp = self.mh2.covered_bp

    def build_prefetch_result(self):
        # unique prefetch values
        self.jaccard = self.cmp.jaccard
        self.f_query_match = self.cmp.mh2_containment #db_mh.contained_by(query_mh)
        self.f_match_query = self.cmp.mh1_containment #query_mh.contained_by(db_mh)

    def __post_init__(self):
        self.init_sigcomparison()
        self.build_prefetch_result()

    def prefetchresultdict(self):
        # in prefetch, we shorten all md5's
        self.query_md5 = self.shorten_md5(self.query_md5)
        self.md5 = self.shorten_md5(self.md5)
        self.match_md5 = self.shorten_md5(self.match_md5)
        return self.to_write(columns=self.prefetch_write_cols)

    @property
    def writedict(self):
        return self.prefetchresultdict()



@dataclass
class GatherResult(PrefetchResult):
    gather_querymh: MinHash = None
    gather_result_rank: int = None
    total_abund: int = None
   # orig_query_len includes len(noident_mh), which has been subtracted out of query, right? This subtraction is an issue for query_covered bp,etc!
   # can we set some original query info so we can get the right vals to output?
    orig_query_len: int = None

    # pass in for now bc I'm not calculating correctly...
    remaining_bp: int = None
    average_abund: int = None
    median_abund: int = None
    std_abund: int = None
    f_unique_weighted: int =None

    gather_write_cols = ['intersect_bp', 'f_orig_query', 'f_match', 'f_unique_to_query',
                         'f_unique_weighted','average_abund', 'median_abund', 'std_abund', 'filename', # here we use 'filename'
                         'name', 'md5', 'f_match_orig', 'unique_intersect_bp', 'gather_result_rank',
                         'remaining_bp', 'query_filename', 'query_name', 'query_md5', 'query_bp', 'ksize',
                         'moltype', 'scaled', 'query_n_hashes', 'query_abundance']

    def init_gathersketchcomparison(self):
        # compare remaining gather hashes with match. Force at cmp_scaled. Do we need for force match flatten()?
        self.gather_comparison = FracMinHashComparison(self.gather_querymh, self.match.minhash.flatten(), cmp_scaled=self.cmp_scaled, threshold_bp=self.threshold_bp)
        #self.gather_comparison = FracMinHashComparison(self.gather_querymh, self.match.minhash, cmp_scaled=self.cmp_scaled, threshold_bp=self.threshold_bp)

    def check_gatherresult_input(self):
        # check we have what we need:
        if self.cmp_scaled is None:
            raise ValueError("Error: must provide comparison scaled value ('cmp_scaled') for GatherResult")
        if self.gather_querymh is None:
            raise ValueError("Error: must provide current gather sketch (remaining hashes) for GatherResult")
        if self.gather_result_rank is None:
            raise ValueError("Error: must provide 'gather_result_rank' to GatherResult")
        if not self.total_abund: # catch total_abund = 0 as well
            raise ValueError("Error: must provide sum of all abundances ('total_abund') to GatherResult")

    def build_gather_result(self):
        # build gather specific attributes
        self.f_match_orig = self.cmp.mh2_containment #original match containment
        self.f_match = self.gather_comparison.mh2_containment # unique match containment
        #intersect bp of remaining, unaccounted for hashes with the database match
        self.unique_intersect_bp = self.gather_comparison.intersect_bp
        self.f_orig_query = len(self.cmp.intersect_mh) / self.orig_query_len
        self.f_unique_to_query = len(self.gather_comparison.intersect_mh)/self.orig_query_len

    #def build_gather_result_incorrect(self):
    #    # get current query-weighted intersection
    #    self.query_weighted_unique_intersection = self.gather_comparison.mh1_weighted_intersection
    #    # calculate scores weighted by abundances
    #    #f_unique_weighted = sum((orig_query_abunds[k] for k in intersect_mh.hashes ))
    #    #f_unique_weighted /= sum_abunds
    #    if self.query_abundance:
    #        self.f_unique_weighted =  float(self.query_weighted_unique_intersection.sum_abundances) / self.total_abund
    #    else:
    #        self.f_unique_weighted = self.f_unique_to_query # is this right?
    #    self.remaining_bp = self.query_bp - self.gather_comparison.intersect_bp
    #    #self.remaining_bp = self.gather_comparison.mh1_bp - self.gather_comparison.intersect_bp
    #    self.average_abund = self.query_weighted_unique_intersection.mean_abundance
    #    self.median_abund = self.query_weighted_unique_intersection.median_abundance
    #    self.std_abund = self.query_weighted_unique_intersection.std_abundance

    def __post_init__(self):
        self.check_gatherresult_input()
        self.init_sigcomparison() # initialize original sketch vs match sketch comparison (inherited from PrefetchResult)
        self.init_gathersketchcomparison() # initialize remaining gather sketch vs match sketch comparison
        self.build_gather_result() # build gather-specific attributes

    def gatherresultdict(self):
        # for gather, we only shorten the query_md5
        self.query_md5 = self.shorten_md5(self.query_md5)
        #self.md5 = self.shorten_md5(self.md5)
        #self.match_md5 = self.shorten_md5(self.match_md5)
        return self.to_write(columns=self.gather_write_cols)

    @property
    def writedict(self):
        return self.gatherresultdict()

    # we can write prefetch results from a GatherResult
    @property
    def prefetchwritedict(self):
        self.build_prefetch_result()
        return self.prefetchresultdict()



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

    x = []
    for (score, match, filename) in results:
        x.append(SearchResult(query, match,
                               similarity=score,
                               filename = filename))
    return x


def search_databases_with_abund_query(query, databases, **kwargs):
    results = []
    found_md5 = set()

    if kwargs.get('do_containment') or kwargs.get('do_max_containment'):
        raise TypeError("containment searches cannot be done with abund sketches")

    for db in databases:
        search_iter = db.search_abund(query, **kwargs)
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
                               similarity=score,  # this is actually cosine sim (abund). do we want to specify this in SearchResult somehow?
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
        result = counter.peek(query.minhash, threshold_bp)
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
                 threshold_bp=0, ignore_abundance=False, noident_mh=None):
        # track original query information for later usage?
        track_abundance = query.minhash.track_abundance and not ignore_abundance
        self.orig_query = query
        self.orig_query_bp = len(query.minhash) * query.minhash.scaled
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

        query_mh = query_mh.to_mutable()
        query_mh.remove_many(noident_mh)

        orig_query_mh = query_mh.flatten()
        query.minhash = orig_query_mh.to_mutable()

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

        # calculate intersection with query hashes:
        unique_intersect_bp = scaled * len(intersect_mh)
        intersect_orig_mh = orig_query_mh & found_mh
        intersect_bp = scaled * len(intersect_orig_mh)

        # calculate fractions wrt first denominator - genome size
        assert intersect_mh.contained_by(found_mh) == 1.0  ## TODO: add this check to GatherResult?
        f_match = len(intersect_mh) / len(found_mh)
        f_orig_query = len(intersect_orig_mh) / orig_query_len

        # calculate fractions wrt second denominator - metagenome size
        assert intersect_mh.contained_by(orig_query_mh) == 1.0
        f_unique_to_query = len(intersect_mh) / orig_query_len

        # calculate fraction of subject match with orig query
        f_match_orig = found_mh.contained_by(orig_query_mh)

        # calculate scores weighted by abundances
        f_unique_weighted = sum((orig_query_abunds[k] for k in intersect_mh.hashes ))
        f_unique_weighted /= sum_abunds

        # calculate stats on abundances, if desired.
        average_abund, median_abund, std_abund = None, None, None
        if track_abundance:
            intersect_abunds = (orig_query_abunds[k] for k in intersect_mh.hashes )
            intersect_abunds = list(intersect_abunds)

            average_abund = np.mean(intersect_abunds)
            median_abund = np.median(intersect_abunds)
            std_abund = np.std(intersect_abunds)

        # construct a new query, subtracting hashes found in previous one.
        new_query_mh = query_mh.to_mutable()
        new_query_mh.remove_many(found_mh)
        new_query = SourmashSignature(new_query_mh)

        remaining_bp = scaled * len(new_query_mh)
        #remaining_bp = new_query_mh.bp

        # compute weighted_missed for remaining query hashes
        query_hashes = set(query_mh.hashes) - set(found_mh.hashes)
        weighted_missed = sum((orig_query_abunds[k] for k in query_hashes))
        weighted_missed += self.noident_query_sum_abunds
        weighted_missed /= sum_abunds

        # build a GatherResult
        result = GatherResult(self.orig_query, best_match,
                              cmp_scaled=scaled,
                              filename=filename,
                              gather_result_rank=self.result_n,
                              total_abund= sum_abunds,
                              gather_querymh=query.minhash,
                              ignore_abundance= not track_abundance,
                              threshold_bp=threshold_bp,
                              orig_query_len=orig_query_len,
                              f_unique_weighted=f_unique_weighted,
                              average_abund=average_abund,
                              median_abund=median_abund,
                              std_abund=std_abund,
                              remaining_bp=remaining_bp,
                              )

        # temp dev: make sure these vals are correct
        assert result.remaining_bp == remaining_bp
        assert result.f_unique_weighted == f_unique_weighted
        assert result.average_abund == average_abund
        assert result.median_abund == median_abund
        assert result.std_abund == std_abund


        self.result_n += 1
        self.query = new_query
        self.orig_query_mh = orig_query_mh

        return result, weighted_missed


###
### prefetch code
###

def calculate_prefetch_info(query, match, scaled, threshold_bp):
#def calculate_prefetch_info(query, match, threshold_bp):
    """
    For a single query and match, calculate all search info and return a PrefetchResult.
    """
    # base intersections on downsampled minhashes
#    query_mh = query.minhash

   # scaled = max(scaled, match.minhash.scaled)
 #   query_mh = query_mh.downsample(scaled=scaled)
 #   db_mh = match.minhash.flatten().downsample(scaled=scaled)

    # calculate db match intersection with query hashes:
#    intersect_mh = query_mh & db_mh
#    threshold = threshold_bp / scaled
#    assert len(intersect_mh) >= threshold

    #f_query_match = db_mh.contained_by(query_mh)
    #f_match_query = query_mh.contained_by(db_mh)
    #max_containment = max(f_query_match, f_match_query)

    # build a PrefetchResult
    result = PrefetchResult(query,match, threshold_bp=threshold_bp)
    assert result.pass_threshold

    return result


def prefetch_database(query, database, threshold_bp):
    """
    Find all matches to `query_mh` >= `threshold_bp` in `database`.
    """
    scaled = query.minhash.scaled
    assert scaled

    # iterate over all signatures in database, find matches
    for result in database.prefetch(query, threshold_bp):
        #result = calculate_prefetch_info(query, result.signature, threshold_bp)
        result = calculate_prefetch_info(query, result.signature, scaled, threshold_bp)
        #result = PrefetchResult(query=query, match=result.signature, search_scaled=scaled, threshold_bp=threshold_bp)
        yield result
