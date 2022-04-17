"""
Code for searching collections of signatures.
"""
from collections import namedtuple
from enum import Enum
from multiprocessing.sharedctypes import Value
import numpy as np
from dataclasses import dataclass

from .signature import SourmashSignature, MinHash


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

def shorten_md5(md5):
    return md5[:8]



@dataclass
class FracMinHashComparison:
    """Class for standard comparison between two scaled minhashes"""
    mh1: MinHash
    mh2: MinHash
    cmp_scaled: int = None # scaled value for this comparison (defaults to maximum scaled between the two sigs)
    flatten: bool = False # optionally ignore abundances
    threshold_bp: int = 0 # default 0 bp threshold

    def _downsample_to_cmp_scaled(self, use_abundances=True):
        """
        Downsample and/or flatten minhashes for comparison
        """
        if self.cmp_scaled is None: # record the scaled we're doing this comparison on 
            self.cmp_scaled = max(self.mh1.scaled, self.mh2.scaled)
        if use_abundances:
            self.mh1_cmp = self.mh1.downsample(scaled=self.cmp_scaled)
            self.mh2_cmp = self.mh2.downsample(scaled=self.cmp_scaled)
        else:
            self.mh1_cmp = self.mh1.flatten().downsample(scaled=self.cmp_scaled)
            self.mh2_cmp = self.mh2.flatten().downsample(scaled=self.cmp_scaled)
 
    @property
    def pass_threshold(self):
        return self.intersect_bp >= self.threshold_bp

    @property
    def intersect_mh(self):
        # flatten and intersect
        return self.mh1_cmp.flatten().intersection(self.mh2_cmp.flatten())

    @property
    def intersect_bp(self):
        return len(self.intersect_mh) * self.cmp_scaled

    @property
    def mh1_weighted_intersection(self):
         # map mh1 hash abundances to all intersection hashes.
        if not self.mh1.track_abundance:
            return self.intersect_mh # flattened should have abund 1 for every hash, right?
        else:
            intersect_w_abunds = self.intersect_mh.inflate(self.mh1)
            return intersect_w_abunds

    @property
    def mh2_weighted_intersection(self):
         # map mh2 hash abundances to all intersection hashes.
        if not self.mh2.track_abundance:
            return  self.intersect_mh # flattened should have abund 1 for every hash, right?
        else:
            intersect_w_abunds = self.intersect_mh.inflate(self.mh2)
            return intersect_w_abunds

    # these four work on a single minhash-- add to MinHash instead?
    def sum_mh_abundances(self, mh):
        return sum(v for v in mh.hashes.values())

    ### WHAT DO WE ACTUALLY WANT TO RETURN HERE? Both None and 0 cause issues in different tests...
    def mean_mh_abundances(self, mh):
#        return np.mean((v for v in mh.hashes.values()))
        if not mh.track_abundance:
        #    raise ValueError("Error: Cannot get mean abundance for non-abundance signatures.")
            return "" 
            #return 1.0
        return np.mean(mh.hashes.values())

    def median_mh_abundances(self, mh):
        if not mh.track_abundance:
            #return 1.0
            return ""
        return np.median(mh.hashes.values())
#        return np.median(v for v in mh.hashes.values())

    def std_mh_abundances(self, mh):
        if not mh.track_abundance:
            #return 0.0
            return ""
        return np.std(mh.hashes.values())
#        return np.std(v for v in mh.hashes.values())

    @property
    def jaccard(self):
        return self.mh1_cmp.jaccard(self.mh2_cmp)

    @property
    def angular_similarity(self):
        if not self.mh1_cmp.track_abundance and self.mh2_cmp.track_abundance:
            raise ValueError("Error: Angular (cosine) similarity requires both sketches to track hash abundance.")
        return self.mh1_cmp.angular_similarity(self.mh2_cmp)

    @property
    def cosine_similarity(self):
        return self.angular_similarity

    @property
    def mh1_containment(self):
        return self.mh1_cmp.contained_by(self.mh2_cmp)

    @property
    def mh2_containment(self):
        return self.mh2_cmp.contained_by(self.mh1_cmp)

    @property
    def max_containment(self):
        return self.mh1_cmp.max_containment(self.mh2_cmp)

    @property
    def avg_containment(self):
        return np.mean(self.mh1_containment, self.mh2_containment)

    def check_sketch_compatibility(self):
        # do we need this check + error? Minhash functions should already complain appropriately...
        if not self.mh1.scaled and self.mh2.scaled:
            raise ValueError("Error: Both sketches must be 'scaled'.")
        k1 = self.mh1.ksize
        k2 = self.mh2.ksize
        if k1 != k2:
            raise ValueError(f"Error: Invalid Comparison, ksizes: {k1}, {k2}). Must compare sketches of the same ksize.")
        m1 = self.mh2.moltype
        m2= self.mh2.moltype
        if m1 != m2:
            raise ValueError(f"Error: Invalid Comparison, moltypes: {m1}, {m2}). Must compare sketches of the same moltype.")
        # switch to the following? no,bc checks scaled too
        #if not self.mh1.is_compatible(self.mh2):
        #    raise ValueError("Error: Cannot compare incompatible sketches.")

    def init_sketchcomparison(self):
        "Initialize ScaledComparison using values from provided FracMinHashes"
        self.check_sketch_compatibility()
        self.ksize = self.mh1.ksize
        self.moltype= self.mh1.moltype
        # get original scaled values -- maybe just use .scaled? no need to store separately...
        self.mh1_scaled = self.mh1.scaled
        self.mh2_scaled = self.mh2.scaled
        # handle abundance
        self.mh1_abundance = self.mh1.track_abundance
        self.mh2_abundance = self.mh2.track_abundance
        # for these, do we want the originals, or the cmp_scaled versions?? (or both?)
        # TODO: add to MinHash?, better there than here.
        self.mh1_n_hashes = len(self.mh1)
        self.mh2_n_hashes = len(self.mh2)
        self.mh1_bp = self.mh1_n_hashes * self.mh1_scaled
        self.mh2_bp = self.mh2_n_hashes * self.mh2_scaled
        # does this enable all abund options we want?
        self.compare_abundances = all([self.mh1_abundance, self.mh2_abundance]) # only use abunds if both have them
        if self.flatten:
            self.compare_abundances = False # force flatten
        self._downsample_to_cmp_scaled(use_abundances=self.compare_abundances)
        # build comparison minhashes at cmp_scaled
        #sevlf._downsample_andor_flatten(use_abundances=self.compare_abundances)
        #self.intersect_bp = len(self.intersect_mh) * self.cmp_scaled
    
    def __post_init__(self):
        self.init_sketchcomparison()


@dataclass
class SignatureComparison:
    """Base class for sourmash search results."""
    query: SourmashSignature
    match: SourmashSignature
    cmp_scaled: int = None # allow ability to pass in a specific scaled value for comparison
    threshold_bp: int = 0 # optional bp threshold for comparison

    write_cols = ['similarity', 'md5', 'filename', 'name',
                  'query_filename', 'query_name', 'query_md5']

    def get_signature_info(self):
        # get info about each signature
        self.match_name = self.match.name
        self.query_name = self.query.name
      #  self.match_filename = self.match.filename
#        self.query_filename = self.query.filename
        self.match_md5 = self.match.md5sum()
        self.query_md5 = self.query.md5sum()

    def get_fracminhash_comparison_info(self):
        self.sketch_comparison = FracMinHashComparison(self.query.minhash, self.match.minhash, cmp_scaled=self.cmp_scaled, threshold_bp=self.threshold_bp)
        # populate fracminhashcomparison info into this class. kinda repetitive, better way to import all sketchcomparison info?
        self.ksize = self.sketch_comparison.ksize
        self.moltype = self.sketch_comparison.moltype
        self.cmp_scaled = self.sketch_comparison.cmp_scaled
        self.query_scaled = self.sketch_comparison.mh1.scaled
        self.match_scaled = self.sketch_comparison.mh2.scaled
        self.query_abundance = self.sketch_comparison.mh1.track_abundance
        self.match_abundance = self.sketch_comparison.mh2.track_abundance
        self.compare_abundances = self.sketch_comparison.compare_abundances
        
        # actually, maybe add these to the minhash implementation to calculate as needed?
        self.query_n_hashes = self.sketch_comparison.mh1_n_hashes
        self.match_n_hashes = self.sketch_comparison.mh2_n_hashes 
        self.query_bp = self.sketch_comparison.mh1_bp
        self.match_bp = self.sketch_comparison.mh2_bp

        # automatically make some comparisons
        self.intersect_bp = self.sketch_comparison.intersect_bp

    def init_sigcomparison(self):
        ## get minhash sketch info / comparison
        self.get_signature_info()
        self.get_fracminhash_comparison_info()
        #deprecate these at some point, in favor of match_md5, etc?
        self.md5 = self.match_md5
        self.name = self.match_name

    def __post_init__(self):
        self.init_sigcomparison()

    @property
    def match_filename(self):
        return self.match.filename

    @property
    def query_filename(self):
        return self.query.filename

    #@property
    #def do_intersect_bp(self):
    #    return self.sketch_comparison.intersect_bp

    @property
    def query_weighted_intersection(self):
        return self.sketch_comparison.mh1_weighted_intersection

    @property
    def match_weighted_intersection(self):
        return self.sketch_comparison.mh2_weighted_intersection

    @property
    def query_weighted_intersection_sum_abunds(self):
        return self.sketch_comparison.sum_mh_abundances(self.query_weighted_intersection)

    @property
    def query_weighted_intersection_mean_abunds(self):
        return self.sketch_comparison.mean_mh_abundances(self.query_weighted_intersection)

    @property
    def query_weighted_intersection_median_abunds(self):
        return self.sketch_comparison.median_mh_abundances(self.query_weighted_intersection)

    @property
    def query_weighted_intersection_std_abunds(self):
        return self.sketch_comparison.std_mh_abundances(self.query_weighted_intersection)

    @property
    def jaccard(self):
        return self.sketch_comparison.jaccard

    @property
    def query_containment(self):
        return self.sketch_comparison.mh1_containment

    @property
    def match_containment(self):
        return self.sketch_comparison.mh2_containment

    @property
    def max_containment(self):
        return self.sketch_comparison.max_containment

    @property
    def avg_containment(self):
        return self.sketch_comparison.avg_containment

    @property
    def angular_similarity(self):
        return self.sketch_comparison.angular_similarity

    @property
    def pass_threshold(self):
        return self.sketch_comparison.pass_threshold
    
    def to_write(self, columns=write_cols):
        info = {k: v for k, v in self.__dict__.items()
                if k in columns and v is not None}
        return info


@dataclass
class SearchResult(SignatureComparison):
    """Base class for sourmash search results."""
    similarity: float = None
    filename: str = None
#    comparison_type: str = None # could specify jaccard vs containment vs angular sim?

    #columns for standard SearchResult output
    search_write_cols = ['similarity', 'md5', 'filename', 'name',
                         'query_filename', 'query_name', 'query_md5']

    def __post_init__(self):
        # sometimes filename is not set in sig (match_filename is None), 
        # and `search` is able to pass in the filename...
        if self.filename is None and self.match_filename is not None:
            self.filename = self.match_filename
        self.init_sigcomparison()

    def searchresultdict(self):
        self.query_md5 = shorten_md5(self.query_md5)
        return self.to_write(columns=self.search_write_cols)

    @property
    def writedict(self):
        return self.searchresultdict()

@dataclass
class PrefetchResult(SignatureComparison):

    prefetch_write_cols = ['intersect_bp', 'jaccard', 'max_containment', 'f_query_match',
                           'f_match_query', 'match', 'match_filename', 'match_name',
                           'match_md5', 'match_bp', 'query', 'query_filename', 'query_name',
                           'query_md5', 'query_bp', 'ksize', 'moltype', 'num', 'scaled',
                           'query_n_hashes', 'query_abundance']

    def build_prefetch_result(self):
    ## ARE THESE THE RIGHT DIRECTION???
    #f_query_match = db_mh.contained_by(query_mh)
    #f_match_query = query_mh.contained_by(db_mh)
       self.f_match_query = self.query_containment
       self.f_query_match = self.match_containment


    def __post_init__(self):
        # do we ever pass filename in prefetch?
        self.init_sigcomparison()
        self.build_prefetch_result()

    def prefetchresultdict(self):
        self.query_md5 = shorten_md5(self.query_md5)
        self.md5 = shorten_md5(self.md5)
        self.match_md5 = shorten_md5(self.match_md5)
        return self.to_write(columns=self.prefetch_write_cols)
    
    @property
    def writedict(self):
        return self.prefetchresultdict()


@dataclass
class GatherResult(PrefetchResult):
#class GatherResult(SignatureComparison):
    current_gathersketch: MinHash = None
    gather_result_rank: int = None
    total_abund: int = None
    filename: str = None

    gather_write_cols = ['intersect_bp', 'f_orig_query', 'f_match', 'f_unique_to_query',
                         'f_unique_weighted','average_abund', 'median_abund', 'std_abund', 'filename',
                         'name', 'md5', 'match', 'f_match_orig', 'unique_intersect_bp', 'gather_result_rank',
                         'remaining_bp', 'query_filename', 'query_name', 'query_md5', 'query_bp', 'ksize',
                         'moltype', 'num', 'scaled', 'query_n_hashes', 'query_abundance']

    def init_gathersketchcomparison(self):
        # compare remaining gather hashes with match. Force at cmp_scaled.
        self.gather_comparison = FracMinHashComparison(self.current_gathersketch, self.match.minhash, cmp_scaled=self.cmp_scaled, threshold_bp=self.threshold_bp)

    def __post_init__(self):
        if self.cmp_scaled is None:
            raise ValueError("Error: must provide comparison scaled value ('cmp_scaled') for GatherResult")
        if self.current_gathersketch is None:
            raise ValueError("Error: must provide current gather sketch (remaining hashes) for GatherResult") # todo: fix this description
        if self.gather_result_rank is None:
            raise ValueError("Error: must provide 'gather_result_rank' to GatherResult")
        if self.total_abund is None:
            raise ValueError("Error: must provide 'sum_abunds' to GatherResult")
        self.init_sigcomparison() # initialize original sketch vs match sketch comparison
        self.init_gathersketchcomparison() # initialize remaining gather sketch vs match sketch comparison
        self.build_gather_result()
        if self.filename is None and self.match_filename is not None:
            self.filename = self.match_filename

    def gatherresultdict(self):
        # todo: standardize when we shorten the md5s between each type of result
        self.query_md5 = shorten_md5(self.query_md5)
        #self.md5 = shorten_md5(self.md5)
        #self.match_md5 = shorten_md5(self.match_md5)
        return self.to_write(columns=self.gather_write_cols)

    @property
    def writedict(self):
        return self.gatherresultdict()
    
    def build_gather_result(self):
        self.f_orig_query = self.query_containment
#       #containment of remaining, unaccounted for hashes with the database match
        self.f_unique_to_query =self.gather_comparison.mh1_containment
        # get current query-weighted intersection
        self.current_weighted_intersection =  self.gather_comparison.mh1_weighted_intersection
        self.f_unique_weighted =  self.gather_comparison.sum_mh_abundances(self.current_weighted_intersection)/ self.total_abund
        #original match containment
        self.f_match_orig = self.match_containment
        # unique match containment
        self.f_match = self.gather_comparison.mh2_containment
        #intersect bp of remaining, unaccounted for hashes with the database match
        self.unique_intersect_bp = self.gather_comparison.intersect_bp
        #current gather sketch bp - intersect bp OR PASS THIS IN?
        self.remaining_bp = self.query_bp - self.gather_comparison.intersect_bp
#        self.remaining_bp = self.gather_comparison.mh1_bp - self.gather_comparison.intersect_bp
        self.average_abund = self.gather_comparison.mean_mh_abundances(self.current_weighted_intersection)
        self.median_abund = self.gather_comparison.median_mh_abundances(self.current_weighted_intersection)
        self.std_abund = self.gather_comparison.std_mh_abundances(self.current_weighted_intersection)
        # how to we get properties to build? e.g. self.jaccard?
    
 

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
#    search_scaled = query.minhash.scaled # is this true?
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
    search_scaled = query.minhash.scaled ## is this true? or do i need max of match/query scaled?
    for (score, match, filename) in results:
        x.append(SearchResult(query, match,
                               search_scaled=search_scaled,
                               similarity=score,  # this is actually cosine sim? (abund)
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
        assert intersect_mh.contained_by(found_mh) == 1.0
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

        # compute weighted_missed for remaining query hashes
        query_hashes = set(query_mh.hashes) - set(found_mh.hashes)
        weighted_missed = sum((orig_query_abunds[k] for k in query_hashes))
        weighted_missed += self.noident_query_sum_abunds
        weighted_missed /= sum_abunds

        # build a GatherResult
        result = GatherResult(self.orig_query, best_match, # query or original query??
                              cmp_scaled=scaled,
                              filename=filename,
                              gather_result_rank=self.result_n,
                              total_abund= sum_abunds,
                              current_gathersketch=query_mh, # do we need to pass in found_mh too?
#                              intersect_bp=intersect_bp,
#                              unique_intersect_bp=unique_intersect_bp,
#                              f_orig_query=f_orig_query,
#                              f_match=f_match,
#                              f_match_orig=f_match_orig,
#                              f_unique_to_query=f_unique_to_query,
#                              f_unique_weighted=f_unique_weighted,
#                              average_abund=average_abund,
#                              median_abund=median_abund,
#                              std_abund=std_abund,
#                              remaining_bp=remaining_bp
                              )
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
