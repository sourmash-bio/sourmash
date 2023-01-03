"""
Sketch Comparison Classes
"""
import numpy as np
from dataclasses import dataclass

from .signature import MinHash

@dataclass
class BaseMinHashComparison:
    """Class for standard comparison between two MinHashes"""
    mh1: MinHash
    mh2: MinHash
    ignore_abundance: bool = False # optionally ignore abundances
    jaccard_ani_untrustworthy: bool = False

    def downsample_and_handle_ignore_abundance(self, cmp_num=None, cmp_scaled=None):
        """
        Downsample and/or flatten minhashes for comparison
        """
        if self.ignore_abundance:
            self.mh1_cmp = self.mh1.flatten()
            self.mh2_cmp = self.mh2.flatten()
        else:
            self.mh1_cmp = self.mh1
            self.mh2_cmp = self.mh2
        if cmp_scaled is not None:
            self.mh1_cmp = self.mh1_cmp.downsample(scaled=cmp_scaled)
            self.mh2_cmp = self.mh2_cmp.downsample(scaled=cmp_scaled)
        elif cmp_num is not None:
            self.mh1_cmp = self.mh1_cmp.downsample(num=cmp_num)
            self.mh2_cmp = self.mh2_cmp.downsample(num=cmp_num)
        else:
            raise ValueError("Error: must pass in a comparison scaled or num value.")

    def check_compatibility_and_downsample(self, cmp_num=None, cmp_scaled=None):
        if not any([(self.mh1.num and self.mh2.num), (self.mh1.scaled and self.mh2.scaled)]):
            raise TypeError("Error: Both sketches must be 'num' or 'scaled'.")

        #need to downsample first because is_compatible checks scaled (though does not check num)
        self.downsample_and_handle_ignore_abundance(cmp_num=cmp_num, cmp_scaled=cmp_scaled)
        if not self.mh1_cmp.is_compatible(self.mh2_cmp):
            raise TypeError("Error: Cannot compare incompatible sketches.")
        self.ksize = self.mh1.ksize
        self.moltype = self.mh1.moltype

    @property
    def intersect_mh(self):
        # flatten and intersect
        return self.mh1_cmp.flatten().intersection(self.mh2_cmp.flatten())

    @property
    def jaccard(self):
        return self.mh1_cmp.jaccard(self.mh2_cmp)

    def estimate_jaccard_ani(self, jaccard=None):
        jinfo = self.mh1_cmp.jaccard_ani(self.mh2_cmp, jaccard=jaccard)
        # propagate params
        self.jaccard_ani = jinfo.ani
        if jinfo.p_exceeds_threshold:
            self.potential_false_negative = True
        self.jaccard_ani_untrustworthy = jinfo.je_exceeds_threshold

    @property
    def angular_similarity(self):
        # Note: this currently throws TypeError if self.ignore_abundance.
        return self.mh1_cmp.angular_similarity(self.mh2_cmp)

    @property
    def cosine_similarity(self):
        return self.angular_similarity
    
@dataclass
class NumMinHashComparison(BaseMinHashComparison):
    """Class for standard comparison between two num minhashes"""
    cmp_num: int = None

    def __post_init__(self):
        "Initialize NumMinHashComparison using values from provided MinHashes"
        if self.cmp_num is None: # record the num we're doing this comparison on
            self.cmp_num = min(self.mh1.num, self.mh2.num)
        self.check_compatibility_and_downsample(cmp_num=self.cmp_num)

    @property
    def size_may_be_inaccurate(self):
        return False # not using size estimation, can ignore

@dataclass
class FracMinHashComparison(BaseMinHashComparison):
    """Class for standard comparison between two scaled minhashes"""
    cmp_scaled: int = None # optionally force scaled value for this comparison
    threshold_bp: int = 0
    estimate_ani_ci: bool = False
    ani_confidence: float = 0.95
#    pfn_threshold: float = 1e-3

    def __post_init__(self):
        "Initialize ScaledComparison using values from provided FracMinHashes"
        if self.cmp_scaled is None:
            # comparison scaled defaults to maximum scaled between the two sigs
            self.cmp_scaled = max(self.mh1.scaled, self.mh2.scaled)
        self.check_compatibility_and_downsample(cmp_scaled=self.cmp_scaled)
        self.potential_false_negative = False

    @property
    def pass_threshold(self):
        return self.total_unique_intersect_hashes >= self.threshold_bp

    @property
    def size_may_be_inaccurate(self):
        # if either size estimation may be inaccurate
        # NOTE: do we want to do this at original scaled instead?
        if not self.mh1_cmp.size_is_accurate() or not self.mh2_cmp.size_is_accurate():
            return True
        return False

    @property
    def total_unique_intersect_hashes(self):
        """
        approx equal to intersect_bp
        To get true bp estimates, we would need to add `(k-1)`. However, this complicates
        the iterative gather algorithm, so let's stick with hashes.
        """
        return len(self.intersect_mh) * self.cmp_scaled # + (ksize-1) #for bp estimation

    @property
    def mh1_containment_in_mh2(self):
        return self.mh1_cmp.contained_by(self.mh2_cmp)

    def estimate_ani_from_mh1_containment_in_mh2(self, containment = None):
        # build result once
        m1_cani = self.mh1_cmp.containment_ani(self.mh2_cmp,
                                            containment=containment,
                                            confidence=self.ani_confidence,
                                            estimate_ci=self.estimate_ani_ci)
#                                            prob_threshold=self.pfn_threshold)
        # propagate params
        self.ani_from_mh1_containment_in_mh2 = m1_cani.ani
        if m1_cani.p_exceeds_threshold:
            # only update if True
            self.potential_false_negative = True
        if self.estimate_ani_ci:
            self.ani_from_mh1_containment_in_mh2_low = m1_cani.ani_low
            self.ani_from_mh1_containment_in_mh2_high = m1_cani.ani_high

    @property
    def mh2_containment_in_mh1(self):
        return self.mh2_cmp.contained_by(self.mh1_cmp)

    def estimate_ani_from_mh2_containment_in_mh1(self, containment=None):
        m2_cani =  self.mh2_cmp.containment_ani(self.mh1_cmp,
                                            containment=containment,
                                            confidence=self.ani_confidence,
                                            estimate_ci=self.estimate_ani_ci)
#                                            prob_threshold=self.pfn_threshold)
        self.ani_from_mh2_containment_in_mh1 = m2_cani.ani
        if m2_cani.p_exceeds_threshold:
            self.potential_false_negative = True
        if self.estimate_ani_ci:
            self.ani_from_mh2_containment_in_mh1_low = m2_cani.ani_low
            self.ani_from_mh2_containment_in_mh1_high = m2_cani.ani_high
    
    @property
    def max_containment(self):
        return self.mh1_cmp.max_containment(self.mh2_cmp)

    def estimate_max_containment_ani(self, max_containment=None):
        mc_ani_info = self.mh1_cmp.max_containment_ani(self.mh2_cmp,
                                                max_containment=max_containment,
                                                confidence=self.ani_confidence,
                                                estimate_ci=self.estimate_ani_ci)
#                                                prob_threshold=self.pfn_threshold)
        # propagate params
        self.max_containment_ani = mc_ani_info.ani
        if mc_ani_info.p_exceeds_threshold:
            self.potential_false_negative = True
        if self.estimate_ani_ci:
            self.max_containment_ani_low = mc_ani_info.ani_low
            self.max_containment_ani_high = mc_ani_info.ani_high

    @property
    def avg_containment(self):
        return self.mh1_cmp.avg_containment(self.mh2_cmp)

    @property
    def avg_containment_ani(self):
        "Returns single average_containment_ani value. Sets self.potential_false_negative internally."
        self.estimate_ani_from_mh1_containment_in_mh2()
        self.estimate_ani_from_mh2_containment_in_mh1()
        if any([self.ani_from_mh1_containment_in_mh2 is None, self.ani_from_mh2_containment_in_mh1 is None]):
            return None
        else:
            return (self.ani_from_mh1_containment_in_mh2 + self.ani_from_mh2_containment_in_mh1)/2

    def estimate_all_containment_ani(self):
        "Estimate all containment ANI values."
        self.estimate_ani_from_mh1_containment_in_mh2()
        self.estimate_ani_from_mh2_containment_in_mh1()
        if any([self.ani_from_mh1_containment_in_mh2 is None, self.ani_from_mh2_containment_in_mh1 is None]):
#            self.estimate_max_containment_ani()
            self.max_containment_ani = None
        else:
            self.max_containment_ani = max([self.ani_from_mh1_containment_in_mh2, self.ani_from_mh2_containment_in_mh1])

    def weighted_intersection(self, from_mh=None, from_abundD={}):
         # map abundances to all intersection hashes.
        abund_mh = self.intersect_mh.copy_and_clear()
        abund_mh.track_abundance = True
        # if from_mh is provided, it takes precedence over from_abund dict
        if from_mh is not None and from_mh.track_abundance:
            from_abundD = from_mh.hashes
        if from_abundD:
            # this sets any hash not present in abundD to 1. Is that desired? Or should we return 0?
            abunds = {k: from_abundD.get(k, 1) for k in self.intersect_mh.hashes }
            abund_mh.set_abundances(abunds)
            return abund_mh
        # if no abundances are passed in, return intersect_mh
        # future note: do we want to return 1 as abundance instead?
        return self.intersect_mh
