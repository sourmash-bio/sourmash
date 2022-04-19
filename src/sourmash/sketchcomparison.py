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

    def check_comparison_compatibility(self):
        # do we need this check + error? Minhash functions should already complain appropriately...
        k1 = self.mh1.ksize
        k2 = self.mh2.ksize
        if k1 != k2:
            raise ValueError(f"Error: Invalid Comparison, ksizes: {k1}, {k2}). Must compare sketches of the same ksize.")
        self.ksize = self.mh1.ksize
        m1 = self.mh2.moltype
        m2= self.mh2.moltype
        if m1 != m2:
            raise ValueError(f"Error: Invalid Comparison, moltypes: {m1}, {m2}). Must compare sketches of the same moltype.")
        self.moltype= self.mh1.moltype
        # check num, scaled
        if not any([(self.mh1.num and self.mh2.num), (self.mh1.scaled and self.mh2.scaled)]):
            raise ValueError("Error: Both sketches must be 'num' or 'scaled'.")

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

    @property
    def intersect_mh(self):
        # flatten and intersect
        return self.mh1_cmp.flatten().intersection(self.mh2_cmp.flatten())

    @property
    def jaccard(self):
        return self.mh1_cmp.jaccard(self.mh2_cmp)

    @property
    def angular_similarity(self):
        return self.mh1_cmp.angular_similarity(self.mh2_cmp)
        # do we want to shield against error here? Or let TypeError through?
        #if not (self.mh1_cmp.track_abundance and self.mh2_cmp.track_abundance):
        #    return self.mh1_cmp.angular_similarity(self.mh2_cmp)
        #else:
        #    return ""

    @property
    def cosine_similarity(self):
        return self.angular_similarity


@dataclass
class NumMinHashComparison(BaseMinHashComparison):
    """Class for standard comparison between two scaled minhashes"""
    cmp_num: int = None

    def __post_init__(self):
        "Initialize NumMinHashComparison using values from provided MinHashes"
        if self.cmp_num is None: # record the num we're doing this compa#rison on
            self.cmp_num = min(self.mh1.num, self.mh2.num)
        self.check_comparison_compatibility()
        self.downsample_and_handle_ignore_abundance(cmp_num=self.cmp_num)

@dataclass
class FracMinHashComparison(BaseMinHashComparison):
    """Class for standard comparison between two scaled minhashes"""
    cmp_scaled: int = None # scaled value for this comparison (defaults to maximum scaled between the two sigs)
    threshold_bp: int = 0

    def __post_init__(self):
        "Initialize ScaledComparison using values from provided FracMinHashes"
        if self.cmp_scaled is None: # record the scaled we're doing this comparison on
            self.cmp_scaled = max(self.mh1.scaled, self.mh2.scaled)
        self.check_comparison_compatibility()
        self.downsample_and_handle_ignore_abundance(cmp_scaled=self.cmp_scaled)
        # for these, do we want the originals, or the cmp_scaled versions?? (or both?). do we need them at all?
        self.mh1_scaled = self.mh1.scaled
        self.mh2_scaled = self.mh2.scaled
        self.mh1_n_hashes = len(self.mh1)
        self.mh2_n_hashes = len(self.mh2)
        self.mh1_bp = self.mh1.covered_bp
        self.mh2_bp = self.mh2.covered_bp

    @property
    def pass_threshold(self):
        return self.intersect_bp >= self.threshold_bp

    @property
    def intersect_bp(self):
        return len(self.intersect_mh) * self.cmp_scaled

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
        return np.mean([self.mh1_containment, self.mh2_containment])

    @property
    def mh1_weighted_intersection(self):
         # map mh1 hash abundances to all intersection hashes.
        if self.ignore_abundance or not self.mh1.track_abundance:
            return self.intersect_mh # or do we want to set all abundances to 1?
        else:
        # current functionality inflates by original minhash (not downsampled)-- is that what we want??
            return self.intersect_mh.inflate(self.mh1)

    @property
    def mh2_weighted_intersection(self):
         # map mh2 hash abundances to all intersection hashes.
        if self.ignore_abundance or not self.mh2.track_abundance:
            return self.intersect_mh  # or do we want to set all abundances to 1?
        else:
        # current functionality inflates by original minhash (not downsampled)-- is that what we want??
            return self.intersect_mh.inflate(self.mh2)