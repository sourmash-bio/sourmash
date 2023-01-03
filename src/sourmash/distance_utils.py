"""
Utilities for jaccard/containment --> distance estimation
Equations from: https://github.com/KoslickiLab/mutation-rate-ci-calculator
Reference: https://doi.org/10.1101/2022.01.11.475870
"""
from dataclasses import dataclass, field
from scipy.optimize import brentq
from scipy.stats import norm as scipy_norm
from scipy.stats import binom
import numpy as np
from math import log, exp

from .logging import notify

def check_distance(dist):
    if not 0 <= dist <= 1:
        raise ValueError(f"Error: distance value {dist :.4f} is not between 0 and 1!")
    else:
        return dist

def check_prob_threshold(val, threshold=1e-3):
    """
    Check likelihood of no shared hashes based on chance alone (false neg).
    If too many exceed threshold, recommend user lower their scaled value.
    # !! when using this, keep count and recommend user lower scaled val
    """
    exceeds_threshold = False
    if threshold is not None and val > threshold:
        exceeds_threshold = True
    return val, exceeds_threshold

def check_jaccard_error(val, threshold=1e-4):
    exceeds_threshold = False
    if threshold is not None and val > threshold:
        exceeds_threshold = True
    return val, exceeds_threshold

@dataclass
class ANIResult:
    """Base class for distance/ANI from k-mer containment."""
    dist: float
    p_nothing_in_common: float
    p_threshold: float = 1e-3
    size_is_inaccurate: bool = False
    p_exceeds_threshold: bool = field(init=False)

    def check_dist_and_p_threshold(self):
        # check values
        self.dist = check_distance(self.dist)
        self.p_nothing_in_common, self.p_exceeds_threshold = check_prob_threshold(self.p_nothing_in_common, self.p_threshold)

    def __post_init__(self):
        self.check_dist_and_p_threshold()

    @property
    def ani(self):
        if self.size_is_inaccurate:
            return None
        return 1 - self.dist


@dataclass
class jaccardANIResult(ANIResult):
    """Class for distance/ANI from jaccard (includes jaccard_error)."""
    jaccard_error: float = None
    je_threshold: float = 1e-4

    def __post_init__(self):
        # check values
        self.check_dist_and_p_threshold()
        # check jaccard error
        if self.jaccard_error is not None:
            self.jaccard_error, self.je_exceeds_threshold = check_jaccard_error(self.jaccard_error, self.je_threshold)
        else:
            raise ValueError("Error: jaccard_error cannot be None.")

    @property
    def ani(self):
        # if jaccard error is too high (exceeds threshold), do not trust ANI estimate
        if self.je_exceeds_threshold or self.size_is_inaccurate:
            return None
        return 1 - self.dist


@dataclass
class ciANIResult(ANIResult):
    """
    Class for distance/ANI from containment: with confidence intervals.

    Set CI defaults to None, just in case CI can't be estimated for given sample.
    """
    dist_low: float = None
    dist_high: float = None

    def __post_init__(self):
        # check values
        self.check_dist_and_p_threshold()

        if self.dist_low is not None and self.dist_high is not None:
            self.dist_low = check_distance(self.dist_low)
            self.dist_high = check_distance(self.dist_high)

    @property
    def ani_low(self):
        if self.dist_high is None or self.size_is_inaccurate:
            return None
        return 1 - self.dist_high

    @property
    def ani_high(self):
        if self.dist_low is None or self.size_is_inaccurate:
            return None
        return 1 - self.dist_low


def r1_to_q(k, r1):
    r1 = float(r1)
    q = 1 - (1 - r1) ** k
    return float(q)


def var_n_mutated(L, k, r1, *, q=None):
    # there are computational issues in the variance formula that we solve here
    # by the use of higher-precision arithmetic; the problem occurs when r is
    # very small; for example, with L=10,k=2,r1=1e-6 standard precision
    # gives varN<0 which is nonsense; by using the mpf type, we get the correct
    # answer which is about 0.000038.
    if r1 == 0:
        return 0.0
    r1 = float(r1)
    if q == None:  # we assume that if q is provided, it is correct for r1
        q = r1_to_q(k, r1)
    varN = (
        L * (1 - q) * (q * (2 * k + (2 / r1) - 1) - 2 * k)
        + k * (k - 1) * (1 - q) ** 2
        + (2 * (1 - q) / (r1**2)) * ((1 + (k - 1) * (1 - q)) * r1 - q)
    )
    if varN < 0.0:  # this seems to happen only with super tiny test data
        raise ValueError("Error: varN <0.0!")
    return float(varN)


def exp_n_mutated(L, k, r1):
    q = r1_to_q(k, r1)
    return L * q


def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2


def probit(p):
    return scipy_norm.ppf(p)


def handle_seqlen_nkmers(ksize, *, sequence_len_bp=None, n_unique_kmers=None):
    if n_unique_kmers is not None:
        return n_unique_kmers
    elif sequence_len_bp is None:
        # both are None, raise ValueError
        raise ValueError("Error: distance estimation requires input of either 'sequence_len_bp' or 'n_unique_kmers'")
    else:
        n_unique_kmers = sequence_len_bp - (ksize - 1)
        return n_unique_kmers


def set_size_chernoff(set_size, scaled, *, relative_error=0.05):
    """
    Computes the probability that the estimate: sketch_size * scaled deviates from the true
    set_size by more than relative_error. This relies on the fact that the sketch_size
    is binomially distributed with parameters sketch_size and 1/scale. The two-sided Chernoff
    bounds are used. This is depreciated in favor of set_size_exact_prob due to the later
    being accurate even for very small set sizes
    @param set_size: The number of distinct k-mers in the given set
    @param relative_error: the desired relative error (defaults to 5%)
    @return: float (the upper bound probability)
    """
    upper_bound = 1 - 2 * np.exp(- relative_error**2*set_size/(scaled * 3))
    return upper_bound


def set_size_exact_prob(set_size, scaled, *, relative_error=0.05):
    """
    Computes the exact probability that the estimate: sketch_size * scaled deviates from the true
    set_size by more than relative_error. This relies on the fact that the sketch_size
    is binomially distributed with parameters sketch_size and 1/scale. The CDF of the binomial distribution
    is used.
    @param set_size: The number of distinct k-mers in the given set
    @param relative_error: the desired relative error (defaults to 5%)
    @return: float (the upper bound probability)
    """
    # Need to check if the edge case is an integer or not. If not, don't include it in the equation
    pmf_arg = -set_size/scaled * (relative_error - 1)
    if pmf_arg == int(pmf_arg):
        prob = binom.cdf(set_size/scaled * (relative_error + 1), set_size, 1/scaled) - \
               binom.cdf(-set_size/scaled * (relative_error - 1), set_size, 1/scaled) + \
               binom.pmf(-set_size/scaled * (relative_error - 1), set_size, 1/scaled)
    else:
        prob = binom.cdf(set_size / scaled * (relative_error + 1), set_size, 1 / scaled) - \
               binom.cdf(-set_size / scaled * (relative_error - 1), set_size, 1 / scaled)
    return prob


def get_expected_log_probability(n_unique_kmers, ksize, mutation_rate, scaled_fraction):
    """helper function
    Note that scaled here needs to be between 0 and 1
    (e.g. scaled 1000 --> scaled_fraction 0.001)
    """
    exp_nmut = exp_n_mutated(n_unique_kmers, ksize, mutation_rate)
    try:
        return (n_unique_kmers - exp_nmut) * log(1.0 - scaled_fraction)
    except:
        return float("-inf")


def get_exp_probability_nothing_common(
    mutation_rate, ksize, scaled, *, n_unique_kmers=None, sequence_len_bp=None
):
    """
    Given parameters, calculate the expected probability that nothing will be common
    between a fracminhash sketch of a original sequence and a fracminhash sketch of a mutated
    sequence. If this is above a threshold, we should suspect that the two sketches may have
    nothing in common. The threshold needs to be set with proper insights.

    Arguments: n_unique_kmers, ksize, mutation_rate, scaled
    Returns: float - expected likelihood that nothing is common between sketches
    """
    n_unique_kmers = handle_seqlen_nkmers(ksize, sequence_len_bp=sequence_len_bp,n_unique_kmers=n_unique_kmers)
    f_scaled = 1.0 / float(scaled)
    if mutation_rate == 1.0:
        return 1.0
    elif mutation_rate == 0.0:
        return 0.0
    return exp(
        get_expected_log_probability(n_unique_kmers, ksize, mutation_rate, f_scaled)
    )


def containment_to_distance(
    containment,
    ksize,
    scaled,
    *,
    n_unique_kmers=None,
    sequence_len_bp=None,
    confidence=0.95,
    estimate_ci=False,
    prob_threshold=1e-3,
):
    """
    Containment --> distance CI (one step)
    """
    sol1, sol2, point_estimate = None, None, None
    n_unique_kmers = handle_seqlen_nkmers(ksize, sequence_len_bp = sequence_len_bp, n_unique_kmers=n_unique_kmers)
    if containment == 0:
        #point_estimate = 1.0
        point_estimate = sol1 = sol2 = 1.0
    elif containment == 1:
        #point_estimate = 0.0
        point_estimate = sol1 = sol2 = 0.0
    else:
        point_estimate = 1.0 - containment ** (1.0 / ksize)
        if estimate_ci:
            try:
                alpha = 1 - confidence
                z_alpha = probit(1 - alpha / 2)
                f_scaled = (
                    1.0 / scaled
                )  # these use scaled as a fraction between 0 and 1

                bias_factor = 1 - (1 - f_scaled) ** n_unique_kmers

                term_1 = (1.0 - f_scaled) / (
                    f_scaled * n_unique_kmers**3 * bias_factor**2
                )
                term_2 = lambda pest: n_unique_kmers * exp_n_mutated(
                    n_unique_kmers, ksize, pest
                ) - exp_n_mutated_squared(n_unique_kmers, ksize, pest)
                term_3 = lambda pest: var_n_mutated(n_unique_kmers, ksize, pest) / (
                    n_unique_kmers**2
                )

                var_direct = lambda pest: term_1 * term_2(pest) + term_3(pest)

                f1 = (
                    lambda pest: (1 - pest) ** ksize
                    + z_alpha * np.sqrt(var_direct(pest))
                    - containment
                )
                f2 = (
                    lambda pest: (1 - pest) ** ksize
                    - z_alpha * np.sqrt(var_direct(pest))
                    - containment
                )

                sol1 = brentq(f1, 0.0000001, 0.9999999)
                sol2 = brentq(f2, 0.0000001, 0.9999999)

            except ValueError as exc:
                # afaict, this only happens with extremely small test data
                notify(
                    "WARNING: Cannot estimate ANI confidence intervals from containment. Do your sketches contain enough hashes?"
                )
                notify(str(exc))
                sol1 = sol2 = None

    # Do this here, so that we don't need to reconvert distance <--> identity later.
    prob_nothing_in_common = get_exp_probability_nothing_common(
        point_estimate, ksize, scaled, n_unique_kmers=n_unique_kmers
    )
    return ciANIResult(point_estimate, prob_nothing_in_common, dist_low=sol2, dist_high=sol1, p_threshold=prob_threshold)


def jaccard_to_distance(
    jaccard,
    ksize,
    scaled,
    *,
    n_unique_kmers=None,
    sequence_len_bp=None,
    prob_threshold=1e-3,
    err_threshold=1e-4,
):
    """
    Given parameters, calculate point estimate for mutation rate from jaccard index.
    Uses formulas derived mathematically to compute the point estimate. The formula uses
    approximations, therefore a tiny error is associated with it. A lower bound of that error
    is also returned. A high error indicates that the point estimate cannot be trusted.
    Threshold of the error is open to interpretation, but suggested that > 10^-4 should be
    handled with caution.

    Note that the error is NOT a mutation rate, and therefore cannot be considered in
    something like mut.rate +/- error.

    Arguments: jaccard, ksize, scaled, n_unique_kmers
    # Returns: tuple (point_estimate_of_mutation_rate, lower_bound_of_error)

    # Returns: JaccardANIResult

    Note: point estimate does not consider impact of scaled, but p_nothing_in_common can be
    useful for determining whether scaled is sufficient for these comparisons.
    """
    error_lower_bound = None
    n_unique_kmers = handle_seqlen_nkmers(ksize, sequence_len_bp=sequence_len_bp, n_unique_kmers=n_unique_kmers)
    if jaccard == 0:
        point_estimate = 1.0
        error_lower_bound = 0.0
    elif jaccard == 1:
        point_estimate = 0.0
        error_lower_bound = 0.0
    else:
        point_estimate = 1.0 - (2.0 * jaccard / float(1 + jaccard)) ** (
            1.0 / float(ksize)
        )

        exp_n_mut = exp_n_mutated(n_unique_kmers, ksize, point_estimate)
        var_n_mut = var_n_mutated(n_unique_kmers, ksize, point_estimate)
        error_lower_bound = (
            1.0 * n_unique_kmers * var_n_mut / (n_unique_kmers + exp_n_mut) ** 3
        )
    prob_nothing_in_common = get_exp_probability_nothing_common(
        point_estimate, ksize, scaled, n_unique_kmers=n_unique_kmers
    )
    return jaccardANIResult(point_estimate, prob_nothing_in_common, jaccard_error=error_lower_bound, p_threshold=prob_threshold, je_threshold=err_threshold)
