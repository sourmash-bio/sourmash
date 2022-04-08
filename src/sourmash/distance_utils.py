"""
Utility functions for jaccard/containment --> distance estimates
From https://github.com/KoslickiLab/mutation-rate-ci-calculator
https://doi.org/10.1101/2022.01.11.475870
"""
from scipy.optimize import brentq
from scipy.stats import norm as scipy_norm
from scipy.special import hyp2f1
from numpy import sqrt
from math import log, exp

from .logging import notify


def r1_to_q(k,r1):
    r1 = float(r1)
    q = 1-(1-r1)**k
    return float(q)


def var_n_mutated(L,k,r1,q=None):
	# there are computational issues in the variance formula that we solve here
	# by the use of higher-precision arithmetic; the problem occurs when r is
	# very small; for example, with L=10,k=2,r1=1e-6 standard precision
	# gives varN<0 which is nonsense; by using the mpf type, we get the correct
	# answer which is about 0.000038.
	if (r1 == 0): return 0.0
	r1 = float(r1)
	if (q == None): # we assume that if q is provided, it is correct for r1
	    q = r1_to_q(k,r1)
	varN = L*(1-q)*(q*(2*k+(2/r1)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
         + (2*(1-q)/(r1**2))*((1+(k-1)*(1-q))*r1-q)
	if (varN<0.0): # this seems to happen only with super tiny test data
         raise ValueError('Error: varN <0.0!')
	return float(varN)


def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q


def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2


def probit(p):
    return scipy_norm.ppf(p)


def sequence_len_to_n_kmers(sequence_len_bp, ksize):
    n_kmers = sequence_len_bp - (ksize-1)
    return n_kmers


def distance_to_identity(dist,d_low=None,d_high=None):
    """
    ANI = 1-distance
    """
    if not 0 <= dist <= 1:
        raise ValueError(f"Error: distance value {dist} is not between 0 and 1!")
    ident = 1-dist
    result = ident
    id_low,id_high=None,None
    if d_low is not None and d_high is not None:
        if (0<=d_low<=1) and (0<=d_high<=1):
            id_high = 1-d_low
            id_low = 1-d_high
            result = [ident, id_low, id_high]
    return result


def get_expected_log_probability(n_unique_kmers, ksize, mutation_rate, scaled_fraction):
    '''helper function
       Note that scaled here needs to be between 0 and 1
       (e.g. scaled 1000 --> scaled_fraction 0.001)
    '''
    exp_nmut = exp_n_mutated(n_unique_kmers, ksize, mutation_rate)
    try:
        return (n_unique_kmers - exp_nmut) * log(1.0 - scaled_fraction)
    except:
        return float('-inf')


def get_exp_probability_nothing_common(mutation_rate, ksize, scaled, n_unique_kmers=None, sequence_len_bp=None):
    """
    Given parameters, calculate the expected probability that nothing will be common
    between a fracminhash sketch of a original sequence and a fracminhash sketch of a mutated
    sequence. If this is above a threshold, we should suspect that the two sketches may have
    nothing in common. The threshold needs to be set with proper insights.

    Arguments: n_unique_kmers, ksize, mutation_rate, scaled
    Returns: float - expected likelihood that nothing is common between sketches
    """
    # NTP note: do we have any checks for ksize >=1 in sourmash_args? The rest should be taken care of.
    # assert 0.0 <= mutation_rate <= 1.0 and ksize >= 1 and scaled >= 1
    if sequence_len_bp and not n_unique_kmers:
        n_unique_kmers = sequence_len_to_n_kmers(sequence_len_bp, ksize)
    f_scaled = 1.0 / float(scaled)
    if mutation_rate == 1.0:
        return 1.0
    elif mutation_rate == 0.0:
        return 0.0
    return exp( get_expected_log_probability(n_unique_kmers, ksize, mutation_rate, f_scaled) )


def containment_to_distance(containment, ksize, scaled, n_unique_kmers=None, sequence_len_bp=None, confidence=0.95, return_identity=False, return_ci = False, prob_threshold = 10.0**(-3)):
    """
    Containment --> distance CI (one step)
    """
    sol1,sol2,point_estimate=None,None,None
    if sequence_len_bp and not n_unique_kmers:
        n_unique_kmers = sequence_len_to_n_kmers(sequence_len_bp, ksize)
    if containment <= 0.0001:
        sol2 = sol1 = point_estimate = 1.0
    elif containment >= 0.9999:
        sol1 = sol2 = point_estimate = 0.0
    else:
        point_estimate = 1.0-containment**(1.0/ksize)
        if return_ci:
            try:
                alpha = 1 - confidence
                z_alpha = probit(1-alpha/2)
                f_scaled = 1.0/scaled # these use scaled as a fraction between 0 and 1

                bias_factor = 1 - (1 - f_scaled) ** n_unique_kmers

                term_1 = (1.0-f_scaled) / (f_scaled * n_unique_kmers**3 * bias_factor**2)
                term_2 = lambda pest: n_unique_kmers * exp_n_mutated(n_unique_kmers, ksize, pest) - exp_n_mutated_squared(n_unique_kmers, ksize, pest)
                term_3 = lambda pest: var_n_mutated(n_unique_kmers, ksize, pest) / (n_unique_kmers**2)

                var_direct = lambda pest: term_1 * term_2(pest) + term_3(pest)

                f1 = lambda pest: (1-pest)**ksize + z_alpha * sqrt(var_direct(pest)) - containment
                f2 = lambda pest: (1-pest)**ksize - z_alpha * sqrt(var_direct(pest)) - containment

                sol1 = brentq(f1, 0.0000001, 0.9999999)
                sol2 = brentq(f2, 0.0000001, 0.9999999)

            except ValueError as exc:
                # afaict, this only happens with extremely small test data
                notify("WARNING: Cannot estimate ANI from containment. Do your sketches contain enough hashes?")
                notify(str(exc))
                return None, None, None

    # Do this here, so that we don't need to reconvert distance <--> identity later.
    prob_nothing_in_common = get_exp_probability_nothing_common(point_estimate, ksize, scaled, n_unique_kmers=n_unique_kmers)
    if prob_nothing_in_common >= prob_threshold:
        # TO DO: keep count; suggest user decrease scaled value. If that is unsuccessful, maybe decrease ksize
        notify('WARNING: These sketches may have no hashes in common based on chance alone.')
       
    if return_identity:
        if any([sol1 is None, sol2 is None]):
            point_estimate = distance_to_identity(point_estimate)
        else:
            point_estimate,sol2,sol1 = distance_to_identity(point_estimate,sol2,sol1)

    if return_ci:
        return point_estimate, sol2, sol1, prob_nothing_in_common

    return point_estimate, prob_nothing_in_common


def jaccard_to_distance(jaccard, ksize, scaled, n_unique_kmers=None, sequence_len_bp=None, return_identity=False, prob_threshold = 10.0**(-3), err_threshold=10.0**(-4.0)):
    """
    Given parameters, calculate point estimate for mutation rate from jaccard index.
    First checks if parameters are valid (checks are not exhaustive). Then uses formulas
    derived mathematically to compute the point estimate. The formula uses approximations,
    therefore a tiny error is associated with it. A lower bound of that error is also returned.
    A high error indicates that the point estimate cannot be trusted. Threshold of the error
    is open to interpretation, but suggested that > 10^-4 should be handled with caution.

    Note that the error is NOT a mutation rate, and therefore cannot be considered in
    something like mut.rate +/- error.

    Arguments: jaccard, ksize, scaled, n_unique_kmers
    # Returns: tuple (point_estimate_of_mutation_rate, lower_bound_of_error)

    # Returns: point_estimate_of_mutation_rate

    Note: point estimate does not consider impact of scaled, but p_nothing_in_common can be
    useful for determining whether scaled is sufficient for these comparisons.
    """
    error_lower_bound = None
    if sequence_len_bp and not n_unique_kmers:
        n_unique_kmers = sequence_len_to_n_kmers(sequence_len_bp, ksize)
    if jaccard <= 0.0001:
        point_estimate = 1.0
        error_lower_bound = 0.0
    elif jaccard >= 0.9999:
        point_estimate = 0.0
        error_lower_bound = 0.0
    else:
        point_estimate = 1.0 - ( 2.0 * jaccard / float(1+jaccard) ) ** (1.0/float(ksize))

        exp_n_mut = exp_n_mutated(n_unique_kmers, ksize, point_estimate)
        var_n_mut = var_n_mutated(n_unique_kmers, ksize, point_estimate)
        error_lower_bound = 1.0 * n_unique_kmers * var_n_mut / (n_unique_kmers + exp_n_mut)**3
    
    if error_lower_bound is not None and error_lower_bound > err_threshold:
        notify(f"WARNING: Error on Jaccard distance point estimate is too high ({error_lower_bound}).")
    
    prob_nothing_in_common = get_exp_probability_nothing_common(point_estimate, ksize, scaled, n_unique_kmers=n_unique_kmers)
    if prob_nothing_in_common >= prob_threshold:
        # to do: keep count and recommend user lower scaled val
        notify('WARNING: These sketches may have no hashes in common based on chance alone.')

    if return_identity:
        point_estimate = distance_to_identity(point_estimate)

    return point_estimate, prob_nothing_in_common, error_lower_bound