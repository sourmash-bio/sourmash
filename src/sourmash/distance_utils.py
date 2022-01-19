"""
Utility functions for jaccard/containment --> distance estimates
From https://github.com/KoslickiLab/mutation-rate-ci-calculator
https://doi.org/10.1101/2022.01.11.475870
"""
from scipy.optimize import brentq
from scipy.stats import norm as scipy_norm
from scipy.special import hyp2f1
from numpy import sqrt

from .logging import notify

#FROM  mrcc.kmer_mutation_formulas_thm5
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
	if (varN<0.0):
         raise ValueError('Error: varN <0.0!')
	return float(varN)


def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q


# FROM mrcc.third_moment_calculator
def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2


def third_moment_nmut_exact(L,k,p):
    t1 = (L * (-2 + 3*L) * p**2 + 3 * (1 - p)**(2*k) * (2 + (-1 + k - L) * p * (2 + k * p - L * p)) - (1 - p)**k * (6 + p * (-6 + L * (-6 + p + 6 * L * p))))/(p**2)
    t2 = (-2 + 2 * k - L) * (-1 + 2 * k - L) * (2 * k - L) * (-1 + (1 - p)**k)**3
    t3 = (1/(p**3))*(-6 * (-1 + k)**2 * (k - L) * p**3 + 6 * (1 - p)**(3 * k) * (2 + (-2 + 2 * k - L) * p) + (1 - p)**(2 * k) * (-12 + 6 * (2 * k + L) * p + 6 * (4 * k**2 + 2 * (1 + L) - 3 * k * (2 + L)) * p**2 - (-1 + k) * k * (-2 + 4 * k - 3 * L) * p**3) + 6 * (-1 + k) * (1 - p)**k * p * (-2 + p * (2 - k + 2 * L + (k * (-2 + 3 * k - 3 * L) + L) * p)))
    t4 = 6 * (-1 + (1 - p)**k) * ((k + k**2 - 2 * k * L + (-1 + L) * L) * (-1 + 2 * (1 - p)**k) * hyp2f1(1, 2 + k - L, k - L, 1) + (k + k**2 - 2 * k * L + (-1 + L) * L) * (1 - p)**k * (-1 + p) * hyp2f1(1, 2 + k - L, k - L, 1 - p) - (-2 * k + 4 * k**2 + L - 4 * k * L + L**2) * ((-1 + 2 * (1 - p)**k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1)- (-1 + p)**(2 * k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1 - p)))
    return t1+t2+t3+t4


def exp_n_mutated_cubed(L, k, p):
    return third_moment_nmut_exact(L, k, p)


#from mrcc.p_from_scaled_containment
def probit(p):
    return scipy_norm.ppf(p)


def sequence_len_to_n_kmers(sequence_len_bp, ksize):
    n_kmers = sequence_len_bp - (ksize-1)
    return n_kmers


def containment_to_distance(containment, ksize, scaled, n_unique_kmers=None, sequence_len_bp=None, confidence=0.95, return_identity=False):
    """
    Containment --> distance CI (one step)
    """
    if sequence_len_bp and not n_unique_kmers:
        n_unique_kmers = sequence_len_to_n_kmers(sequence_len_bp, ksize)
    if containment <= 0.0001: # changed from 0.0, to mirror jaccard_to_distance_CI_one_step
        sol2 = sol1 = point_estimate = 1.0
    elif containment >= 0.9999:  # changed from 1.0, to mirror jaccard_to_distance_CI_one_step
        sol1 = sol2 = point_estimate = 0.0
    else:
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

        point_estimate = 1.0-containment**(1.0/ksize)
    if return_identity:
        point_estimate,sol2,sol1 = distance_to_identity(point_estimate,sol2,sol1)
        #distance_using_CImidpoint = (sol2+sol1)/2.0
    return point_estimate,sol2,sol1


#from from mrcc.p_from_scaled_jaccard
def variance_scaled_jaccard(L, p, k, s):
    exp_n_mut = exp_n_mutated(L, k, p)
    exp_n_mut_squared = exp_n_mutated_squared(L, k, p)
    exp_n_mut_cubed = exp_n_mutated_cubed(L, k, p)
    bias_factor = 1 - (1 - s) ** ( int(L + exp_n_mut) )

    factor1 = (1-s)/(s * bias_factor**2)
    factor2 = (2 * L * exp_n_mut - 2 * exp_n_mut_squared) / (L ** 3 + 3*L*exp_n_mut_squared + 3*L*L*exp_n_mut + exp_n_mut_cubed)
    term1 = factor1 * factor2
    term2 = (L**2 - 2 * L * exp_n_mut + exp_n_mut_squared) / (L**2 + 2 * L * exp_n_mut + exp_n_mut_squared)
    term3 = ((L - exp_n_mut) / (L + exp_n_mut))**2

    return term1 + term2 - term3


def jaccard_to_distance(jaccard, ksize, scaled, n_unique_kmers=None, sequence_len_bp=None, confidence=0.95, return_identity=False):
    """
    Scaled Jaccard to distance estimate and CI (one step)
    """
    if sequence_len_bp and not n_unique_kmers:
        n_unique_kmers = sequence_len_to_n_kmers(sequence_len_bp, ksize)
    if jaccard <= 0.0001:
        sol2 = sol1 = point_estimate = 1.0
    elif jaccard >= 0.9999:
        sol1 = sol2 = point_estimate = 0.0
    else:
        alpha = 1 - confidence
        z_alpha = probit(1-alpha/2)
        f_scaled = 1.0/scaled # these use scaled as a fraction between 0 and 1

        var_direct = lambda pest: variance_scaled_jaccard(n_unique_kmers, pest, ksize, f_scaled)

        f1 = lambda pest: 2.0/(2- (1-pest)**ksize ) - 1 + z_alpha * sqrt(var_direct(pest)) - jaccard
        f2 = lambda pest: 2.0/(2- (1-pest)**ksize ) - 1 - z_alpha * sqrt(var_direct(pest)) - jaccard

        sol1 = brentq(f1, 0.0001, 0.9999)
        sol2 = brentq(f2, 0.0001, 0.9999)

        point_estimate = 1.0 - (2.0 - 2.0/(jaccard + 1))**(1.0/ksize)
    if return_identity:
        point_estimate,sol2,sol1 = distance_to_identity(point_estimate,sol2,sol1)

    return point_estimate,sol2,sol1


def distance_to_identity(dist,d_low=None,d_high=None):
    """
    ANI = 1-distance
    """
    for d in [dist,d_low,d_high]:
        if not 0 <= d <= 1:
            raise ValueError(f"Error: distance value {d} is not between 0 and 1!")
    id = 1-dist
    id_low,id_high=None,None
    if d_low is not None: # need to be explicit so will work on 0 value
        id_high = 1-d_low
    if d_high is not None: # need to be explicit so will work on 0 value
        id_low = 1-d_high
    return id,id_low,id_high