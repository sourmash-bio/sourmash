"""
Tests for distance utils.
"""
import pytest
import numpy as np
from sourmash.distance_utils import (containment_to_distance, get_exp_probability_nothing_common,
                                    handle_seqlen_nkmers, jaccard_to_distance,
                                    ANIResult, ciANIResult, jaccardANIResult, var_n_mutated,
                                    set_size_chernoff, set_size_exact_prob)

def test_aniresult():
    res = ANIResult(0.4, 0.1)
    assert res.dist == 0.4
    assert res.ani == 0.6
    assert res.p_nothing_in_common == 0.1
    assert res.p_exceeds_threshold ==True
    # check that they're equivalent
    res2 = ANIResult(0.4, 0.1)
    assert res == res2
    res3 = ANIResult(0.5, 0)
    assert res != res3
    assert res3.p_exceeds_threshold ==False

def test_aniresult_bad_distance():
    """
    Fail if distance is not between 0 and 1.
    """
    with pytest.raises(Exception) as exc:
        ANIResult(1.1, 0.1)
    print("\n", str(exc.value))
    assert "distance value 1.1000 is not between 0 and 1!" in str(exc.value)
    with pytest.raises(Exception) as exc:
        ANIResult(-0.1, 0.1)
    print("\n", str(exc.value))
    assert "distance value -0.1000 is not between 0 and 1!" in str(exc.value)


def test_jaccard_aniresult():
    res = jaccardANIResult(0.4, 0.1, jaccard_error=0.03)
    assert res.dist == 0.4
    assert res.ani == None
    assert res.p_nothing_in_common == 0.1
    assert res.jaccard_error == 0.03
    assert res.p_exceeds_threshold ==True
    assert res.je_exceeds_threshold ==True
    res3 = jaccardANIResult(0.4, 0.1, jaccard_error=0.03, je_threshold=0.1)
    assert res3.je_exceeds_threshold ==False
    assert res3.ani == 0.6


def test_jaccard_aniresult_nojaccarderror():
    #jaccard error is None
    with pytest.raises(Exception) as exc:
        jaccardANIResult(0.4, 0.1, None)
    print("\n", str(exc.value))
    assert "Error: jaccard_error cannot be None." in str(exc.value)


def test_ci_aniresult():
    res = ciANIResult(0.4, 0.1, dist_low=0.3,dist_high=0.5)
    print(res)
    assert res.dist == 0.4
    assert res.ani == 0.6
    assert res.p_nothing_in_common == 0.1
    assert res.ani_low == 0.5
    assert res.ani_high == 0.7
    res2 = ciANIResult(0.4, 0.1, dist_low=0.3,dist_high=0.5)
    assert res == res2
    res3 = ciANIResult(0.4, 0.2, dist_low=0.3, dist_high=0.5)
    assert res != res3


def test_containment_to_distance_zero():
    contain = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    res = containment_to_distance(contain,ksize,scaled, n_unique_kmers=nkmers, estimate_ci=True)
    print(res)
    # check results
    exp_dist,exp_low,exp_high,pnc = 1.0,1.0,1.0,1.0
    exp_id, exp_idlow,exp_idhigh,pnc = 0.0,0.0,0.0,1.0
    assert res.dist == exp_dist
    assert res.dist_low == exp_low
    assert res.dist_high == exp_high
    assert res.p_nothing_in_common == pnc
    assert res.ani == exp_id
    assert res.ani_low == exp_idlow
    assert res.ani_high == exp_idhigh
    # check without returning ci
    res2 = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print(res2)
    exp_res = ciANIResult(dist=1.0, dist_low=1.0, dist_high=1.0, p_nothing_in_common=1.0, p_threshold=0.001)
    assert res2 == exp_res


def test_containment_to_distance_one():
    contain = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,estimate_ci=True)
    print(res)
    exp_dist, exp_low,exp_high,pnc = 0.0,0.0,0.0,0.0
    exp_id, exp_idlow,exp_idhigh,pnc = 1.0,1.0,1.0,0.0
    assert res.dist == exp_dist
    assert res.dist_low == exp_low
    assert res.dist_high == exp_high
    assert res.p_nothing_in_common == pnc
    assert res.ani == exp_id
    assert res.ani_low == exp_idlow
    assert res.ani_high == exp_idhigh

    # check without returning ci
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    assert res.dist == exp_dist
    assert res.ani == exp_id
    assert res.p_nothing_in_common == pnc
    assert res.ani_low == 1.0
    assert res.ani_high == 1.0


def test_containment_to_distance_scaled1():
    contain = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.032468221476108394
    assert res.ani == 0.9675317785238916
    assert res.dist_low == 0.028709912966405623
    assert res.ani_high == 0.9712900870335944
    assert res.dist_high == 0.03647860197289783
    assert res.ani_low == 0.9635213980271021
    assert res.p_nothing_in_common == 0.0
    # without returning ci
    res2 = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    assert (res2.dist,res2.ani,res2.p_nothing_in_common) == (0.032468221476108394, 0.9675317785238916, 0.0)
    assert (res2.dist,res2.ani,res2.p_nothing_in_common) == (res.dist, res.ani, res.p_nothing_in_common)


def test_containment_to_distance_scaled100():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.07158545548052564
    assert res.dist_low == 0.05320779238601372
    assert res.dist_high == 0.09055547672455365
    assert res.p_nothing_in_common == 4.3171247410658655e-05
    assert res.p_exceeds_threshold == False


def test_containment_to_distance_scaled100_2():
    contain = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    res= containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.032468221476108394
    assert res.dist_low == 0.023712063916639017
    assert res.dist_high == 0.04309960543965866
    assert res.p_exceeds_threshold == False


def test_containment_to_distance_k10():
    contain = 0.5
    scaled = 100
    nkmers = 10000
    ksize=10
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.06696700846319259
    assert res.dist_low == 0.04982777541057476
    assert res.dist_high == 0.08745108232411622
    assert res.p_exceeds_threshold == False


def test_containment_to_distance_confidence():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    confidence=0.99
    res = containment_to_distance(contain,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers, estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.07158545548052564
    assert res.dist_low == 0.04802880300938562
    assert res.dist_high == 0.09619930040790341
    assert res.p_exceeds_threshold == False
    confidence=0.90
    res2 = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,confidence=confidence, estimate_ci=True)
    print(res2)
    # check results
    assert res2.dist == res.dist
    assert res2.dist_low == 0.05599435479247415
    assert res2.dist_high == 0.08758718871990222
    assert res.p_exceeds_threshold == False


def test_nkmers_to_bp_containment():
    containment = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = handle_seqlen_nkmers(ksize, sequence_len_bp= bp_len)
    print("nkmers_from_bp:", nkmers)
    confidence=0.99
    kmer_res = containment_to_distance(containment,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers,estimate_ci=True)
    bp_res = containment_to_distance(containment,ksize,scaled,confidence=confidence,sequence_len_bp=bp_len,estimate_ci=True)
    print(f"\nkDIST: {kmer_res}")
    print(f"\nbpDIST:,{bp_res}")
    # check results
    assert kmer_res==bp_res
    assert kmer_res.dist == 0.07158545548052564
    assert kmer_res.dist_low == 0.04802880300938562
    assert kmer_res.dist_high == 0.09619930040790341


def test_jaccard_to_distance_zero():
    jaccard = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    res= jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    # check results
    assert res.dist == 1.0
    assert res.ani == 0.0
    assert res.p_nothing_in_common == 1.0
    assert res.jaccard_error == 0.0


def test_jaccard_to_distance_one():
    jaccard = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    res= jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    # check results
    assert res.dist == 0.0
    assert res.ani == 1.0
    assert res.p_nothing_in_common == 0.0
    assert res.jaccard_error == 0.0


def test_jaccard_to_distance_scaled():
    # scaled value doesn't impact point estimate or jaccard error, just p_nothing_in_common
    jaccard = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    res = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    # check results
    assert round(res.dist, 3) == round(0.019122659390482077, 3)
    assert res.ani == None
    assert res.p_exceeds_threshold == False
    assert res.jaccard_error == 0.00018351337045518042
    assert res.je_exceeds_threshold ==True
    scaled = 100
    res2 = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res2)
    assert res2.dist == res.dist
    assert res2.jaccard_error == res.jaccard_error
    assert res2.p_nothing_in_common != res.p_nothing_in_common
    assert res2.p_exceeds_threshold ==False


def test_jaccard_to_distance_k31():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=31
    res = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    # check results
    assert res.je_exceeds_threshold ==True
    assert res.ani == None
    assert res.p_exceeds_threshold == False
    res2 = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers, err_threshold=0.1)
    assert res2.je_exceeds_threshold == False
    assert res2.ani == 0.9870056455892898


def test_jaccard_to_distance_k31_2():
    jaccard = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    res = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    # check results
    assert res.ani == 0.9464928391768298
    assert res.p_exceeds_threshold == False
    assert res.je_exceeds_threshold == False


def test_nkmers_to_bp_jaccard():
    jaccard = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = handle_seqlen_nkmers(ksize, sequence_len_bp= bp_len)
    print("nkmers_from_bp:", nkmers)
    kmer_res = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    bp_res = jaccard_to_distance(jaccard,ksize,scaled,sequence_len_bp=bp_len)
    print(f"\nkmer_res: {kmer_res}")
    print(f"\nbp_res: {bp_res}")
    # check results
    assert kmer_res == bp_res
    assert kmer_res.dist == 0.0535071608231702
    assert kmer_res.p_exceeds_threshold == False
    assert kmer_res.je_exceeds_threshold == False


def test_exp_prob_nothing_common():
    dist = 0.25
    ksize = 31
    scaled = 10
    bp_len = 1000030
    nkmers = handle_seqlen_nkmers(ksize, sequence_len_bp= bp_len)
    print("nkmers_from_bp:", nkmers)

    nkmers_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,n_unique_kmers=nkmers)
    print(f"prob nothing in common: {nkmers_pnc}")
    bp_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,sequence_len_bp=bp_len)
    assert nkmers_pnc == bp_pnc == 7.437016945722123e-07


def test_containment_to_distance_tinytestdata_var0():
    """
    tiny test data to trigger the following:
    WARNING: Cannot estimate ANI confidence intervals from containment. Do your sketches contain enough hashes?
    Error: varN <0.0!
    """
    contain = 0.9
    scaled = 1
    nkmers = 4
    ksize=31
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers, estimate_ci=True)
    print(res)
    # check results
    assert res.dist == 0.003392957179023992
    assert res.dist_low == None
    assert res.dist_high == None
    assert res.ani_low == None
    assert res.ani_high == None
    assert res.p_exceeds_threshold == False


def test_var_n_mutated():
    # check 0
    r = 0
    ksize = 31
    nkmers = 200
    var_n_mut = var_n_mutated(nkmers,ksize,r)
    print(f"var_n_mutated: {var_n_mut}")
    assert var_n_mut == 0
    # check var 0.0 valuerror
    r = 10
    ksize = 31
    nkmers = 200
    with pytest.raises(ValueError) as exc:
        var_n_mut = var_n_mutated(nkmers,ksize,r)
    assert "Error: varN <0.0!" in str(exc)
    # check successful
    r = 0.4
    ksize = 31
    nkmers = 200000
    var_n_mut = var_n_mutated(nkmers,ksize,r)
    print(f"var_n_mutated: {var_n_mut}")
    assert var_n_mut == 0.10611425440741508


def test_handle_seqlen_nkmers():
    bp_len = 10030
    ksize=31
    # convert seqlen to nkmers
    nkmers = handle_seqlen_nkmers(ksize, sequence_len_bp= bp_len)
    assert nkmers == 10000
    # if nkmers is provided, just use that
    nkmers = handle_seqlen_nkmers(ksize, sequence_len_bp= bp_len, n_unique_kmers= bp_len)
    assert nkmers == 10030
    # if neither seqlen or nkmers provided, complain
    with pytest.raises(ValueError) as exc:
        nkmers = handle_seqlen_nkmers(ksize)
    assert("Error: distance estimation requires input of either 'sequence_len_bp' or 'n_unique_kmers'") in str(exc)


def test_set_size_chernoff():
    eps = 10**(-6)
    rel_error = 0.01
    set_size = 1000000
    s = 1/0.1  # I'm used to using a scale value between 0 and 1
    value_from_mathematica = 0.928652
    assert np.abs(set_size_chernoff(set_size, s, relative_error=rel_error) - value_from_mathematica) < eps

    rel_error = 0.05
    set_size = 10000
    s = 1
    value_from_mathematica = 0.999519
    assert np.abs(set_size_chernoff(set_size, s, relative_error=rel_error) - value_from_mathematica) < eps

    rel_error = 0.001
    set_size = 10
    s = 1/.01
    value_from_mathematica = -1
    assert np.abs(set_size_chernoff(set_size, s, relative_error=rel_error) - value_from_mathematica) < eps


def test_set_size_exact_prob():
    # values obtained from Mathematica
    # specifically: Probability[Abs[X*s - n]/n <= delta,
    #   X \[Distributed] BinomialDistribution[n, 1/s]] // N
    set_size = 100
    scaled = 2
    relative_error = 0.05
    prob = set_size_exact_prob(set_size, scaled, relative_error=relative_error)
    true_prob = 0.382701
    np.testing.assert_array_almost_equal(true_prob, prob, decimal=3)

    set_size = 200
    scaled = 5
    relative_error = 0.15
    prob = set_size_exact_prob(set_size, scaled, relative_error=relative_error)
    true_prob = 0.749858
    np.testing.assert_array_almost_equal(true_prob, prob, decimal=3)

    set_size = 10
    scaled = 10
    relative_error = 0.10
    prob = set_size_exact_prob(set_size, scaled, relative_error=relative_error)
    true_prob = 0.38742
    np.testing.assert_array_almost_equal(true_prob, prob, decimal=3)

    set_size = 1000
    scaled = 10
    relative_error = 0.10
    prob = set_size_exact_prob(set_size, scaled, relative_error=relative_error)
    true_prob = 0.73182
    np.testing.assert_array_almost_equal(true_prob, prob, decimal=3)
