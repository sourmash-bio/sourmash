"""
Tests for distance utils.
"""
import pytest
from sourmash.distance_utils import (containment_to_distance, get_exp_probability_nothing_common,
                                    handle_seqlen_nkmers, jaccard_to_distance,
                                    ANIResult, ciANIResult, jaccardANIResult)

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
    assert res.ani == 0.6
    assert res.p_nothing_in_common == 0.1
    assert res.jaccard_error == 0.03
    assert res.p_exceeds_threshold ==True
    assert res.je_exceeds_threshold ==True
    res2 = jaccardANIResult(0.4, 0.1, jaccard_error=0.03, je_threshold=0.1)
    assert res2.je_exceeds_threshold ==False


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
    res = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print(res)
    exp_res = ciANIResult(dist=1.0, p_nothing_in_common=1.0, p_threshold=0.001)
    assert res == exp_res


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
    nkmers = handle_seqlen_nkmers(bp_len,ksize)
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
    assert res.dist == 0.019122659390482077
    assert res.ani == 0.9808773406095179
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
    assert res.ani == 0.9870056455892898
    assert res.p_exceeds_threshold == False
    assert res.je_exceeds_threshold ==True
    res2 = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers, err_threshold=0.1)
    assert res2.ani == res.ani
    assert res2.je_exceeds_threshold == False


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
    nkmers = handle_seqlen_nkmers(bp_len,ksize)
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
    nkmers = handle_seqlen_nkmers(bp_len,ksize)
    print("nkmers_from_bp:", nkmers)

    nkmers_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,n_unique_kmers=nkmers)
    print(f"prob nothing in common: {nkmers_pnc}")
    bp_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,sequence_len_bp=bp_len)
    assert nkmers_pnc == bp_pnc == 7.437016945722123e-07
