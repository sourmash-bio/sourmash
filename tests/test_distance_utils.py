"""
Tests for distance utils.
"""
import pytest
from sourmash.distance_utils import containment_to_distance, get_exp_probability_nothing_common, \
                                    jaccard_to_distance, sequence_len_to_n_kmers, ANIResult, \
                                    ciANIResult, jaccardANIResult

def test_aniresult():
    res = ANIResult(0.4, 0.1)
    assert res.dist == 0.4
    assert res.ani == 0.6
    assert res.p_nothing_in_common == 0.1

def test_aniresult_bad_distance():
    """
    Fail if distance is not between 0 and 1.
    """
    with pytest.raises(Exception) as exc:
        ANIResult(1.1, 0.1)
    print("\n", str(exc.value))
    assert "distance value 1.1 is not between 0 and 1!" in str(exc.value)
    with pytest.raises(Exception) as exc:
        jaccardANIResult(1.1, 0.1)
    print("\n", str(exc.value))
    assert "distance value 1.1 is not between 0 and 1!" in str(exc.value)

def test_jaccard_aniresult():
    res = jaccardANIResult(0.4, 0.1, jaccard_error=0.003)
    assert res.dist == 0.4
    assert res.ani == 0.6
    assert res.p_nothing_in_common == 0.1
    assert res.jaccard_error == 0.003

#def test_distance_to_identity():
#    ident,id_low,id_high = distance_to_identity(0.5,d_low=0.4,d_high=0.6)
#    assert ident == 0.5
#    assert id_low == 0.4
#    assert id_high ==0.6


#def test_distance_to_identity_just_point_estimate():
#    ident = distance_to_identity(0.4)
#    assert ident == 0.6


#def test_distance_to_identity_fail_distance():
#    """
#    Fail if distance is not between 0 and 1.
#    """
#    with pytest.raises(Exception) as exc:
#        ident,id_low,id_high = distance_to_identity(1.1,d_low=0.4,d_high=0.6)
#    assert "distance value 1.1 is not between 0 and 1!" in str(exc.value)
#    with pytest.raises(Exception) as exc:
#        ident,id_low,id_high = distance_to_identity(-0.1,d_low=0.4,d_high=0.6)
#    assert "distance value -0.1 is not between 0 and 1!" in str(exc.value)


#def test_distance_to_identity_fail_incorrect_input():
#    """
#    Ignore 2nd input, unless put both dist_low and dist_high
#    """
#    with pytest.raises(Exception) as exc: 
#        distance_to_identity(0.4,d_low=0.5)
#    print("\n", str(exc.value))
#    assert "Error: `distance_to_identity` expects either one value (distance) or three values" in str(exc.value)


def test_containment_to_distance_zero():
    contain = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high,p_not_in_common= containment_to_distance(contain,ksize,scaled, n_unique_kmers=nkmers, return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist,exp_low,exp_high,pnc = 1.0,1.0,1.0,1.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert p_not_in_common == pnc
    # check without returning ci
    result = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    assert dist == exp_dist
    assert p_not_in_common == pnc
    assert result.dist == exp_dist
    assert result.dist_low == None
    assert result.dist_high == None
    assert result.p_nothing_in_common == pnc
    # return identity instead
    exp_id, exp_idlow,exp_idhigh,pnc = 0.0,0.0,0.0,1.0
    result = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True, return_ci=True)
    assert result.ani == exp_id
    assert result.ci_low == exp_idlow
    assert result.ci_high == exp_idhigh
    assert result.p_nothing_in_common == pnc


def test_containment_to_distance_one():
    contain = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.0,0.0,0.0,0.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert p_not_in_common == pnc
    # check without returning ci
    dist,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    assert dist == exp_dist
    assert p_not_in_common == pnc
    # return identity instead
    exp_id, exp_idlow,exp_idhigh,pnc = 1.0,1.0,1.0,0.0
    ident,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True, return_ci=True)
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh
    assert p_not_in_common == pnc


def test_containment_to_distance_scaled1():
    contain = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high},{p_not_in_common}")
    # check results
    assert (dist,low,high,p_not_in_common) == (0.032468221476108394,0.028709912966405623,0.03647860197289783,0.0)
    # without returning ci
    dist,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    assert (dist,p_not_in_common) == (0.032468221476108394,0.0)
    # return identity instead
    ident,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True,return_ci=True)
    print(f"{ident},{low},{high}")
    assert (ident,low,high,p_not_in_common) == (0.9675317785238916,0.9635213980271021,0.9712900870335944,0.0)
    # without returning ci
    ident,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    assert (ident,p_not_in_common) == (0.9675317785238916,0.0)


def test_containment_to_distance_2():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high},{p_not_in_common}")
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.07158545548052564,0.05320779238601372,0.09055547672455365,4.3171247410658655e-05
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert p_not_in_common == pnc


def test_containment_to_distance_scaled100():
    contain = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.032468221476108394,0.023712063916639017,0.04309960543965866,1.4995915609979772e-22
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert pnc == p_not_in_common


def test_containment_to_distance_k10():
    contain = 0.5
    scaled = 100
    nkmers = 10000
    ksize=10
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high},{p_not_in_common}")
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.06696700846319259,0.04982777541057476,0.08745108232411622,1.499591560997956e-22
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert pnc == p_not_in_common


def test_containment_to_distance_confidence():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    confidence=0.99
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers, return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high},{p_not_in_common}")
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.07158545548052564,0.04802880300938562,0.09619930040790341,4.3171247410658655e-05
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert pnc == p_not_in_common
    confidence=0.90
    dist,low,high,p_not_in_common = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,confidence=confidence, return_ci=True)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high},{p_not_in_common}")
    # check results
    exp_dist, exp_low,exp_high,pnc = 0.07158545548052564,0.05599435479247415,0.08758718871990222,4.3171247410658655e-05
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    assert pnc == p_not_in_common


def test_nkmers_to_bp_containment():
    containment = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = sequence_len_to_n_kmers(bp_len,ksize)
    print("nkmers_from_bp:", bp_len)
    confidence=0.99
    kmer_dist = containment_to_distance(containment,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers,return_ci=True)
    bp_dist = containment_to_distance(containment,ksize,scaled,confidence=confidence,sequence_len_bp=bp_len,return_ci=True)
    print(f"\nkDIST: {kmer_dist}")
    print(f"\nbpDIST:,{bp_dist}")
    # check results
    assert kmer_dist == (0.07158545548052564, 0.04802880300938562, 0.09619930040790341,4.3171247410658655e-05)
    assert bp_dist ==  (0.07158545548052564, 0.04802880300938562, 0.09619930040790341,4.3171247410658655e-05)
    assert kmer_dist==bp_dist



def test_jaccard_to_distance_zero():
    jaccard = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist, p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(f"{dist},{p_not_in_common},{errlb}")
    # check results
    exp_dist, pnc,err = 1.0,1.0,0.0
    assert dist == exp_dist
    assert p_not_in_common == pnc
    assert errlb == err
    # return identity instead
    exp_id=0.0
    ident,low,high = jaccard_to_distance(jaccard,ksize,scaled,return_identity=True,n_unique_kmers=nkmers)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert p_not_in_common == pnc
    assert errlb == err


def test_jaccard_to_distance_one():
    jaccard = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist, p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(f"{dist},{p_not_in_common},{errlb}")
    # check results
    exp_dist, pnc,err = 0.0,0.0,0.0
    assert dist == exp_dist
    assert p_not_in_common == pnc
    assert errlb == err
    # return identity instead
    exp_id=1.0
    ident,low,high = jaccard_to_distance(jaccard,ksize,scaled,return_identity=True,n_unique_kmers=nkmers)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert p_not_in_common == pnc
    assert errlb == err


def test_jaccard_to_distance_scaled1():
    jaccard = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist, p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print(f"{dist},{p_not_in_common},{errlb}")
    # check results
    assert (dist,p_not_in_common,errlb) == (0.019122659390482077, 0.0, 0.00018351337045518042) 
    # return identity instead
    ident, p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    print(f"{ident},{p_not_in_common},{errlb}")
    assert (ident,p_not_in_common, errlb) == (0.9808773406095179, 0.0, 0.00018351337045518042)


def test_jaccard_to_distance_scaled100():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    ident,p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers, return_identity=True)
    print(f"{ident},{p_not_in_common},{errlb}")
    # check results
    assert (ident,p_not_in_common, errlb) == (0.9808773406095179, 7.967045858822289e-30, 0.00018351337045518042) 


def test_jaccard_to_distance_k31():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=31
    ident,p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers, return_identity=True)
    print(f"{ident},{p_not_in_common},{errlb}")
    # check results
    assert (ident,p_not_in_common, errlb) == (0.9870056455892898, 7.967045858822514e-30, 0.00027078924329882424) 


def test_jaccard_to_distance_k31_2():
    jaccard = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    ident,p_not_in_common, errlb = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers, return_identity=True)
    print(f"{ident},{p_not_in_common},{errlb}")
    # check results
    assert (ident,p_not_in_common, errlb) == (0.9464928391768298, 1.1587511478551202e-08, 5.588721845727948e-05) 


def test_nkmers_to_bp_jaccard():
    jaccard = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = sequence_len_to_n_kmers(bp_len,ksize)
    print("nkmers_from_bp:", nkmers)
    kmer_dist = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    bp_dist = jaccard_to_distance(jaccard,ksize,scaled,sequence_len_bp=bp_len)
    print(f"\nk_jDIST: {kmer_dist}")
    print(f"\nbp_jDIST: {bp_dist}")
    # check results
    assert kmer_dist == (0.0535071608231702, 1.1587511478551202e-08, 5.588721845727948e-05)
    assert bp_dist ==  (0.0535071608231702, 1.1587511478551202e-08, 5.588721845727948e-05)
    assert kmer_dist==bp_dist


def test_exp_prob_nothing_common():
    dist = 0.25
    ksize = 31
    scaled = 10
    bp_len = 1000030
    nkmers = sequence_len_to_n_kmers(bp_len,ksize)
    print("nkmers_from_bp:", nkmers)

    nkmers_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,n_unique_kmers=nkmers)
    print(f"prob nothing in common: {nkmers_pnc}")
    bp_pnc = get_exp_probability_nothing_common(dist,ksize,scaled,sequence_len_bp=bp_len)
    assert nkmers_pnc == bp_pnc == 7.437016945722123e-07
