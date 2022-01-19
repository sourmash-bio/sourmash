"""
Tests for distance utils.
"""

import pytest
from sourmash.distance_utils import containment_to_distance, jaccard_to_distance, distance_to_identity, sequence_len_to_n_kmers

def test_distance_to_identity():
    id,id_low,id_high = distance_to_identity(0.5,0.4,0.6)
    assert id == 0.5
    assert id_low == 0.4
    assert id_high ==0.6


def test_distance_to_identity_fail():
    with pytest.raises(Exception) as exc:
        id,id_low,id_high = distance_to_identity(1.1,0.4,0.6)
    assert "distance value 1.1 is not between 0 and 1!" in str(exc.value)
    with pytest.raises(Exception) as exc:
        id,id_low,id_high = distance_to_identity(-0.1,0.4,0.6)
    assert "distance value -0.1 is not between 0 and 1!" in str(exc.value)


def test_containment_to_distance_zero():
    contain = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,ksize,scaled, n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist, exp_low,exp_high = 1.0,1.0,1.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 0.0,0.0,0.0
    ident,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_containment_to_distance_one():
    contain = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist, exp_low,exp_high = 0.0,0.0,0.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 1.0,1.0,1.0
    ident,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_containment_to_distance_scaled1():
    contain = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.032468221476108394,0.028709912966405623,0.03647860197289783
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 0.9675317785238916,0.9635213980271021,0.9712900870335944
    ident,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_containment_to_distance_2():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    dist,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.07158545548052564,0.05320779238601372,0.09055547672455365
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_containment_to_distance_scaled100():
    contain = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.032468221476108394,0.023712063916639017,0.04309960543965866
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_containment_to_distance_k10():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=10
    dist,low,high = containment_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.06696700846319259,0.04982777541057476,0.08745108232411622
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high

def test_containment_to_distance_confidence():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    confidence=0.99
    dist,low,high = containment_to_distance(contain,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.07158545548052564,0.04802880300938562,0.09619930040790341
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high

    confidence=0.90
    dist,low,high = containment_to_distance(contain,ksize,scaled,n_unique_kmers=nkmers,confidence=confidence)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.07158545548052564,0.05599435479247415,0.08758718871990222
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_nkmers_to_bp_containment():
    containment = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = sequence_len_to_n_kmers(bp_len,ksize)
    print("nkmers_from_bp:", bp_len)
    confidence=0.99
    kmer_dist = containment_to_distance(containment,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers)
    bp_dist = containment_to_distance(containment,ksize,scaled,confidence=confidence,sequence_len_bp=bp_len)
    print(f"\nkDIST:", kmer_dist)
    print(f"\nbpDIST:",bp_dist)
    # check results
    assert kmer_dist == (0.07158545548052564, 0.04802880300938562, 0.09619930040790341)
    assert bp_dist ==  (0.07158545548052564, 0.04802880300938562, 0.09619930040790341)
    assert kmer_dist==bp_dist


def test_jaccard_to_distance_zero():
    jaccard = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 1.0,1.0,1.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 0.0,0.0,0.0
    ident,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_jaccard_to_distance_one():
    jaccard = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0,0.0,0.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 1.0,1.0,1.0
    ident,low,high = jaccard_to_distance(jaccard,ksize,scaled,return_identity=True,n_unique_kmers=nkmers)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_jaccard_to_distance_scaled1():
    jaccard = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.019122659390482077,0.017549816764592427,0.02085752009333101
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high
    # return identity instead
    exp_id, exp_idlow,exp_idhigh = 0.9808773406095179,0.979142479906669,0.9824501832354076
    ident,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers,return_identity=True)
    print(f"{ident},{low},{high}")
    assert ident == exp_id
    assert low == exp_idlow
    assert high == exp_idhigh


def test_jaccard_to_distance_scaled100():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.019122659390482077,0.014156518319318126,0.025095471100086125
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_k31():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=31
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.012994354410710174,0.009559840649526866,0.017172145712325716
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_k31_2():
    jaccard = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0535071608231702,0.040739632227821516,0.06673746391115623
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_confidence():
    jaccard = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    confidence=0.99
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0535071608231702,0.03702518586582857,0.07080999238232429
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high

    confidence=0.90
    dist,low,high = jaccard_to_distance(jaccard,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0535071608231702,0.042708361726101415,0.0646280650023921
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_nkmers_to_bp():
    jaccard = 0.1
    scaled = 100
    bp_len = 10030
    ksize=31
    nkmers = sequence_len_to_n_kmers(bp_len,ksize)
    print("nkmers_from_bp:", bp_len)
    confidence=0.99
    kmer_dist = jaccard_to_distance(jaccard,ksize,scaled,confidence=confidence,n_unique_kmers=nkmers)
    bp_dist = jaccard_to_distance(jaccard,ksize,scaled,confidence=confidence,sequence_len_bp=bp_len)
    print(f"\nkDIST:", kmer_dist)
    print(f"\nbpDIST:",bp_dist)
    # check results
    assert kmer_dist == (0.0535071608231702, 0.03702518586582857, 0.07080999238232429)
    assert bp_dist ==  (0.0535071608231702, 0.03702518586582857, 0.07080999238232429)
    assert kmer_dist==bp_dist
