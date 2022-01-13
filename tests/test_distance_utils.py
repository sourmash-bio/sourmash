"""
Tests for distance utils.
"""
from sourmash.distance_utils import containment_to_distance, jaccard_to_distance


def test_containment_to_distance_zero():
    contain = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist, exp_low,exp_high = 1.0,1.0,1.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_containment_to_distance_one():
    contain = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    # check results
    exp_dist, exp_low,exp_high = 0.0,0.0,0.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_containment_to_distance_scaled1():
    contain = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.032468221476108394,0.028709912966405623,0.03647860197289783
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_containment_to_distance_2():
    contain = 0.1
    scaled = 100
    nkmers = 10000
    ksize=31
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled)
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
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled)
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
    dist,low,high = containment_to_distance(jaccard,nkmers,ksize,scaled)
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
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled,confidence)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.07158545548052564,0.04802880300938562,0.09619930040790341
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high

    confidence=0.90
    dist,low,high = containment_to_distance(contain,nkmers,ksize,scaled,confidence)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.07158545548052564,0.05599435479247415,0.08758718871990222
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_zero():
    jaccard = 0
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 1.0,1.0,1.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_one():
    jaccard = 1
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0,0.0,0.0
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_scaled1():
    jaccard = 0.5
    scaled = 1
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.019122659390482077,0.017549816764592427,0.02085752009333101
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high


def test_jaccard_to_distance_scaled100():
    jaccard = 0.5
    scaled = 100
    nkmers = 10000
    ksize=21
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
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
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
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
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled)
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
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled,confidence)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0535071608231702,0.03702518586582857,0.07080999238232429
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high

    confidence=0.90
    dist,low,high = jaccard_to_distance(jaccard,nkmers,ksize,scaled,confidence)
    print("\nDIST:", dist)
    print("CI:", low, " - ", high)
    print(f"{dist},{low},{high}")
    # check results
    exp_dist, exp_low,exp_high = 0.0535071608231702,0.042708361726101415,0.0646280650023921
    assert dist == exp_dist
    assert low == exp_low
    assert high == exp_high