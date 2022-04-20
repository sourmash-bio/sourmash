"""
Tests for the 'SketchComparison' classes.
"""

import numpy as np
import pytest


from sourmash.minhash import MinHash
from sourmash.sketchcomparison import BaseMinHashComparison, FracMinHashComparison, NumMinHashComparison

import sourmash_tst_utils as utils

# can we parameterize scaled too (so don't need separate downsample tests?)
def test_FracMinHashComparison(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    # build FracMinHashComparison 
    cmp = FracMinHashComparison(a, b)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_scaled == 1
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.mh1_containment == a.contained_by(b)
    assert cmp.mh2_containment == b.contained_by(a)
    assert cmp.avg_containment == np.mean([a.contained_by(b), b.contained_by(a)])
    assert cmp.max_containment == a.max_containment(b)
    assert cmp.jaccard == a.jaccard(b) == b.jaccard(a)
    intersect_mh = a.flatten().intersection(b.flatten())
    assert cmp.intersect_mh == intersect_mh == b.flatten().intersection(a.flatten())
    assert cmp.intersect_bp == 4
    assert cmp.pass_threshold # default threshold is 0; this should pass
    if track_abundance:
        assert cmp.angular_similarity == a.angular_similarity(b) == b.angular_similarity(a)
        assert cmp.cosine_similarity == a.angular_similarity(b) == b.angular_similarity(a)
        assert cmp.weighted_intersection(from_mh=cmp.mh1).hashes == intersect_mh.inflate(a).hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2).hashes == intersect_mh.inflate(b).hashes
        assert cmp.weighted_intersection(from_abundD=a_values).hashes == intersect_mh.inflate(a).hashes
        assert cmp.weighted_intersection(from_abundD=b_values).hashes == intersect_mh.inflate(b).hashes
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        assert cmp.weighted_intersection(from_mh=cmp.mh1).hashes == intersect_mh.hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2).hashes == intersect_mh.hashes
    

def test_FracMinHashComparison_downsample(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    cmp_scaled = 2
    ds_a = a.downsample(scaled=cmp_scaled)
    ds_b = b.downsample(scaled=cmp_scaled)

    # build FracMinHashComparison 
    cmp = FracMinHashComparison(a, b, cmp_scaled = cmp_scaled)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.mh1_cmp == ds_a
    assert cmp.mh2_cmp == ds_b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_scaled == cmp_scaled
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.mh1_containment == ds_a.contained_by(ds_b)
    assert cmp.mh2_containment == ds_b.contained_by(ds_a)
    assert cmp.avg_containment == np.mean([ds_a.contained_by(ds_b), ds_b.contained_by(ds_a)])
    assert cmp.max_containment == ds_a.max_containment(ds_b)
    assert cmp.jaccard == ds_a.jaccard(ds_b) == ds_b.jaccard(ds_a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    assert cmp.intersect_bp == 8
    assert cmp.pass_threshold # default threshold is 0; this should pass
    if track_abundance:
        assert cmp.angular_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.cosine_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.weighted_intersection(from_mh=cmp.mh1_cmp).hashes == intersect_mh.inflate(ds_a).hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2_cmp).hashes == intersect_mh.inflate(ds_b).hashes
        assert cmp.weighted_intersection(from_abundD=cmp.mh1_cmp.hashes).hashes == intersect_mh.inflate(ds_a).hashes
        assert cmp.weighted_intersection(from_abundD=cmp.mh2_cmp.hashes).hashes == intersect_mh.inflate(ds_b).hashes
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        assert cmp.weighted_intersection(from_mh=cmp.mh1_cmp).hashes == intersect_mh.hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2_cmp).hashes == intersect_mh.hashes


def test_FracMinHashComparison_autodownsample(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=2, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    cmp_scaled = 2
    ds_a = a.downsample(scaled=cmp_scaled)
    ds_b = b.downsample(scaled=cmp_scaled)

    # build FracMinHashComparison 
    cmp = FracMinHashComparison(a, b)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.mh1_cmp == ds_a
    assert cmp.mh2_cmp == ds_b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_scaled == cmp_scaled
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.mh1_containment == ds_a.contained_by(ds_b)
    assert cmp.mh2_containment == ds_b.contained_by(ds_a)
    assert cmp.avg_containment == np.mean([ds_a.contained_by(ds_b), ds_b.contained_by(ds_a)])
    assert cmp.max_containment == ds_a.max_containment(ds_b)
    assert cmp.jaccard == ds_a.jaccard(ds_b) == ds_b.jaccard(ds_a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    assert cmp.intersect_bp == 8
    assert cmp.pass_threshold # default threshold is 0; this should pass
    if track_abundance:
        assert cmp.angular_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.cosine_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.weighted_intersection(from_mh=cmp.mh1_cmp).hashes == intersect_mh.inflate(ds_a).hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2_cmp).hashes == intersect_mh.inflate(ds_b).hashes
        assert cmp.weighted_intersection(from_abundD=a_values).hashes == intersect_mh.inflate(a).hashes
        assert cmp.weighted_intersection(from_abundD=b_values).hashes == intersect_mh.inflate(b).hashes
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        assert cmp.weighted_intersection(from_mh=cmp.mh1_cmp).hashes == intersect_mh.hashes
        assert cmp.weighted_intersection(from_mh=cmp.mh2_cmp).hashes == intersect_mh.hashes


def test_FracMinHashComparison_ignore_abundance(track_abundance):
    # build FracMinHash Comparison and check values
    a = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    intersection_w_abund = {1:8, 3:5, 5:3, 8:3}
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())
    
    cmp_scaled = 2
    ds_a = a.flatten().downsample(scaled=cmp_scaled)
    ds_b = b.flatten().downsample(scaled=cmp_scaled)

    # build FracMinHashComparison 
    cmp = FracMinHashComparison(a, b, cmp_scaled = cmp_scaled, ignore_abundance=True)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == True
    assert cmp.cmp_scaled == cmp_scaled
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.mh1_containment == a.contained_by(b)
    assert cmp.mh2_containment == b.contained_by(a)
    assert cmp.avg_containment == np.mean([a.contained_by(b), b.contained_by(a)])
    assert cmp.max_containment == a.max_containment(b)
    assert cmp.jaccard == a.jaccard(b) == b.jaccard(a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    assert cmp.intersect_bp == 8
    assert cmp.pass_threshold # default threshold is 0; this should pass
    # with ignore_abundance = True, all of these should not be usable. Do we want errors, or ""/None?
    with pytest.raises(TypeError) as exc:
        cmp.angular_similarity
    print(str(exc))
    assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
    with pytest.raises(TypeError) as exc:
        cmp.cosine_similarity
    print(str(exc))
    assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
    assert not cmp.mh1_cmp.track_abundance
    assert not cmp.mh2_cmp.track_abundance
    assert cmp.weighted_intersection(from_mh=cmp.mh1_cmp).hashes == intersect_mh.hashes
    assert cmp.weighted_intersection(from_mh=cmp.mh2_cmp).hashes == intersect_mh.hashes


def test_FracMinHashComparison_fail_threshold(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(0, 21, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=1, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    cmp_scaled = 2
    ds_a = a.flatten().downsample(scaled=cmp_scaled)
    ds_b = b.flatten().downsample(scaled=cmp_scaled)

    # build FracMinHashComparison
    cmp = FracMinHashComparison(a, b, cmp_scaled = cmp_scaled, threshold_bp=10)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_scaled == cmp_scaled
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.mh1_containment == a.contained_by(b)
    assert cmp.mh2_containment == b.contained_by(a)
    assert cmp.avg_containment == np.mean([a.contained_by(b), b.contained_by(a)])
    assert cmp.max_containment == a.max_containment(b)
    assert cmp.jaccard == a.jaccard(b) == b.jaccard(a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    assert cmp.intersect_bp == 8
    assert not cmp.pass_threshold # threshold is 10; this should fail


def test_FracMinHashComparison_incompatible_ksize(track_abundance):
    a = MinHash(0, 31, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 21, scaled=2, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    with pytest.raises(TypeError) as exc:
        FracMinHashComparison(a, b)
    print(str(exc))
    assert "Error: Invalid Comparison, ksizes: 31, 21. Must compare sketches of the same ksize." in str(exc)


def test_FracMinHashComparison_incompatible_moltype(track_abundance):
    a = MinHash(0, 31, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 31, scaled=2, is_protein=True, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    with pytest.raises(TypeError) as exc:
        FracMinHashComparison(a, b)
    print(str(exc))
    assert "Error: Invalid Comparison, moltypes: DNA, protein. Must compare sketches of the same moltype." in str(exc)



def test_FracMinHashComparison_incompatible_sketchtype(track_abundance):
    a = MinHash(0, 31, scaled=1, track_abundance=track_abundance)
    b = MinHash(10, 31, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    with pytest.raises(TypeError) as exc:
        FracMinHashComparison(a, b)
    print(str(exc))
    assert "Error: Both sketches must be 'num' or 'scaled'." in str(exc)


def test_FracMinHashComparison_redownsample_without_scaled(track_abundance):
    a = MinHash(0, 31, scaled=1, track_abundance=track_abundance)
    b = MinHash(0, 31, scaled=10, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    cmp = FracMinHashComparison(a, b)
    assert cmp.cmp_scaled == 10

    with pytest.raises(ValueError) as exc:
        # try to redownsample without passing in cmp_num
        cmp.downsample_and_handle_ignore_abundance()
    print(str(exc))
    assert "Error: must pass in a comparison scaled or num value." in str(exc)


def test_NumMinHashComparison(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(10, 21, scaled=0, track_abundance=track_abundance)
    b = MinHash(10, 21, scaled=0, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())
    
    assert a.num and b.num and not a.scaled and not b.scaled
    
    # build NumMinHashComparison 
    cmp = NumMinHashComparison(a, b)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_num == 10
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.jaccard == a.jaccard(b) == b.jaccard(a)
    intersect_mh = a.flatten().intersection(b.flatten())
    assert cmp.intersect_mh == intersect_mh == b.flatten().intersection(a.flatten())
    if track_abundance:
        assert cmp.angular_similarity == a.angular_similarity(b) == b.angular_similarity(a)
        assert cmp.cosine_similarity == a.angular_similarity(b) == b.angular_similarity(a)
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)


def test_NumMinHashComparison_downsample(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(10, 21, scaled=0, track_abundance=track_abundance)
    b = MinHash(10, 21, scaled=0, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())
    
    assert a.num and b.num and not a.scaled and not b.scaled

    cmp_num = 5
    ds_a = a.downsample(num=cmp_num)
    ds_b = b.downsample(num=cmp_num)
    # build NumMinHashComparison 
    cmp = NumMinHashComparison(a, b, cmp_num = cmp_num)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_num == cmp_num
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.jaccard == ds_a.jaccard(ds_b) == ds_b.jaccard(ds_a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    if track_abundance:
        assert cmp.angular_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.cosine_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)


def test_NumMinHashComparison_autodownsample(track_abundance):
    # build FracMinHash Comparison and check values 
    a = MinHash(10, 21, scaled=0, track_abundance=track_abundance)
    b = MinHash(5, 21, scaled=0, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())
    
    assert a.num and b.num and not a.scaled and not b.scaled

    cmp_num = 5
    ds_a = a.downsample(num=cmp_num)
    ds_b = b.downsample(num=cmp_num)
    # build NumMinHashComparison 
    cmp = NumMinHashComparison(a, b)
    assert cmp.mh1 == a
    assert cmp.mh2 == b
    assert cmp.ignore_abundance == False
    assert cmp.cmp_num == cmp_num
    assert cmp.ksize == 21
    assert cmp.moltype == "DNA"
    assert cmp.jaccard == ds_a.jaccard(ds_b) == ds_b.jaccard(ds_a)
    intersect_mh = ds_a.flatten().intersection(ds_b.flatten())
    assert cmp.intersect_mh == intersect_mh == ds_b.flatten().intersection(ds_a.flatten())
    if track_abundance:
        assert cmp.angular_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
        assert cmp.cosine_similarity == ds_a.angular_similarity(ds_b) == ds_b.angular_similarity(ds_a)
    else:
        with pytest.raises(TypeError) as exc:
            cmp.angular_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)
        with pytest.raises(TypeError) as exc:
            cmp.cosine_similarity
        print(str(exc))
        assert "Error: Angular (cosine) similarity requires both sketches to track hash abundance." in str(exc)


def test_NumMinHashComparison_incompatible_ksize(track_abundance):
    a_num = MinHash(20, 31, track_abundance=track_abundance)
    b_num = MinHash(10, 21, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a_num.set_abundances(a_values)
        b_num.set_abundances(b_values)
    else:
        a_num.add_many(a_values.keys())
        b_num.add_many(b_values.keys())

    # build NumMinHashComparison
    with pytest.raises(TypeError) as exc:
        NumMinHashComparison(a_num, b_num)
    print(str(exc))
    assert "Error: Invalid Comparison, ksizes: 31, 21. Must compare sketches of the same ksize." in str(exc)


def test_NumMinHashComparison_incompatible_moltype(track_abundance):
    a_num = MinHash(20, 31, track_abundance=track_abundance)
    b_num = MinHash(10, 31, is_protein=True, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a_num.set_abundances(a_values)
        b_num.set_abundances(b_values)
    else:
        a_num.add_many(a_values.keys())
        b_num.add_many(b_values.keys())

    with pytest.raises(TypeError) as exc:
        NumMinHashComparison(a_num, b_num)
    print(str(exc))
    assert "Error: Invalid Comparison, moltypes: DNA, protein. Must compare sketches of the same moltype." in str(exc)


def test_NumMinHashComparison_incompatible_sketchtype(track_abundance):
    a = MinHash(0, 31, scaled=1, track_abundance=track_abundance)
    b = MinHash(10, 31, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    with pytest.raises(TypeError) as exc:
        NumMinHashComparison(a, b)
    print(str(exc))
    assert "Error: Both sketches must be 'num' or 'scaled'." in str(exc)


def test_NumMinHashComparison_redownsample_without_num(track_abundance):
    a = MinHash(10, 31, track_abundance=track_abundance)
    b = MinHash(5, 31, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }

    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    cmp = NumMinHashComparison(a, b)

    with pytest.raises(ValueError) as exc:
        # try to redownsample without passing in cmp_num
        cmp.downsample_and_handle_ignore_abundance()
    print(str(exc))
    assert "Error: must pass in a comparison scaled or num value." in str(exc)
