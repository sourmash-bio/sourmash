import os

import pytest

import sourmash
from sourmash.signature import SourmashSignature, save_signatures, \
    load_signatures, load_one_signature, FrozenSourmashSignature
import sourmash_tst_utils as utils
from sourmash.minhash import MinHash, FrozenMinHash
from sourmash_tst_utils import SourmashCommandFailed


def test_minhash_copy(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig = SourmashSignature(e, name='foo')
    f = e.copy()
    assert e == f


def test_sig_copy(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')
    sig2 = sig1.copy()
    assert sig1 == sig2


def test_sig_copy_frozen(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')
    sig2 = sig1.copy()
    assert sig1 == sig2
    with pytest.raises(TypeError) as e:
        sig2.minhash.add_hash(5)
    assert 'FrozenMinHash does not support modification' in str(e.value)


def test_sig_copy_frozen_mutable(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')
    sig1.minhash = sig1.minhash.to_mutable()
    sig2 = sig1.copy()
    assert sig1 == sig2
    with pytest.raises(TypeError) as e:
        sig2.minhash.add_hash(5)
    assert 'FrozenMinHash does not support modification' in str(e.value)


def test_compare(track_abundance):
    # same content, same name -> equal
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')

    f = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    f.add_kmer("AT" * 10)
    sig2 = SourmashSignature(f, name='foo')

    assert e == f


def test_compare_ne(track_abundance):
    # same content, different names -> different
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')

    f = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    f.add_kmer("AT" * 10)
    sig2 = SourmashSignature(f, name='bar')

    assert sig1 != sig2


def test_compare_ne2(track_abundance):
    # same content, different filename -> different
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo', filename='a')

    f = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    f.add_kmer("AT" * 10)
    sig2 = SourmashSignature(f, name='foo', filename='b')

    assert sig1 != sig2
    assert sig2 != sig1


def test_compare_ne2_reverse(track_abundance):
    # same content, one has filename, other does not -> different
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig1 = SourmashSignature(e, name='foo')

    f = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    f.add_kmer("AT" * 10)
    sig2 = SourmashSignature(f, filename='b')

    assert sig2 != sig1
    assert sig1 != sig2


def test_hashable(track_abundance):
    # check: can we use signatures as keys in dictionaries and sets?
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)

    sig = SourmashSignature(e)

    x = set()
    x.add(sig)


def test_str(track_abundance):
    # signatures should be printable
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)

    sig = SourmashSignature(e)

    print(sig)
    assert repr(sig) == "SourmashSignature('', 59502a74)"

    sig._name = 'fizbar'
    assert repr(sig) == 'SourmashSignature(\'fizbar\', 59502a74)'


def test_roundtrip(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig = SourmashSignature(e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0
    assert isinstance(sig, SourmashSignature)
    assert not isinstance(sig, FrozenSourmashSignature)
    assert isinstance(sig2, FrozenSourmashSignature)

    assert isinstance(e, MinHash)
    assert isinstance(sig.minhash, FrozenMinHash)
    assert isinstance(sig2.minhash, FrozenMinHash)


def test_roundtrip_mutable_frozen(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig = SourmashSignature(e)
    assert isinstance(sig.minhash, FrozenMinHash)
    sig.minhash = sig.minhash.to_mutable()

    sig2 = sig.to_frozen()
    assert isinstance(sig2.minhash, FrozenMinHash)


def test_load_signature_ksize_nonint(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    sig = SourmashSignature(e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s, ksize='20'))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_empty(track_abundance):
    # edge case, but: empty minhash? :)
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)

    sig = SourmashSignature(e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 0
    assert sig2.similarity(sig) == 0


def test_roundtrip_scaled(track_abundance):
    e = MinHash(n=0, ksize=20, track_abundance=track_abundance,
                         max_hash=10)
    e.add_hash(5)
    sig = SourmashSignature(e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert e.scaled == e2.scaled

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_seed(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance,
                         seed=10)
    e.add_hash(5)
    sig = SourmashSignature(e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert e.seed == e2.seed

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_similarity_downsample(track_abundance):
    e = MinHash(n=0, ksize=20, track_abundance=track_abundance,
                         max_hash=2**63)
    f = MinHash(n=0, ksize=20, track_abundance=track_abundance,
                         max_hash=2**2)

    e.add_hash(1)
    e.add_hash(5)
    assert len(e.hashes) == 2

    f.add_hash(1)
    f.add_hash(5)                 # should be discarded due to max_hash
    assert len(f.hashes) == 1

    ee = SourmashSignature(e)
    ff = SourmashSignature(f)

    with pytest.raises(ValueError) as e:       # mismatch in max_hash
        ee.similarity(ff)

    assert 'mismatch in scaled; comparison fail' in str(e.value)

    x = ee.similarity(ff, downsample=True)
    assert round(x, 1) == 1.0


def test_add_sequence_bad_dna(track_abundance):
    # test add_sequence behavior on bad DNA
    mh = MinHash(n=1, ksize=21)
    sig = SourmashSignature(mh)

    with pytest.raises(ValueError) as e:
        sig.add_sequence("N" * 21, force=False)

    assert 'invalid DNA character in input k-mer: NNNNNNNNNNNNNNNNNNNNN' in str(e.value)


def test_md5(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_hash(5)
    sig = SourmashSignature(e)
    assert sig.md5sum() == 'eae27d77ca20db309e056e3d2dcd7d69', sig.md5sum()


def test_str_1(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature(e, name='foo')
    assert str(sig) == 'foo'


def test_str_2(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature(e, filename='foo.txt')
    assert str(sig) == 'foo.txt'


def test_str_3(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature(e, name='foo',
                            filename='foo.txt')
    assert str(sig) == 'foo'


def test_name_4(track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature(e)
    assert str(sig) == sig.md5sum()[:8]


def test_save_load_multisig(track_abundance):
    e1 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature(e1)

    e2 = MinHash(n=1, ksize=25, track_abundance=track_abundance)
    sig2 = SourmashSignature(e2)

    x = save_signatures([sig1, sig2])
    y = list(load_signatures(x))

    print(x)

    assert len(y) == 2
    assert sig1 in y                      # order not guaranteed, note.
    assert sig2 in y
    assert sig1 != sig2


def test_load_one_fail_nosig(track_abundance):
    x = save_signatures([])
    print((x,))
    with pytest.raises(ValueError):
        y = load_one_signature(x)


def test_load_one_succeed(track_abundance):
    e1 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature(e1)

    x = save_signatures([sig1])

    y = load_one_signature(x)
    assert sig1 == y


def test_load_one_fail_multisig(track_abundance):
    e1 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature(e1)

    e2 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig2 = SourmashSignature(e2)

    x = save_signatures([sig1, sig2])

    with pytest.raises(ValueError):
        y = load_one_signature(x)


def test_save_minified(track_abundance):
    e1 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature(e1, name="foo")

    e2 = MinHash(n=1, ksize=25, track_abundance=track_abundance)
    sig2 = SourmashSignature(e2, name="bar baz")

    x = save_signatures([sig1, sig2])
    assert b'\n' not in x
    assert len(x.split(b'\n')) == 1

    y = list(load_signatures(x))
    assert len(y) == 2
    assert any(sig.name == 'foo' for sig in y)
    assert any(sig.name == 'bar baz' for sig in y)


def test_load_minified(track_abundance):
    sigfile = utils.get_test_data('genome-s10+s11.sig')
    sigs = load_signatures(sigfile)

    minified = save_signatures(sigs)
    with open(sigfile, 'r') as f:
        orig_file = f.read()
    assert len(minified) < len(orig_file)
    assert b'\n' not in minified


def test_load_compressed(track_abundance):
    e1 = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature(e1)

    x = save_signatures([sig1], compression=5)

    y = load_one_signature(x)
    assert sig1 == y

    sigfile = utils.get_test_data('genome-s10+s11.sig.gz')
    sigs = load_signatures(sigfile)


def test_binary_fp(tmpdir, track_abundance):
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)

    path = tmpdir.join("1.sig")
    with open(str(path), 'wb') as fp:
        sig = SourmashSignature(e)
        s = save_signatures([sig], fp)


def test_load_signatures_no_file_do_raise(tmpdir):
    path = tmpdir.join("dne.sig")
    siglist = load_signatures(path, do_raise=True)
    with pytest.raises(Exception):
        list(siglist)


def test_load_signatures_no_file_do_not_raise(tmpdir):
    path = tmpdir.join("dne.sig")
    siglist = load_signatures(path)
    siglist = list(siglist)
    assert not siglist


def test_max_containment():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 5))

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    assert ss1.contained_by(ss2) == 1/4
    assert ss2.contained_by(ss1) == 1/2
    assert ss1.max_containment(ss2) == 1/2
    assert ss2.max_containment(ss1) == 1/2


def test_max_containment_empty():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    assert ss1.contained_by(ss2) == 0
    assert ss2.contained_by(ss1) == 0
    assert ss1.max_containment(ss2) == 0
    assert ss2.max_containment(ss1) == 0


def test_max_containment_equal():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=False)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=False)

    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 3, 4))

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    assert ss1.contained_by(ss2) == 1
    assert ss2.contained_by(ss1) == 1
    assert ss1.max_containment(ss2) == 1
    assert ss2.max_containment(ss1) == 1


def test_containment_ANI():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2, ksize=31)

    s1_cont_s2 = ss1.containment_ani(ss2, estimate_ci =True)
    s2_cont_s1 = ss2.containment_ani(ss1, estimate_ci =True)
    print("\nss1 contained by ss2", s1_cont_s2)
    print("ss2 contained by ss1", s2_cont_s1)

    assert (round(s1_cont_s2.ani,3), s1_cont_s2.ani_low, s1_cont_s2.ani_high) == (1.0,1.0,1.0)
    assert (round(s2_cont_s1.ani,3), round(s2_cont_s1.ani_low,3), round(s2_cont_s1.ani_high,3)) == (0.966, 0.965, 0.967)

    s1_mc_s2 = ss1.max_containment_ani(ss2, estimate_ci =True)
    s2_mc_s1 = ss2.max_containment_ani(ss1, estimate_ci =True)
    print("mh1 max containment", s1_mc_s2)
    print("mh2 max containment", s2_mc_s1)
    s1_mc_s2.size_is_inaccurate = False
    s2_mc_s1.size_is_inaccurate = False
    assert s1_mc_s2 == s2_mc_s1
    assert (round(s1_mc_s2.ani, 3), round(s1_mc_s2.ani_low, 3), round(s1_mc_s2.ani_high, 3)) == (1.0,1.0,1.0)


def test_containment_ANI_precalc_containment():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2, ksize=31)
    # precalc containments and assert same results
    s1c = ss1.contained_by(ss2)
    s2c = ss2.contained_by(ss1)
    mc = max(s1c, s2c)

    assert ss1.containment_ani(ss2, estimate_ci=True) ==  ss1.containment_ani(ss2, containment=s1c, estimate_ci=True)
    assert ss2.containment_ani(ss1) ==  ss2.containment_ani(ss1, containment=s2c)
    assert ss1.max_containment_ani(ss2) ==  ss2.max_containment_ani(ss1)
    assert ss1.max_containment_ani(ss2) ==  ss1.max_containment_ani(ss2, max_containment=mc)
    assert ss1.max_containment_ani(ss2) ==  ss2.max_containment_ani(ss1, max_containment=mc)


def test_avg_containment():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2, ksize=31)
    # check average_containment_ani
    ac_s1 = ss1.avg_containment(ss2)
    ac_s2 = ss2.avg_containment(ss1)
    assert ac_s1 == ac_s2 == (ss1.contained_by(ss2) + ss2.contained_by(ss1))/2 == 0.6619979467456603


def test_avg_containment_ani():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2, ksize=31)
    # check average_containment_ani
    ac_s1 = ss1.avg_containment_ani(ss2)
    ac_s2 = ss2.avg_containment_ani(ss1)
    assert ac_s1 == ac_s2 == (ss1.containment_ani(ss2).ani + ss2.containment_ani(ss1).ani)/2 


def test_containment_ANI_downsample():
    f2 = utils.get_test_data('2+63.fa.sig')
    f3 = utils.get_test_data('47+63.fa.sig')
    ss2 = sourmash.load_one_signature(f2, ksize=31)
    ss3 = sourmash.load_one_signature(f3, ksize=31)
    # check that downsampling works properly
    print(ss2.minhash.scaled)

    ss2 = ss2.to_mutable()
    ss2.minhash = ss2.minhash.downsample(scaled=2000)
    assert ss2.minhash.scaled != ss3.minhash.scaled
    ds_s3c = ss2.containment_ani(ss3, downsample=True)
    ds_s4c = ss3.containment_ani(ss2, downsample=True)
    mc_w_ds_1 =  ss2.max_containment_ani(ss3, downsample=True)
    mc_w_ds_2 =  ss3.max_containment_ani(ss2, downsample=True)

    with pytest.raises(ValueError) as e:
        ss2.containment_ani(ss3)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    with pytest.raises(ValueError) as e:
        ss2.max_containment_ani(ss3)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    ss3 = ss3.to_mutable()
    ss3.minhash = ss3.minhash.downsample(scaled=2000)
    assert ss2.minhash.scaled == ss3.minhash.scaled
    ds_s3c_manual = ss2.containment_ani(ss3)
    ds_s4c_manual = ss3.containment_ani(ss2)
    ds_mc_manual =  ss2.max_containment_ani(ss3)
    assert ds_s3c == ds_s3c_manual
    assert ds_s4c == ds_s4c_manual
    assert mc_w_ds_1 == mc_w_ds_2 == ds_mc_manual


def test_jaccard_ANI():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2)

    print("\nJACCARD_ANI", ss1.jaccard_ani(ss2))

    s1_jani_s2 = ss1.jaccard_ani(ss2)
    s2_jani_s1 = ss2.jaccard_ani(ss1)

    assert s1_jani_s2 == s2_jani_s1
    assert (s1_jani_s2.ani, s1_jani_s2.p_nothing_in_common, s1_jani_s2.jaccard_error) == (0.9783711630110239, 0.0, 3.891666770716877e-07)


def test_jaccard_ANI_untrustworthy():
    f1 = utils.get_test_data('2.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2)

    print("\nJACCARD_ANI", ss1.jaccard_ani(ss2))

    s1_jani_s2 = ss1.jaccard_ani(ss2, err_threshold=1e-7)

    # since size is inaccurate on 2.fa.sig, need to override to be able to get ani
    s1_jani_s2.size_is_inaccurate = False

    assert s1_jani_s2.ani == None
    assert s1_jani_s2.je_exceeds_threshold==True
    assert s1_jani_s2.je_threshold == 1e-7


def test_jaccard_ANI_precalc_jaccard():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2)
    # precalc jaccard and assert same result
    jaccard = ss1.jaccard(ss2)
    print("\nJACCARD_ANI", ss1.jaccard_ani(ss2,jaccard=jaccard))

    assert ss1.jaccard_ani(ss2) == ss1.jaccard_ani(ss2, jaccard=jaccard) == ss2.jaccard_ani(ss1, jaccard=jaccard)
    wrong_jaccard = jaccard - 0.1
    assert ss1.jaccard_ani(ss2) != ss1.jaccard_ani(ss2, jaccard=wrong_jaccard)


def test_jaccard_ANI_downsample():
    f1 = utils.get_test_data('47+63.fa.sig')
    f2 = utils.get_test_data('2+63.fa.sig')
    ss1 = sourmash.load_one_signature(f1, ksize=31)
    ss2 = sourmash.load_one_signature(f2)

    print(ss1.minhash.scaled)
    ss1 = ss1.to_mutable()
    ss1.minhash = ss1.minhash.downsample(scaled=2000)
    assert ss1.minhash.scaled != ss2.minhash.scaled
    with pytest.raises(ValueError) as e:
        ss1.jaccard_ani(ss2)
        assert "ValueError: mismatch in scaled; comparison fail" in e

    ds_s1c = ss1.jaccard_ani(ss2, downsample=True)
    ds_s2c = ss2.jaccard_ani(ss1, downsample=True)

    ss2 = ss2.to_mutable()
    ss2.minhash = ss2.minhash.downsample(scaled=2000)
    assert ss1.minhash.scaled == ss2.minhash.scaled
    ds_j_manual = ss1.jaccard_ani(ss2)
    assert ds_s1c == ds_s2c == ds_j_manual


def test_frozen_signature_update_1(track_abundance):
    # setting .name should fail on a FrozenSourmashSignature
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    ss = SourmashSignature(e, name='foo').to_frozen()

    with pytest.raises(ValueError):
        ss.name = 'foo2'


def test_frozen_signature_update_2(track_abundance):
    # setting .minhash should fail on a FrozenSourmashSignature
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    e2 = e.copy_and_clear()
    ss = SourmashSignature(e, name='foo').to_frozen()

    with pytest.raises(ValueError):
        ss.minhash = e2


def test_frozen_signature_update_3(track_abundance):
    # setting .minhash should succeed with update() context manager
    e = MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_kmer("AT" * 10)
    ss = SourmashSignature(e, name='foo').to_frozen()

    with ss.update() as ss2:
        ss2.name = 'foo2'

    assert ss2.name == 'foo2'
