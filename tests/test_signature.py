import os

import pytest

import sourmash
from sourmash.signature import SourmashSignature, save_signatures, \
    load_signatures, load_one_signature
import sourmash_tst_utils as utils
from sourmash import MinHash
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
