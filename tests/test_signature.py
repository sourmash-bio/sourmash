from __future__ import print_function, unicode_literals

import pytest

import sourmash_lib
from sourmash_lib.signature import SourmashSignature, save_signatures, \
    load_signatures, load_one_signature


def test_hashable(track_abundance):
    # check: can we use signatures as keys in dictionaries and sets?
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add("AT" * 10)

    sig = SourmashSignature('', e)

    x = set()
    x.add(sig)


def test_str(track_abundance):
    # signatures should be printable
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add("AT" * 10)

    sig = SourmashSignature('', e)

    print(sig)
    assert str(sig) == 'SourmashSignature(59502a74)'
    assert repr(sig) == 'SourmashSignature(59502a74)'

    sig.d['name'] = 'fizbar'
    assert str(sig) == 'SourmashSignature(\'fizbar\', 59502a74)'
    assert repr(sig) == 'SourmashSignature(\'fizbar\', 59502a74)'


def test_roundtrip(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add("AT" * 10)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_load_signature_ksize_nonint(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add("AT" * 10)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s, ksize='20'))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_empty(track_abundance):
    # edge case, but: empty minhash? :)
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)

    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 0
    assert sig2.similarity(sig) == 0


def test_roundtrip_max_hash(track_abundance):
    e = sourmash_lib.MinHash(n=0, ksize=20, track_abundance=track_abundance,
                             max_hash=10)
    e.add_hash(5)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert e.max_hash == e2.max_hash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_seed(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance,
                             seed=10)
    e.add_hash(5)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert e.seed == e2.seed

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_empty_email(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add("AT" * 10)
    sig = SourmashSignature('', e)
    s = save_signatures([sig])
    print(s)
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.minhash

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_similarity_downsample(track_abundance):
    e = sourmash_lib.MinHash(n=0, ksize=20, track_abundance=track_abundance,
                             max_hash=2**63)
    f = sourmash_lib.MinHash(n=0, ksize=20, track_abundance=track_abundance,
                             max_hash=2**2)

    e.add_hash(1)
    e.add_hash(5)
    assert len(e.get_mins()) == 2

    f.add_hash(1)
    f.add_hash(5)                 # should be discarded due to max_hash
    assert len(f.get_mins()) == 1

    ee = SourmashSignature('', e)
    ff = SourmashSignature('', f)

    with pytest.raises(ValueError):       # mismatch in max_hash
        ee.similarity(ff)

    x = ee.similarity(ff, downsample=True)
    assert round(x, 1) == 1.0


def test_md5(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    e.add_hash(5)
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig._save())
    assert sig.md5sum() == 'eae27d77ca20db309e056e3d2dcd7d69', sig.md5sum()


def test_name(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature('titus@idyll.org', e, name='foo')
    assert sig.name() == 'foo'


def test_name_2(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature('titus@idyll.org', e, filename='foo.txt')
    assert sig.name() == 'foo.txt'


def test_name_3(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature('titus@idyll.org', e, name='foo',
                            filename='foo.txt')
    assert sig.name() == 'foo'


def test_name_4(track_abundance):
    e = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig = SourmashSignature('titus@idyll.org', e)
    assert sig.name() == sig.md5sum()[:8]


def test_save_load_multisig(track_abundance):
    e1 = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature('titus@idyll.org', e1)

    e2 = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig2 = SourmashSignature('titus2@idyll.org', e2)

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
    e1 = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature('titus@idyll.org', e1)

    x = save_signatures([sig1])

    y = load_one_signature(x)
    assert sig1 == y


def test_load_one_fail_multisig(track_abundance):
    e1 = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig1 = SourmashSignature('titus@idyll.org', e1)

    e2 = sourmash_lib.MinHash(n=1, ksize=20, track_abundance=track_abundance)
    sig2 = SourmashSignature('titus2@idyll.org', e2)

    x = save_signatures([sig1, sig2])

    with pytest.raises(ValueError):
        y = load_one_signature(x)
