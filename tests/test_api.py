from __future__ import print_function, unicode_literals
import sourmash


def test_sourmash_signature_api():
    e = sourmash.MinHash(n=1, ksize=20)
    sig = sourmash.SourmashSignature(e)

    s = sourmash.save_signatures([sig])
    sig_x1 = sourmash.load_one_signature(s)
    sig_x2 = list(sourmash.load_signatures(s))[0]

    assert sig_x1 == sig
    assert sig_x2 == sig
