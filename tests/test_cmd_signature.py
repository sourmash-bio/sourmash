"""
Tests for the 'sourmash signature' command line.
"""
from __future__ import print_function, unicode_literals

from . import sourmash_tst_utils as utils
import sourmash

## command line tests


def test_run_sourmash_signature_cmd():
    status, out, err = utils.runscript('sourmash', ['signature'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_sig_cmd():
    status, out, err = utils.runscript('sourmash', ['sig'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


@utils.in_tempdir
def test_sig_merge_1(c):
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('sig', 'merge', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_intersect_1(c):
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63-intersect.fa.sig')
    c.run_sourmash('sig', 'intersect', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_intersect_sig = sourmash.load_one_signature(sig47and63)
    actual_intersect_sig = sourmash.load_one_signature(out)

    print(test_intersect_sig.minhash)
    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_intersect_sig.minhash
