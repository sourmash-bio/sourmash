"""
Tests for `sourmash prefetch` command-line and API functionality.
"""

import sourmash_tst_utils as utils


@utils.in_tempdir
def test_prefetch_basic(c):
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err
