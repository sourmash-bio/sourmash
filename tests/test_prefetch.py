"""
Tests for `sourmash prefetch` command-line and API functionality.
"""
import os
import csv
import pytest

import sourmash_tst_utils as utils
import sourmash


@utils.in_tempdir
def test_prefetch_basic(c):
    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "all sketches will be downsampled to scaled=1000" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


@utils.in_tempdir
def test_prefetch_csv_out(c):
    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    expected_intersect_bp = [2529000, 5177000]
    with open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_intersect_bp):
            assert int(row['intersect_bp']) == expected


@utils.in_tempdir
def test_prefetch_matches(c):
    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    matches_out = c.output('matches.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)

    sigs = sourmash.load_file_as_index(matches_out)

    expected_matches = [sig63, sig47]
    for (match, expected) in zip(sigs.signatures(), expected_matches):
        ss = sourmash.load_one_signature(expected, ksize=31)
        assert match == ss


@utils.in_tempdir
def test_prefetch_no_num_query(c):
    # can't do prefetch with num signatures for query
    sig47 = utils.get_test_data('num/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig47)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0


@utils.in_tempdir
def test_prefetch_no_num_subj(c):
    # can't do prefetch with num signatures for query; no matches!
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('num/63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "ERROR in prefetch_databases:" in c.last_result.err
    assert "no signatures to search" in c.last_result.err
