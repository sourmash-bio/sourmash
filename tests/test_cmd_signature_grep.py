"""
Tests for the 'sourmash signature grep' command line.
"""
import shutil
import os

import pytest

import sourmash_tst_utils as utils
import sourmash
from sourmash_tst_utils import SourmashCommandFailed

## command line tests


def test_grep_1_sig_name(runtmp):
    # search on substring in name
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', 'Shewanella', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella' in ss.name
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_1_sig_name_case_sensitive(runtmp):
    # search on substring in name
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', 'shewanella', sig47)


def test_grep_1_sig_name_case_insensitive(runtmp):
    # search on substring in name, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '-i', 'shewanella', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella' in ss.name
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_2_sig_md5(runtmp):
    # search on substring in md5
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', 'ce52952152f0', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_2_sig_md5_case_sensitive(runtmp):
    # case sensitive no match
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', 'CE52952152f0', sig47)


def test_grep_2_sig_md5_case_insensitive(runtmp):
    # search on substring in md5, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '-i', 'CE52952152f0', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_3_filename(runtmp):
    # filename match
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '47.fa', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert '47.fa' in ss.filename
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_3_filename_regexp(runtmp):
    # search for a regexp on filename
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '^47.fa', sig47)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert '7.fa' in ss.filename
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'
