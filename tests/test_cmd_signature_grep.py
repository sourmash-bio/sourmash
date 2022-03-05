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


def test_grep_1_sig_name_exclude(runtmp):
    # search on substring in name, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    # no matches!
    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', '-v', 'Shewanella', sig47)


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


def test_grep_4_no_manifest(runtmp):
    # fail search when no manifest, by default
    sbt = utils.get_test_data('v6.sbt.zip')

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('sig', 'grep', 'e60265', sbt)

    print(runtmp.last_result.err)
    assert 'ERROR on filename' in runtmp.last_result.err
    assert 'sig grep requires a manifest by default, but no manifest present.' in runtmp.last_result.err


def test_grep_4_no_manifest_ok(runtmp):
    # generate manifest if --no-require-manifest
    sbt = utils.get_test_data('v6.sbt.zip')

    runtmp.run_sourmash('sig', 'grep', 'e60265', sbt, '--no-require-manifest')

    ss = sourmash.load_signatures(runtmp.last_result.out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'e60265' in ss.md5sum()


def test_grep_5_zip_include(runtmp):
    # search zip, include on case sensitive match to name
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_include_case_insensitive(runtmp):
    # search zip, include on case insensitive match to name
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-i', 'os223', allzip)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_exclude(runtmp):
    # search zip, exclude on case-sensitive match
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-v', 'OS185', allzip)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_exclude_case_insensitive(runtmp):
    # search zip, exclude on case-insensitive match
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-vi', 'os185', allzip)

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_6_zip_manifest_csv(runtmp):
    # do --csv and use result as picklist
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip,
                        '--csv', 'match.csv')

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'

    # now run cat with picklist
    runtmp.run_sourmash('sig', 'cat', allzip,
                        '--picklist', 'match.csv::manifest')

    out = runtmp.last_result.out
    ss = sourmash.load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'
