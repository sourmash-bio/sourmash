"""
Tests for the 'sourmash signature' command line.
"""
from __future__ import print_function, unicode_literals
import csv
import shutil
import os

import pytest

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
    # merge of 47 & 63 should be union of mins
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
def test_sig_merge_1_multisig(c):
    # merge of 47 & 63 should be union of mins; here, sigs are in same file.
    multisig = utils.get_test_data('47+63-multisig.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('sig', 'merge', multisig, '--flatten')

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_1_ksize_moltype(c):
    # check ksize, moltype args
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('sig', 'merge', sig47, sig63, '--dna', '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_2(c):
    # merge of 47 with nothing should be 47
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'merge', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_3_abund_ab_ok(c):
    # merge of 47 and 63 with abund should work
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    c.run_sourmash('sig', 'merge', sig47abund, sig63abund)
    actual_merge_sig = sourmash.load_one_signature(c.last_result.out)
    # @CTB: should check that this merge did what we think it should do!


@utils.in_tempdir
def test_sig_merge_3_abund_ab(c):
    # merge of 47 with abund, with 63 without, should fail; and vice versa
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError) as e:
        c.run_sourmash('sig', 'merge', sig47, sig63abund)

    print(c.last_result)
    assert 'incompatible signatures: track_abundance is False in first sig, True in second' in c.last_result.err


@utils.in_tempdir
def test_sig_merge_3_abund_ba(c):
    # merge of 47 without abund, with 63 with, should fail
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError) as e:
        c.run_sourmash('sig', 'merge', sig63abund, sig47)

    print(c.last_result)
    assert 'incompatible signatures: track_abundance is True in first sig, False in second' in c.last_result.err


@utils.in_tempdir
def test_sig_filter_1(c):
    # test basic filtering
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')
    c.run_sourmash('sig', 'filter', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    filtered_sigs = list(sourmash.load_signatures(out))
    filtered_sigs.sort(key=lambda x: x.name())

    assert len(filtered_sigs) == 2

    mh47 = sourmash.load_one_signature(sig47).minhash
    mh63 = sourmash.load_one_signature(sig63).minhash

    assert filtered_sigs[0].minhash == mh47
    assert filtered_sigs[1].minhash == mh63


@utils.in_tempdir
def test_sig_filter_2(c):
    # test basic filtering
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    c.run_sourmash('sig', 'filter', '-m', '2', '-M', '5', sig47)

    # stdout should be new signature
    out = c.last_result.out

    filtered_sig = sourmash.load_one_signature(out)
    test_sig = sourmash.load_one_signature(sig47)

    abunds = test_sig.minhash.get_mins(True)
    abunds = { k: v for (k, v) in abunds.items() if v >= 2 and v <= 5 }
    assert abunds

    assert filtered_sig.minhash.get_mins(True) == abunds


@utils.in_tempdir
def test_sig_filter_3(c):
    # test basic filtering
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    c.run_sourmash('sig', 'filter', '-m', '2', sig47)

    # stdout should be new signature
    out = c.last_result.out

    filtered_sig = sourmash.load_one_signature(out)
    test_sig = sourmash.load_one_signature(sig47)

    abunds = test_sig.minhash.get_mins(True)
    abunds = { k: v for (k, v) in abunds.items() if v >= 2 }
    assert abunds

    assert filtered_sig.minhash.get_mins(True) == abunds


@utils.in_tempdir
def test_sig_merge_flatten(c):
    # merge of 47 without abund, with 63 with, will succeed with --flatten
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')

    c.run_sourmash('sig', 'merge', sig63abund, sig47, '--flatten')

    print(c.last_result)
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_flatten_2(c):
    # merge of 47 with abund, with 63 with, will succeed with --flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')

    c.run_sourmash('sig', 'merge', sig63, sig47abund, '--flatten')

    print(c.last_result)
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_intersect_1(c):
    # intersect of 47 and 63 should be intersection of mins
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


@utils.in_tempdir
def test_sig_intersect_2(c):
    # intersect of 47 with abund and 63 with abund should be same
    # as without abund, i.e. intersect 'flattens'
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')
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


@utils.in_tempdir
def test_sig_intersect_3(c):
    # use --abundances-from to preserve abundances from sig #47
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')
    c.run_sourmash('sig', 'intersect', '--abundances-from', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    actual_intersect_sig = sourmash.load_one_signature(out)

    # actually do an intersection ourselves for the test
    mh47 = sourmash.load_one_signature(sig47).minhash
    mh63 = sourmash.load_one_signature(sig63).minhash
    mh47_abunds = mh47.get_mins(with_abundance=True)
    mh63_mins = set(mh63.get_mins())

    # get the set of mins that are in common
    mh63_mins.intersection_update(mh47_abunds)

    # take abundances from mh47 & create new sig
    mh47_abunds = { k: mh47_abunds[k] for k in mh63_mins }
    test_mh = mh47.copy_and_clear()
    test_mh.set_abundances(mh47_abunds)

    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_intersect_4(c):
    # use --abundances-from to preserve abundances from sig #47
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')
    c.run_sourmash('sig', 'intersect', '--abundances-from', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    actual_intersect_sig = sourmash.load_one_signature(out)

    # actually do an intersection ourselves for the test
    mh47 = sourmash.load_one_signature(sig47).minhash
    mh63 = sourmash.load_one_signature(sig63).minhash
    mh47_abunds = mh47.get_mins(with_abundance=True)
    mh63_mins = set(mh63.get_mins())

    # get the set of mins that are in common
    mh63_mins.intersection_update(mh47_abunds)

    # take abundances from mh47 & create new sig
    mh47_abunds = { k: mh47_abunds[k] for k in mh63_mins }
    test_mh = mh47.copy_and_clear()
    test_mh.set_abundances(mh47_abunds)

    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_intersect_5(c):
    # use --abundances-from to preserve abundances from sig #47
    # make sure that you can't specify a flat sig for --abundances-from
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'intersect', '--abundances-from', sig47, sig63)


@utils.in_tempdir
def test_sig_subtract_1(c):
    # subtract of 63 from 47
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'subtract', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test1_sig = sourmash.load_one_signature(sig47)
    test2_sig = sourmash.load_one_signature(sig63)
    actual_subtract_sig = sourmash.load_one_signature(out)

    mins = set(test1_sig.minhash.get_mins())
    mins -= set(test2_sig.minhash.get_mins())

    assert set(actual_subtract_sig.minhash.get_mins()) == set(mins)


@utils.in_tempdir
def test_sig_subtract_1_multisig(c):
    # subtract of everything from 47
    sig47 = utils.get_test_data('47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'subtract', sig47, multisig, '--flatten')

    # stdout should be new signature
    out = c.last_result.out

    actual_subtract_sig = sourmash.load_one_signature(out)

    assert not set(actual_subtract_sig.minhash.get_mins())


@utils.in_tempdir
def test_sig_subtract_2(c):
    # subtract of 63 from 47 should fail if 47 has abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_subtract_3(c):
    # subtract of 63 from 47 should fail if 63 has abund
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_intersect_2(c):
    # intersect of 47 and nothing should be self
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'intersect', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_intersect_sig = sourmash.load_one_signature(sig47)
    actual_intersect_sig = sourmash.load_one_signature(out)

    print(test_intersect_sig.minhash)
    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_intersect_sig.minhash


@utils.in_tempdir
def test_sig_intersect_2_multisig(c):
    # intersect of all the multisig stuff should be nothing
    sig47 = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'intersect', sig47)

    # stdout should be new signature
    out = c.last_result.out

    actual_intersect_sig = sourmash.load_one_signature(out)

    assert not len(actual_intersect_sig.minhash)


@utils.in_tempdir
def test_sig_rename_1(c):
    # set new name for 47
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'rename', sig47, 'fiz bar')

    # stdout should be new signature
    out = c.last_result.out

    test_rename_sig = sourmash.load_one_signature(sig47)
    actual_rename_sig = sourmash.load_one_signature(out)

    print(test_rename_sig.minhash)
    print(actual_rename_sig.minhash)

    assert actual_rename_sig.minhash == test_rename_sig.minhash
    assert test_rename_sig.name() != actual_rename_sig.name()
    assert actual_rename_sig.name() == 'fiz bar'


@utils.in_tempdir
def test_sig_rename_1_multisig(c):
    # set new name for multiple signatures/files
    multisig = utils.get_test_data('47+63-multisig.sig')
    other_sig = utils.get_test_data('2.fa.sig')
    c.run_sourmash('sig', 'rename', multisig, other_sig, 'fiz bar')

    # stdout should be new signature
    out = c.last_result.out

    n = 0
    for sig in sourmash.load_signatures(out):
        assert sig.name() == 'fiz bar'
        n += 1

    assert n == 9, n


@utils.in_tempdir
def test_sig_rename_2_output_to_same(c):
    # change name of signature "in place", same output file
    sig47 = utils.get_test_data('47.fa.sig')
    inplace = c.output('inplace.sig')
    shutil.copyfile(sig47, inplace)

    print(inplace)

    c.run_sourmash('sig', 'rename', '-d', inplace, 'fiz bar', '-o', inplace)

    actual_rename_sig = sourmash.load_one_signature(inplace)
    assert actual_rename_sig.name() == 'fiz bar'


@utils.in_tempdir
def test_sig_extract_1(c):
    # extract 47 from 47... :)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'extract', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_2(c):
    # extract matches to 47's md5sum from among several
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--md5', '09a0869')

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    print(test_extract_sig.minhash)
    print(actual_extract_sig.minhash)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_3(c):
    # extract nothing (no md5 match)
    sig47 = utils.get_test_data('47.fa.sig')
    with pytest.raises(ValueError) as exc:
        c.run_sourmash('sig', 'extract', sig47, '--md5', 'FOO')


@utils.in_tempdir
def test_sig_extract_4(c):
    # extract matches to 47's name from among several signatures
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--name', 'NC_009665.1')

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    print(test_extract_sig.minhash)
    print(actual_extract_sig.minhash)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_5(c):
    # extract nothing (no name match)
    sig47 = utils.get_test_data('47.fa.sig')
    with pytest.raises(ValueError) as exc:
        c.run_sourmash('sig', 'extract', sig47, '--name', 'FOO')


@utils.in_tempdir
def test_sig_extract_6(c):
    # extract matches to several names from among several signatures
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = sourmash.load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 2


@utils.in_tempdir
def test_sig_flatten_1(c):
    # extract matches to several names from among several signatures & flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'flatten', sig47abund, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = sourmash.load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


@utils.in_tempdir
def test_sig_downsample_1_scaled(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sig47)
    actual_downsample_sig = sourmash.load_one_signature(out)

    test_mh = test_downsample_sig.minhash.downsample_scaled(10000)

    assert actual_downsample_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_downsample_1_scaled_downsample_multisig(c):
    # downsample many scaled signatures in one file
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000', multisig)

    # stdout should be new signatures
    out = c.last_result.out

    for sig in sourmash.load_signatures(out):
        assert sig.minhash.scaled == 10000


@utils.in_tempdir
def test_sig_downsample_1_scaled_to_num(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'downsample', '--num', '500', sig47)

    # stdout should be new signature
    out = c.last_result.out

    actual_downsample_sig = sourmash.load_one_signature(out)
    actual_mins = actual_downsample_sig.minhash.get_mins()
    actual_mins = list(actual_mins)
    actual_mins.sort()

    test_downsample_sig = sourmash.load_one_signature(sig47)
    test_mins = test_downsample_sig.minhash.get_mins()
    test_mins = list(test_mins)
    test_mins.sort()
    test_mins = test_mins[:500]           # take 500 smallest

    assert actual_mins == test_mins


@utils.in_tempdir
def test_sig_downsample_1_scaled_to_num_fail(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', '--num', '50000', sig47)


@utils.in_tempdir
def test_sig_downsample_1_scaled_empty(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', sig47)


@utils.in_tempdir
def test_sig_downsample_2_num(c):
    # downsample a num signature
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')
    c.run_sourmash('sig', 'downsample', '--num', '500',
                   '-k', '21', '--dna', sigs11)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sigs11, ksize=21,
                                                      select_moltype='DNA')
    actual_downsample_sig = sourmash.load_one_signature(out)
    test_mh = test_downsample_sig.minhash.downsample_n(500)

    assert actual_downsample_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_downsample_2_num_to_scaled(c):
    # downsample a num signature and convert it into a scaled sig
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000',
                   '-k', '21', '--dna', sigs11)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sigs11, ksize=21,
                                                      select_moltype='DNA')
    actual_downsample_sig = sourmash.load_one_signature(out)

    test_mins = test_downsample_sig.minhash.get_mins()
    actual_mins = actual_downsample_sig.minhash.get_mins()

    # select those mins that are beneath the new max hash...
    max_hash = actual_downsample_sig.minhash.max_hash
    test_mins_down = { k for k in test_mins if k < max_hash }
    assert test_mins_down == set(actual_mins)


@utils.in_tempdir
def test_sig_downsample_2_num_to_scaled_fail(c):
    # downsample a num signature and FAIL to convert it into a scaled sig
    # because new scaled is too low
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', '--scaled', '100',
                       '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_downsample_2_num_empty(c):
    # downsample a num signature
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_describe_1(c):
    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'describe', sig47)

    out = c.last_result.out
    print(c.last_result)

    expected_output = """\
signature: NC_009665.1 Shewanella baltica OS185, complete genome
source file: 47.fa
md5: 09a08691ce52952152f0e866a59f6261
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 5177
signature license: CC0
""".splitlines()
    for line in expected_output:
        assert line.strip() in out

@utils.in_tempdir
def test_sig_describe_1_hp(c):
    # get basic info on a signature
    testdata = utils.get_test_data('short.fa')
    c.run_sourmash('compute', '-k', '21,30',
                   '--dayhoff', '--hp', '--protein',
                   '--dna',
                   testdata)
    # stdout should be new signature
    computed_sig = os.path.join(c.location, 'short.fa.sig')
    c.run_sourmash('sig', 'describe', computed_sig)

    out = c.last_result.out
    print(c.last_result)

    # Add final trailing slash for this OS
    testdata_dirname = os.path.dirname(testdata) + os.sep
    location = c.location + os.sep

    expected_output = """\
---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: e45a080101751e044d6df861d3d0f3fd
k=21 molecule=protein num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: ef4fa1f3a90f3873187370f1eacc0d9a
k=21 molecule=dayhoff num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0
---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: 20be00d9d577da9faeb77477bf07d3fb
k=21 molecule=hp num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0
---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: 1136a8a68420bd93683e45cdaf109b80
k=21 molecule=DNA num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: 4244d1612598af044e799587132f007e
k=30 molecule=protein num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: 5647819f2eac913e04af51c8d548ad56
k=30 molecule=dayhoff num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: ad1e329dd98b5e32422e9decf298aa5f
k=30 molecule=hp num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: short.fa
source file: short.fa
md5: 71f7c111c01785e5f38efad45b00a0e1
k=30 molecule=DNA num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0
""".splitlines()
    for line in out.splitlines():
        cleaned_line = line.strip().replace(
            testdata_dirname, '').replace(location, '')
        assert cleaned_line in expected_output


@utils.in_tempdir
def test_sig_describe_1_multisig(c):
    # get basic info on multiple signatures in a single file
    sigs = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'describe', sigs)

    out = c.last_result.out
    print(c.last_result)

    expected_output = """\
signature: NC_009665.1 Shewanella baltica OS185, complete genome
signature: NC_009661.1 Shewanella baltica OS185 plasmid pS18501, complete sequence
signature: NC_011663.1 Shewanella baltica OS223, complete genome
signature: NC_011664.1 Shewanella baltica OS223 plasmid pS22301, complete sequence
signature: NC_011668.1 Shewanella baltica OS223 plasmid pS22302, complete sequence
signature: NC_011665.1 Shewanella baltica OS223 plasmid pS22303, complete sequence""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@utils.in_tempdir
def test_sig_describe_2(c):
    # get info in CSV spreadsheet
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'describe', sig47, sig63, '--csv', 'out.csv')

    expected_md5 = ['09a08691ce52952152f0e866a59f6261',
                    '38729c6374925585db28916b82a6f513']

    with open(c.output('out.csv'), 'rt') as fp:
        r = csv.DictReader(fp)

        n = 0

        for row, md5 in zip(r, expected_md5):
            assert row['md5'] == md5
            n += 1

        assert n == 2


@utils.in_tempdir
def test_sig_overlap(c):
    # get overlap details
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'overlap', sig47, sig63)
    out = c.last_result.out

    print(out)

    # md5s
    assert '09a08691ce52952152f0e866a59f6261' in out
    assert '38729c6374925585db28916b82a6f513' in out

    assert 'similarity:                  0.32069' in out
    assert 'number of hashes in common:  2529' in out


@utils.in_tempdir
def test_import_export_1(c):
    # check to make sure we can import what we've exported!
    inp = utils.get_test_data('genome-s11.fa.gz.sig')
    outp = c.output('export.json')

    c.run_sourmash('sig', 'export', inp, '-o', outp, '-k', '21', '--dna')
    c.run_sourmash('sig', 'import', outp)

    original = sourmash.load_one_signature(inp, ksize=21, select_moltype='DNA')
    roundtrip = sourmash.load_one_signature(c.last_result.out)

    assert original.minhash == roundtrip.minhash


@utils.in_tempdir
def test_import_export_2(c):
    # check to make sure we can import a mash JSON dump file.
    # NOTE: msh.json_dump file calculated like so:
    #   mash sketch -s 500 -k 21 ./tests/test-data/genome-s11.fa.gz
    #   mash info -d ./tests/test-data/genome-s11.fa.gz.msh > tests/test-data/genome-s11.fa.gz.msh.json_dump
    #
    sig1 = utils.get_test_data('genome-s11.fa.gz.sig')
    msh_sig = utils.get_test_data('genome-s11.fa.gz.msh.json_dump')

    c.run_sourmash('sig', 'import', msh_sig)
    imported = sourmash.load_one_signature(c.last_result.out)
    compare = sourmash.load_one_signature(sig1, ksize=21, select_moltype='DNA')

    assert imported.minhash == compare.minhash
