"""
Tests for the 'sourmash signature' command line.
"""
import csv
import shutil
import os
import glob

import pytest
import screed

import sourmash_tst_utils as utils
import sourmash
from sourmash.signature import load_signatures
from sourmash.manifest import CollectionManifest
from sourmash_tst_utils import SourmashCommandFailed

## command line tests


def _write_file(runtmp, basename, lines):
    loc = runtmp.output(basename)
    with open(loc, 'wt') as fp:
        fp.write("\n".join(lines))
    return loc


def test_run_sourmash_signature_cmd():
    status, out, err = utils.runscript('sourmash', ['signature'], fail_ok=True)
    assert not 'sourmash: error: argument cmd: invalid choice:' in err
    assert 'Manipulate signature files:' in out
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_sig_cmd():
    status, out, err = utils.runscript('sourmash', ['sig'], fail_ok=True)
    assert not 'sourmash: error: argument cmd: invalid choice:' in err
    assert 'Manipulate signature files:' in out
    assert status != 0                    # no args provided, ok ;)


def test_sig_merge_1_use_full_signature_in_cmd(runtmp):
    c = runtmp

    # merge of 47 & 63 should be union of mins
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('signature', 'merge', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


def test_sig_merge_1_fromfile_picklist(runtmp):
    c = runtmp

    # merge of 47 & 63 should be union of mins
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')

    from_file = _write_file(runtmp, 'list.txt', [sig47, sig63])
    picklist = _write_file(runtmp, 'pl.csv',
                           ['md5short', '09a08691', '38729c63'])

    c.run_sourmash('signature', 'merge', '--from-file', from_file,
                   '--picklist', f'{picklist}:md5short:md5short')

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


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
def test_sig_merge_1_name(c):
    # check name arg
    sig2 = utils.get_test_data('2.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    assignedSigName = 'SIG_NAME'
    outsig = c.output('merged2and63.sig')

    c.run_sourmash('sig', 'merge', sig2, sig63, '--dna', '-k', '31', '-o', "merged2and63.sig", '--name', assignedSigName )

    test_merge_sig = sourmash.load_one_signature(outsig)

    print("outsig", outsig)
    print("xx_test_merge_sig.name", test_merge_sig.name)

    assert assignedSigName == test_merge_sig.name


@utils.in_tempdir
def test_sig_merge_1_ksize_moltype(c):
    # check ksize, moltype args
    sig2 = utils.get_test_data('2.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig2and63 = utils.get_test_data('2+63.fa.sig')
    c.run_sourmash('sig', 'merge', sig2, sig63, '--dna', '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig2and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_1_ksize_moltype_fail(c):
    # check ksize, moltype args
    sig2 = utils.get_test_data('2.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig2and63 = utils.get_test_data('2+63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('sig', 'merge', sig2, sig63)

    assert "ERROR when merging signature" in str(exc.value)


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

    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('sig', 'merge', sig47, sig63abund)

    print(c.last_result)
    assert 'incompatible signatures: track_abundance is False in first sig, True in second' in c.last_result.err


@utils.in_tempdir
def test_sig_merge_3_abund_ba(c):
    # merge of 47 without abund, with 63 with, should fail
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as e:
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

    filtered_sigs = list(load_signatures(out))
    filtered_sigs.sort(key=lambda x: str(x))

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

    abunds = test_sig.minhash.hashes
    abunds = { k: v for (k, v) in abunds.items() if v >= 2 and v <= 5 }
    assert abunds

    assert filtered_sig.minhash.hashes == abunds


@utils.in_tempdir
def test_sig_filter_3(c):
    # test basic filtering
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    c.run_sourmash('sig', 'filter', '-m', '2', sig47)

    # stdout should be new signature
    out = c.last_result.out

    filtered_sig = sourmash.load_one_signature(out)
    test_sig = sourmash.load_one_signature(sig47)

    abunds = test_sig.minhash.hashes
    abunds = { k: v for (k, v) in abunds.items() if v >= 2 }
    assert abunds

    assert filtered_sig.minhash.hashes == abunds


@utils.in_tempdir
def test_sig_filter_3_ksize_select(c):
    # test filtering with ksize selectiong
    psw_mag = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    c.run_sourmash('sig', 'filter', '-m', '2', psw_mag, '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    filtered_sig = sourmash.load_one_signature(out)
    test_sig = sourmash.load_one_signature(psw_mag, ksize=31)

    abunds = test_sig.minhash.hashes
    abunds = { k: v for (k, v) in abunds.items() if v >= 2 }
    assert abunds

    assert filtered_sig.minhash.hashes == abunds


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


def test_sig_intersect_1(runtmp):
    c = runtmp

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


def test_sig_intersect_1_fromfile_picklist(runtmp):
    c = runtmp

    # intersect of 47 and 63 should be intersection of mins
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63-intersect.fa.sig')

    from_file = _write_file(runtmp, 'list.txt', [sig47, sig63])
    picklist = _write_file(runtmp, 'pl.csv',
                           ['md5short', '09a08691', '38729c63'])

    c.run_sourmash('signature', 'intersect', '--from-file', from_file,
                   '--picklist', f'{picklist}:md5short:md5short')

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
    mh47_abunds = mh47.hashes
    mh63_mins = set(mh63.hashes.keys())

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
    mh47_abunds = mh47.hashes
    mh63_mins = set(mh63.hashes.keys())

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

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'intersect', '--abundances-from', sig47, sig63)


@utils.in_tempdir
def test_sig_intersect_6_ksize_fail(c):
    # specify ksize to intersect 2.fa.sig with 47.fa.sig - 2.fa.sig contains
    # multiple ksizes.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'intersect', sig2, sig47)


@utils.in_tempdir
def test_sig_intersect_6_ksize_succeed(c):
    # specify ksize to intersect 2.fa.sig with 47.fa.sig - 2.fa.sig contains
    # multiple ksizes.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')

    c.run_sourmash('sig', 'intersect', '-k', '31', sig2, sig47)

    assert 'loaded and intersected 2 signatures' in c.last_result.err


@utils.in_tempdir
def test_sig_intersect_7(c):
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
def test_sig_intersect_8_multisig(c):
    # intersect of all the multisig stuff should be nothing
    sig47 = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'intersect', sig47)

    # stdout should be new signature
    out = c.last_result.out

    actual_intersect_sig = sourmash.load_one_signature(out)

    assert not len(actual_intersect_sig.minhash)


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

    mins = set(test1_sig.minhash.hashes.keys())
    mins -= set(test2_sig.minhash.hashes.keys())

    assert set(actual_subtract_sig.minhash.hashes.keys()) == set(mins)


@utils.in_tempdir
def test_sig_subtract_1_multisig(c):
    # subtract of everything from 47
    sig47 = utils.get_test_data('47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'subtract', sig47, multisig, '--flatten')

    # stdout should be new signature
    out = c.last_result.out

    actual_subtract_sig = sourmash.load_one_signature(out)

    assert not set(actual_subtract_sig.minhash.hashes.keys())


@utils.in_tempdir
def test_sig_subtract_2(c):
    # subtract of 63 from 47 should fail if 47 has abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_subtract_3(c):
    # subtract of 63 from 47 should fail if 63 has abund
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_subtract_4_ksize_fail(c):
    # subtract of 2 from 47 should fail without -k specified
    sig47 = utils.get_test_data('47.fa.sig')
    sig2 = utils.get_test_data('2.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'subtract', sig47, sig2)


@utils.in_tempdir
def test_sig_subtract_4_ksize_succeed(c):
    # subtract of 2 from 47 should fail without -k specified
    sig47 = utils.get_test_data('47.fa.sig')
    sig2 = utils.get_test_data('2.fa.sig')

    c.run_sourmash('sig', 'subtract', sig47, sig2, '-k', '31')
    assert 'loaded and subtracted 1 signatures' in c.last_result.err


def test_sig_rename_1(runtmp):
    c = runtmp

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
    assert test_rename_sig.name != actual_rename_sig.name
    assert actual_rename_sig.name == 'fiz bar'


def test_sig_rename_1_fromfile_picklist(runtmp):
    c = runtmp

    # set new name for 47
    sig47 = utils.get_test_data('47.fa.sig')

    from_file = _write_file(runtmp, 'list.txt', [sig47])
    picklist = _write_file(runtmp, 'pl.csv', ['md5short', '09a08691'])

    c.run_sourmash('sig', 'rename', '--from-file', from_file, 'fiz bar',
                   '--picklist', f'{picklist}:md5short:md5short')

    # stdout should be new signature
    out = c.last_result.out

    test_rename_sig = sourmash.load_one_signature(sig47)
    actual_rename_sig = sourmash.load_one_signature(out)

    print(test_rename_sig.minhash)
    print(actual_rename_sig.minhash)

    assert actual_rename_sig.minhash == test_rename_sig.minhash
    assert test_rename_sig.name != actual_rename_sig.name
    assert actual_rename_sig.name == 'fiz bar'


@utils.in_tempdir
def test_sig_rename_1_multisig(c):
    # set new name for multiple signatures/files
    multisig = utils.get_test_data('47+63-multisig.sig')
    other_sig = utils.get_test_data('2.fa.sig')
    c.run_sourmash('sig', 'rename', multisig, other_sig, 'fiz bar')

    # stdout should be new signature
    out = c.last_result.out

    n = 0
    for sig in load_signatures(out):
        assert sig.name == 'fiz bar'
        n += 1

    assert n == 9, n


@utils.in_tempdir
def test_sig_rename_1_multisig_ksize(c):
    # set new name for multiple signatures/files; select k=31
    multisig = utils.get_test_data('47+63-multisig.sig')
    other_sig = utils.get_test_data('2.fa.sig')
    c.run_sourmash('sig', 'rename', multisig, other_sig, 'fiz bar', '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    n = 0
    for sig in load_signatures(out):
        assert sig.name == 'fiz bar'
        n += 1

    assert n == 7, n


@utils.in_tempdir
def test_sig_rename_2_output_to_same(c):
    # change name of signature "in place", same output file
    sig47 = utils.get_test_data('47.fa.sig')
    inplace = c.output('inplace.sig')
    shutil.copyfile(sig47, inplace)

    print(inplace)

    c.run_sourmash('sig', 'rename', '-d', inplace, 'fiz bar', '-o', inplace)

    actual_rename_sig = sourmash.load_one_signature(inplace)
    assert actual_rename_sig.name == 'fiz bar'


@utils.in_tempdir
def test_sig_rename_3_file_dne(c):
    # rename on a file that does not exist should fail!
    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('sig', 'rename', 'no-such-sig', 'fiz bar')

    assert "Error while reading signatures from 'no-such-sig'" in c.last_result.err


@utils.in_tempdir
def test_sig_rename_3_file_dne_force(c):
    # rename on a file that does not exist should fail!
    c.run_sourmash('sig', 'rename', 'no-such-sig', 'fiz bar', '-f')
    print(c.last_result.err)

    assert "Error while reading signatures from 'no-such-sig'" in c.last_result.err


@utils.in_thisdir
def test_sig_cat_1(c):
    # cat 47 to 47...
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'cat', sig47)

    # stdout should be same signature
    out = c.last_result.out

    test_cat_sig = sourmash.load_one_signature(sig47)
    actual_cat_sig = sourmash.load_one_signature(out)

    assert actual_cat_sig == test_cat_sig


@utils.in_thisdir
def test_sig_cat_1_no_unique(c):
    # cat 47 to 47... twice
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'cat', sig47, sig47)

    # stdout should be same signature
    out = c.last_result.out

    test_cat_sig = sourmash.load_one_signature(sig47)
    actual_cat_sigs = load_signatures(out)

    for n, sig in enumerate(actual_cat_sigs):
        assert sig == test_cat_sig

    assert n == 1 # two signatures, but enumerate stops at 1.
    assert 'encountered 1 MinHashes multiple times' in c.last_result.err


@utils.in_thisdir
def test_sig_cat_1_unique(c):
    # cat 47 to 47... twice... and get unique
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'cat', sig47, sig47, '--unique')

    # stdout should be same signature
    out = c.last_result.out
    err = c.last_result.err

    test_cat_sig = sourmash.load_one_signature(sig47)
    actual_cat_sigs = load_signatures(out)

    for n, sig in enumerate(actual_cat_sigs):
        assert sig == test_cat_sig

    assert n == 0 # enumerate stops at 0, first sig.
    assert 'encountered 1 MinHashes multiple times' in err
    assert '...and removed the duplicates, because --unique was specified.' in err


@utils.in_thisdir
def test_sig_cat_2(c):
    # cat several
    sig47 = utils.get_test_data('47.fa.sig')
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'cat', sig47, sig47abund, multisig)

    # stdout should be same signatures
    out = c.last_result.out

    siglist = list(load_signatures(out))
    print(len(siglist))

    assert repr(siglist) == """[SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 57e2b22f), SourmashSignature('NC_009661.1 Shewanella baltica OS185 plasmid pS18501, complete sequence', bde81a41), SourmashSignature('NC_011663.1 Shewanella baltica OS223, complete genome', f033bbd8), SourmashSignature('NC_011664.1 Shewanella baltica OS223 plasmid pS22301, complete sequence', 87a9aec4), SourmashSignature('NC_011668.1 Shewanella baltica OS223 plasmid pS22302, complete sequence', 837bf2a7), SourmashSignature('NC_011665.1 Shewanella baltica OS223 plasmid pS22303, complete sequence', 485c3377)]"""


@utils.in_tempdir
def test_sig_cat_2_out(c):
    # cat several
    sig47 = utils.get_test_data('47.fa.sig')
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'cat', sig47, sig47abund, multisig,
                   '-o', 'out.sig')

    # stdout should be same signatures
    out = c.output('out.sig')

    siglist = list(load_signatures(out))
    print(len(siglist))

    assert repr(siglist) == """[SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 57e2b22f), SourmashSignature('NC_009661.1 Shewanella baltica OS185 plasmid pS18501, complete sequence', bde81a41), SourmashSignature('NC_011663.1 Shewanella baltica OS223, complete genome', f033bbd8), SourmashSignature('NC_011664.1 Shewanella baltica OS223 plasmid pS22301, complete sequence', 87a9aec4), SourmashSignature('NC_011668.1 Shewanella baltica OS223 plasmid pS22302, complete sequence', 837bf2a7), SourmashSignature('NC_011665.1 Shewanella baltica OS223 plasmid pS22303, complete sequence', 485c3377)]"""


@utils.in_tempdir
def test_sig_cat_2_out_inplace(c):
    # cat several; check that we can overwrite one of the input files.
    sig47 = utils.get_test_data('47.fa.sig')
    input_sig = c.output('inp.sig')
    shutil.copyfile(sig47, input_sig)

    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')

    # write out to input.
    c.run_sourmash('sig', 'cat', input_sig, sig47abund, multisig,
                   '-o', input_sig)

    # stdout should be same signatures
    out = input_sig

    siglist = list(load_signatures(out))
    print(len(siglist))

    assert repr(siglist) == """[SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 57e2b22f), SourmashSignature('NC_009661.1 Shewanella baltica OS185 plasmid pS18501, complete sequence', bde81a41), SourmashSignature('NC_011663.1 Shewanella baltica OS223, complete genome', f033bbd8), SourmashSignature('NC_011664.1 Shewanella baltica OS223 plasmid pS22301, complete sequence', 87a9aec4), SourmashSignature('NC_011668.1 Shewanella baltica OS223 plasmid pS22302, complete sequence', 837bf2a7), SourmashSignature('NC_011665.1 Shewanella baltica OS223 plasmid pS22303, complete sequence', 485c3377)]"""


@utils.in_tempdir
def test_sig_cat_3_filelist(c):
    # cat using a file list as input
    sig47 = utils.get_test_data('47.fa.sig')
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')

    filelist = c.output("filelist")
    with open(filelist, 'w') as f:
        f.write("\n".join((sig47, sig47abund, multisig)))

    c.run_sourmash('sig', 'cat', filelist,
                   '-o', 'out.sig')

    # stdout should be same signatures
    out = c.output('out.sig')

    # make this a list, not a set, because a set will collapse identical
    # signatures. `sig cat` does not collapse identical signatures, although
    # the pathlist function will ignore duplicate files.
    siglist = list(load_signatures(out))

    # verify the number of signatures matches what we expect to see based
    # on the input files
    all_sigs = []
    all_sigs += list(load_signatures(sig47))
    all_sigs += list(load_signatures(sig47abund))
    all_sigs += list(load_signatures(multisig))

    assert len(all_sigs) == len(siglist)

    # sort the signatures by something deterministic and unique
    siglist.sort(key = lambda x: x.md5sum())

    assert repr(siglist) == """[SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_011665.1 Shewanella baltica OS223 plasmid pS22303, complete sequence', 485c3377), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 57e2b22f), SourmashSignature('NC_011668.1 Shewanella baltica OS223 plasmid pS22302, complete sequence', 837bf2a7), SourmashSignature('NC_011664.1 Shewanella baltica OS223 plasmid pS22301, complete sequence', 87a9aec4), SourmashSignature('NC_009661.1 Shewanella baltica OS185 plasmid pS18501, complete sequence', bde81a41), SourmashSignature('NC_011663.1 Shewanella baltica OS223, complete genome', f033bbd8)]"""


@utils.in_tempdir
def test_sig_cat_4_filelist_with_dbs(c):
    # cat using a file list as input
    sig47 = utils.get_test_data('47.fa.sig')
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sbt = utils.get_test_data('v6.sbt.zip')

    filelist = c.output("filelist")
    with open(filelist, 'w') as f:
        f.write("\n".join((sig47, sig47abund, sbt)))

    c.run_sourmash('sig', 'cat', filelist,
                   '-o', 'out.sig')

    # stdout should be same signatures
    out = c.output('out.sig')

    siglist = list(load_signatures(out))
    print(len(siglist))
    # print("siglist: ",siglist)
    # print("\n")

    # verify the number of signatures matches what we expect to see based
    # on the input files
    all_sigs = []
    all_sigs += list(load_signatures(sig47))
    all_sigs += list(load_signatures(sig47abund))
    all_sigs += list(sourmash.load_file_as_signatures(sbt))

    assert len(all_sigs) == len(siglist)

    # sort the signatures by something deterministic and unique
    siglist.sort(key = lambda x: x.md5sum())

    assert repr(siglist) == """[SourmashSignature('', 0107d767), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('', 4e94e602), SourmashSignature('', 60f7e23c), SourmashSignature('', 6d6e87e1), SourmashSignature('', b59473c9), SourmashSignature('', f0c834bc), SourmashSignature('', f71e7817)]"""


@utils.in_tempdir
def test_sig_cat_5_from_file(c):
    # cat using a file list as input
    sig47 = utils.get_test_data('47.fa.sig')
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sbt = utils.get_test_data('v6.sbt.zip')

    filelist = c.output("filelist")
    with open(filelist, 'w') as f:
        f.write("\n".join((sig47, sig47abund, sbt)))

    c.run_sourmash('sig', 'cat', '--from-file', filelist,
                   '-o', 'out.sig')

    # stdout should be same signatures
    out = c.output('out.sig')

    siglist = list(load_signatures(out))
    print(len(siglist))
    # print("siglist: ",siglist)
    # print("\n")

    # verify the number of signatures matches what we expect to see based
    # on the input files
    all_sigs = []
    all_sigs += list(load_signatures(sig47))
    all_sigs += list(load_signatures(sig47abund))
    all_sigs += list(sourmash.load_file_as_signatures(sbt))

    assert len(all_sigs) == len(siglist)

    # sort the signatures by something deterministic and unique
    siglist.sort(key = lambda x: x.md5sum())

    assert repr(siglist) == """[SourmashSignature('', 0107d767), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691), SourmashSignature('', 4e94e602), SourmashSignature('', 60f7e23c), SourmashSignature('', 6d6e87e1), SourmashSignature('', b59473c9), SourmashSignature('', f0c834bc), SourmashSignature('', f71e7817)]"""


def test_sig_cat_5_from_file_picklist(runtmp):
    c = runtmp

    # cat using a file list as input
    sig47 = utils.get_test_data('47.fa.sig')
    sbt = utils.get_test_data('v6.sbt.zip')

    filelist = c.output("filelist")
    with open(filelist, 'w') as f:
        f.write("\n".join((sig47, sbt)))

    picklist = _write_file(runtmp, 'pl.csv', ['md5short', '09a08691'])

    c.run_sourmash('sig', 'cat', '--from-file', filelist,
                   '--picklist', f'{picklist}:md5short:md5short',
                   '-o', 'out.sig')

    # stdout should be same signatures
    out = c.output('out.sig')

    siglist = list(load_signatures(out))
    print(len(siglist))
    # print("siglist: ",siglist)
    # print("\n")

    # verify the number of signatures matches what we expect to see based
    # on the input files
    all_sigs = []
    all_sigs += list(load_signatures(sig47, ksize=31))

    assert len(all_sigs) == len(siglist)

    # sort the signatures by something deterministic and unique
    siglist.sort(key = lambda x: x.md5sum())

    assert repr(siglist) == """[SourmashSignature('NC_009665.1 Shewanella baltica OS185, complete genome', 09a08691)]"""


def test_sig_split_1(runtmp):
    c = runtmp
    # split 47 into 1 sig :)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'split', sig47)

    outname = '09a08691.k=31.scaled=1000.DNA.dup=0.47.fa.sig'

    assert os.path.exists(c.output(outname))

    test_split_sig = sourmash.load_one_signature(sig47)
    actual_split_sig = sourmash.load_one_signature(c.output(outname))

    assert actual_split_sig == test_split_sig


def test_sig_split_1_fromfile_picklist(runtmp):
    c = runtmp
    # split 47 into 1 sig :)
    sig47 = utils.get_test_data('47.fa.sig')

    from_file = _write_file(runtmp, 'list.txt', [sig47])
    picklist = _write_file(runtmp, 'pl.csv', ['md5short', '09a08691'])

    c.run_sourmash('sig', 'split', '--from-file', from_file,
                   '--picklist', f'{picklist}:md5short:md5short')

    outname = '09a08691.k=31.scaled=1000.DNA.dup=0.47.fa.sig'

    assert os.path.exists(c.output(outname))

    test_split_sig = sourmash.load_one_signature(sig47)
    actual_split_sig = sourmash.load_one_signature(c.output(outname))

    assert actual_split_sig == test_split_sig


@utils.in_tempdir
def test_sig_split_1_overwrite(c):
    # check message about overwriting
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'split', sig47)

    outname = '09a08691.k=31.scaled=1000.DNA.dup=0.47.fa.sig'
    assert os.path.exists(c.output(outname))

    c.run_sourmash('sig', 'split', sig47)

    err = c.last_result.err
    print(err)
    assert '** overwriting existing file ' + outname in err


@utils.in_tempdir
def test_sig_split_2(c):
    # split 47 twice
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'split', sig47, sig47)

    outname1 = '09a08691.k=31.scaled=1000.DNA.dup=0.47.fa.sig'
    outname2 = '09a08691.k=31.scaled=1000.DNA.dup=1.47.fa.sig'

    assert os.path.exists(c.output(outname1))
    assert os.path.exists(c.output(outname2))

    test_split_sig = sourmash.load_one_signature(sig47)

    actual_split_sig = sourmash.load_one_signature(c.output(outname1))
    assert actual_split_sig == test_split_sig

    actual_split_sig = sourmash.load_one_signature(c.output(outname2))
    assert actual_split_sig == test_split_sig


@utils.in_tempdir
def test_sig_split_2_outdir(c):
    # split 47 twice, put in outdir
    sig47 = utils.get_test_data('47.fa.sig')
    outdir = c.output('sigout/')
    c.run_sourmash('sig', 'split', sig47, sig47, '--outdir', outdir)

    outname1 = 'sigout/09a08691.k=31.scaled=1000.DNA.dup=0.47.fa.sig'
    outname2 = 'sigout/09a08691.k=31.scaled=1000.DNA.dup=1.47.fa.sig'

    assert os.path.exists(c.output(outname1))
    assert os.path.exists(c.output(outname2))

    test_split_sig = sourmash.load_one_signature(sig47)

    actual_split_sig = sourmash.load_one_signature(c.output(outname1))
    assert actual_split_sig == test_split_sig

    actual_split_sig = sourmash.load_one_signature(c.output(outname2))
    assert actual_split_sig == test_split_sig


@utils.in_tempdir
def test_sig_split_3_multisig(c):
    # split 47 and 47+63-multisig.sig
    sig47 = utils.get_test_data('47.fa.sig')
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'split', sig47, multisig)

    outlist = ['57e2b22f.k=31.scaled=1000.DNA.dup=0.none.sig',
               'bde81a41.k=31.scaled=1000.DNA.dup=0.none.sig',
               'f033bbd8.k=31.scaled=1000.DNA.dup=0.none.sig',
               '87a9aec4.k=31.scaled=1000.DNA.dup=0.none.sig',
               '837bf2a7.k=31.scaled=1000.DNA.dup=0.none.sig',
               '485c3377.k=31.scaled=1000.DNA.dup=0.none.sig']
    for filename in outlist:
        assert os.path.exists(c.output(filename))


@utils.in_tempdir
def test_sig_split_4_sbt_prot(c):
    # split sbt
    sbt1 = utils.get_test_data('prot/protein.sbt.zip')
    sbt2 = utils.get_test_data('prot/dayhoff.sbt.zip')
    sbt3 = utils.get_test_data('prot/hp.sbt.zip')
    c.run_sourmash('sig', 'split', sbt1, sbt2, sbt3)

    outlist = ['16869d2c.k=19.scaled=100.protein.dup=0.GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
               '120d311c.k=19.scaled=100.protein.dup=0.GCA_001593935.1_ASM159393v1_protein.faa.gz.sig',
               'fbca5e52.k=19.scaled=100.dayhoff.dup=0.GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
               '1cbd888b.k=19.scaled=100.dayhoff.dup=0.GCA_001593935.1_ASM159393v1_protein.faa.gz.sig',
               'ea2a1ad2.k=19.scaled=100.hp.dup=0.GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
               'bb0e6d90.k=19.scaled=100.hp.dup=0.GCA_001593935.1_ASM159393v1_protein.faa.gz.sig']
    for filename in outlist:
        assert os.path.exists(c.output(filename))


@utils.in_tempdir
def test_sig_split_4_lca_prot(c):
    # split lca
    lca1 = utils.get_test_data('prot/protein.lca.json.gz')
    lca2 = utils.get_test_data('prot/dayhoff.lca.json.gz')
    lca3 = utils.get_test_data('prot/hp.lca.json.gz')
    c.run_sourmash('sig', 'split', lca1, lca2, lca3)

    print(c.last_result.out)
    print(c.last_result.err)

    outlist = ['16869d2c.k=19.scaled=100.protein.dup=0.none.sig',
               '120d311c.k=19.scaled=100.protein.dup=0.none.sig',
               'fbca5e52.k=19.scaled=100.dayhoff.dup=0.none.sig',
               '1cbd888b.k=19.scaled=100.dayhoff.dup=0.none.sig',
               'ea2a1ad2.k=19.scaled=100.hp.dup=0.none.sig',
               'bb0e6d90.k=19.scaled=100.hp.dup=0.none.sig']
    for filename in outlist:
        assert os.path.exists(c.output(filename))


@utils.in_tempdir
def test_sig_split_5_no_exist(c):
    # no such file
    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('sig', 'split', 'foo')


def test_sig_split_6_numsigs(runtmp):
    c = runtmp

    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')
    c.run_sourmash('sig', 'split', sigs11)

    print(c.last_result.out)
    print(c.last_result.err)

    outlist = ['1437d8ea.k=21.num=500.DNA.dup=0.genome-s11.fa.gz.sig',
               '37aea787.k=7.num=500.protein.dup=0.genome-s11.fa.gz.sig',
               '68c565be.k=30.num=500.DNA.dup=0.genome-s11.fa.gz.sig',
               '73b6df1c.k=10.num=500.protein.dup=0.genome-s11.fa.gz.sig']

    for filename in outlist:
        assert os.path.exists(c.output(filename))


def test_sig_extract_1(runtmp):
    c = runtmp

    # extract 47 from 47... :)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'extract', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_1(runtmp):
    c = runtmp

    # extract 47 from 47... :)
    sig47 = utils.get_test_data('47.fa.sig')
    from_file = _write_file(runtmp, 'list.txt', [sig47])
    c.run_sourmash('sig', 'extract', '--from-file', from_file)

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
    with pytest.raises(SourmashCommandFailed) as exc:
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
    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('sig', 'extract', sig47, '--name', 'FOO')


@utils.in_tempdir
def test_sig_extract_6(c):
    # extract matches to several names from among several signatures
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 2


@utils.in_tempdir
def test_sig_extract_7(c):
    # extract matches based on ksize
    sig2 = utils.get_test_data('2.fa.sig')
    c.run_sourmash('sig', 'extract', sig2, '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1


@utils.in_tempdir
def test_sig_extract_7_no_ksize(c):
    # extract all three matches when -k not specified
    sig2 = utils.get_test_data('2.fa.sig')
    c.run_sourmash('sig', 'extract', sig2)

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 3


def test_sig_extract_8_picklist_md5(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5full:md5"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

    err = runtmp.last_result.err

    print(err)
    assert "loaded 1 distinct values into picklist." in err
    assert "loaded 1 total that matched ksize & molecule type" in err
    assert "extracted 1 signatures from 2 file(s)" in err
    assert "for given picklist, found 1 matches to 1 distinct values" in err

def test_sig_extract_8_picklist_md5_include(runtmp):
    # extract 47 from 47, using a picklist w/full md5:: explicit include
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5full:md5:include"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

    err = runtmp.last_result.err

    print(err)
    assert "loaded 1 distinct values into picklist." in err
    assert "loaded 1 total that matched ksize & molecule type" in err
    assert "extracted 1 signatures from 2 file(s)" in err
    assert "for given picklist, found 1 matches to 1 distinct values" in err


def test_sig_extract_8_picklist_md5_exclude(runtmp):
    # extract 63 from 47,63 by excluding 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5full:md5:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

    err = runtmp.last_result.err

    print(err)
    assert "loaded 1 distinct values into picklist." in err
    assert "loaded 1 total that matched ksize & molecule type" in err
    assert "extracted 1 signatures from 2 file(s)" in err
    assert "for given picklist, found 1 matches by excluding 1 distinct values" in err


def test_sig_extract_8_picklist_md5_require_all(runtmp):
    # extract 47 from 47, using a picklist w/full md5;
    # confirm that check missing picklist val errors out on --picklist-require
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)
        w.writerow(dict(exactName='', md5full='BAD MD5',
                        md5short='', fullIdent='', nodotIdent=''))

    picklist_arg = f"{picklist_csv}:md5full:md5"
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sig47, sig63,
                        '--picklist', picklist_arg,
                        '--picklist-require-all')

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

    err = runtmp.last_result.err

    print(err)
    assert "loaded 2 distinct values into picklist." in err
    assert "loaded 1 total that matched ksize & molecule type" in err
    assert "extracted 1 signatures from 2 file(s)" in err
    assert "for given picklist, found 1 matches to 2 distinct values" in err
    assert 'WARNING: 1 missing picklist values.' in err
    assert 'ERROR: failing because --picklist-require-all was set' in err


def test_sig_extract_8_picklist_name(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:exactName:name"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_name_exclude(runtmp):
    # exclude 47 based on name
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:exactName:name:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_ident(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:fullIdent:ident"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_ident_exclude(runtmp):
    # exclude 47 based on ident
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:fullIdent:ident:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_ident_dot(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:nodotIdent:identprefix"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_ident_dot_exclude(runtmp):
    # exlude 47 based on identprefix
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:nodotIdent:identprefix:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_md5_short(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5prefix8"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_md5_short_exclude(runtmp):
    # exclude 47 based on md5prefix8
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5prefix8:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_md5_short_alias(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_md5_short_alias_exclude(runtmp):
    # exlude 47 based on md5prefix8 alias, md5short
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short:exclude"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig63)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_sig_extract_8_picklist_md5_short_alias_with_md5_selector_nomatch(runtmp):
    # extract 47 from 47, using a picklist w/full md5 and also md5 selector
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short"
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sig47, sig63,
                        '--picklist', picklist_arg,
                        '--md5', 'XXX') # no match to md5 selector here

    err = runtmp.last_result.err
    assert "no matching signatures to save!" in err


def test_sig_extract_8_picklist_md5_short_alias_with_md5_selector_nomatch_exclude(runtmp):
    # exclude 47 using a picklist w/full md5 and also md5 selector
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short:exclude"
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sig47, sig63,
                        '--picklist', picklist_arg,
                        '--md5', 'XXX') # no match to md5 selector here

    err = runtmp.last_result.err
    assert "no matching signatures to save!" in err


def test_sig_extract_8_picklist_md5_short_alias_with_md5_selector(runtmp):
    # extract 47 from 47, using a picklist w/full md5 and also md5 selector
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short"
    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg,
                    '--md5', '09a08691ce5295215')

    # stdout should be new signature
    out = runtmp.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

def test_sig_extract_8_picklist_md5_short_alias_with_md5_selector_exclude(runtmp):
    # exclude 47, using a picklist w/full md5; but try to select with md5 selector
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='09a08691ce52952152f0e866a59f6261',
               md5short='09a08691ce5295215',
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5short:md5short:exclude"
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist', picklist_arg,
                    '--md5', '09a08691ce5295215')

    # NTP: do we want to emit a more informative "conflicting selectors" type of msg?
    err = runtmp.last_result.err
    print(err)
    assert "loaded 1 distinct values into picklist." in err
    assert "loaded 1 total that matched ksize & molecule type" in err
    assert 'no matching signatures to save!' in err


def test_sig_extract_8_picklist_md5_nomatch(runtmp):
    # use an empty picklist => no match
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5short'])
        w.writeheader()

    picklist_arg = f"{picklist_csv}:md5short:md5prefix8"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist',
                        picklist_arg)

    # stdout should be new signature
    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)
    assert "no matching signatures to save!" in err
    assert runtmp.last_result.status != 0


def test_sig_extract_8_picklist_md5_nomatch_exclude(runtmp):
    # use an empty picklist to exclude => no match => include everything
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5short'])
        w.writeheader()

    picklist_arg = f"{picklist_csv}:md5short:md5prefix8:exclude"

    runtmp.sourmash('sig', 'extract', sig47, sig63, '--picklist',
                        picklist_arg)

    # stdout should be both signatures
    out = runtmp.last_result.out
    extract_siglist = list(load_signatures(out))
    print(len(extract_siglist))
    s47 = sourmash.load_file_as_signatures(sig47)
    s63 = sourmash.load_file_as_signatures(sig63)
    actual_extract_siglist = list(s47) + list(s63)

    assert set(extract_siglist) == set(actual_extract_siglist)

    err = runtmp.last_result.err
    print(err)
    assert runtmp.last_result.status == 0
    assert 'loaded 0 distinct values into picklist.' in err
    assert 'loaded 2 total that matched ksize & molecule type' in err
    assert 'extracted 2 signatures from 2 file(s)' in err
    assert 'for given picklist, found 2 matches by excluding 0 distinct values' in err


def test_sig_extract_9_picklist_md5_ksize_hp_select(runtmp):
    # test with -k and moltype selector
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:md5:md5"

    runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                    picklist_arg, '-k', '19', '--hp')

    # stdout should be new signature
    out = runtmp.last_result.out
    actual_extract_sig = sourmash.load_one_signature(out)

    print(actual_extract_sig.md5sum)
    assert str(actual_extract_sig) == 'GCA_001593925'
    assert actual_extract_sig.md5sum() == 'ea2a1ad233c2908529d124a330bcb672'
    assert actual_extract_sig.minhash.ksize == 19
    assert actual_extract_sig.minhash.moltype == 'hp'


def test_sig_extract_9_picklist_md5_ksize_hp_select_exclude(runtmp):
    # test picklist exclude with -k and moltype selector
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:md5:md5:exclude"

    runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                    picklist_arg, '-k', '19', '--hp')

    # stdout should be new signature
    out = runtmp.last_result.out
    actual_extract_sig = sourmash.load_one_signature(out)
    print(actual_extract_sig.md5sum)

    assert str(actual_extract_sig) == 'GCA_001593935'
    assert actual_extract_sig.md5sum() == 'bb0e6d90df01b7bd5d0956a5f9e3ed12'
    assert actual_extract_sig.minhash.ksize == 19
    assert actual_extract_sig.minhash.moltype == 'hp'


def test_sig_extract_10_picklist_md5_dups_and_empty(runtmp):
    # test empty picklist values, and duplicate picklist values
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))
        w.writerow(dict(md5=''))

    picklist_arg = f"{picklist_csv}:md5:md5"

    runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                    picklist_arg, '-k', '19', '--hp')

    # stdout should be new signature
    out = runtmp.last_result.out
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig.minhash.ksize == 19
    assert actual_extract_sig.minhash.moltype == 'hp'
    assert actual_extract_sig.md5sum() == 'ea2a1ad233c2908529d124a330bcb672'

    err = runtmp.last_result.err
    print(err)

    assert "WARNING: 1 empty values in column 'md5' in picklist file" in err
    assert "WARNING: 1 values in picklist column 'md5' were not distinct" in err


def test_sig_extract_10_picklist_md5_dups_and_empty_exclude(runtmp):
    # test empty picklist values, and duplicate picklist values for exclude
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))
        w.writerow(dict(md5=''))

    picklist_arg = f"{picklist_csv}:md5:md5:exclude"

    runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                    picklist_arg, '-k', '19', '--hp')

    # stdout should be new signature
    out = runtmp.last_result.out
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig.minhash.ksize == 19
    assert actual_extract_sig.minhash.moltype == 'hp'
    assert actual_extract_sig.md5sum() == 'bb0e6d90df01b7bd5d0956a5f9e3ed12'

    err = runtmp.last_result.err
    print(err)

    assert "WARNING: 1 empty values in column 'md5' in picklist file" in err
    assert "WARNING: 1 values in picklist column 'md5' were not distinct" in err


def test_sig_extract_11_picklist_bad_coltype(runtmp):
    # test with invalid picklist coltype
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:md5:BADCOLTYPE"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "invalid picklist column type 'BADCOLTYPE'" in err


def test_sig_extract_11_picklist_bad_coltype_exclude(runtmp):
    # test with invalid picklist coltype
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:md5:BADCOLTYPE:exclude"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "invalid picklist column type 'BADCOLTYPE'" in err


def test_sig_extract_12_picklist_bad_argstr(runtmp):
    # test with invalid argument format to --picklist
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "invalid picklist argument" in err


def test_sig_extract_12_picklist_bad_pickstyle(runtmp):
    # test with invalid argument format to --picklist
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:md5:md5:XXX"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "invalid picklist 'pickstyle' argument, 'XXX': must be 'include' or 'exclude'" in err


def test_sig_extract_12_picklist_bad_colname(runtmp):
    # test with invalid picklist colname
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:BADCOLNAME:md5"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "column 'BADCOLNAME' not in pickfile" in err


def test_sig_extract_12_picklist_bad_colname_exclude(runtmp):
    # test with invalid picklist colname
    sigdir = utils.get_test_data('prot/')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='ea2a1ad233c2908529d124a330bcb672'))

    picklist_arg = f"{picklist_csv}:BADCOLNAME:md5:exclude"

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'extract', sigdir, '--picklist',
                        picklist_arg, '-k', '19', '--hp')

    err = runtmp.last_result.err
    print(err)
    assert "column 'BADCOLNAME' not in pickfile" in err


def test_sig_flatten_1(runtmp):
    c = runtmp

    # extract matches to several names from among several signatures & flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'flatten', sig47abund, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


def test_sig_flatten_1(runtmp):
    c = runtmp

    # extract matches to several names from among several signatures & flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')

    from_file = _write_file(runtmp, 'list.txt', [sig47abund])
    picklist = _write_file(runtmp, 'pl.csv', ['md5short', '09a08691'])

    c.run_sourmash('sig', 'flatten', '--from-file', from_file,
                   '--picklist', f'{picklist}:md5short:md5short')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


@utils.in_tempdir
def test_sig_flatten_1_select_name(c):
    # extract matches to several names from among several signatures & flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'flatten', sig2, sig47abund, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


def test_sig_flatten_1_select_md5(runtmp):
    c = runtmp

    # extract matches to several names from among several signatures & flatten
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'flatten', sig2, sig47abund, '--md5', '09a08691c')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


def test_sig_flatten_2_ksize(runtmp):
    c = runtmp
    # flatten only one signature selected using ksize
    psw_mag = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    c.run_sourmash('sig', 'flatten', psw_mag, '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    siglist = load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1


@utils.in_tempdir
def test_sig_downsample_1_scaled(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sig47)
    actual_downsample_sig = sourmash.load_one_signature(out)

    test_mh = test_downsample_sig.minhash.downsample(scaled=10000)

    assert actual_downsample_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_downsample_1_scaled_downsample_multisig(c):
    # downsample many scaled signatures in one file
    multisig = utils.get_test_data('47+63-multisig.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000', multisig)

    # stdout should be new signatures
    out = c.last_result.out

    for sig in load_signatures(out):
        assert sig.minhash.scaled == 10000


@utils.in_tempdir
def test_sig_downsample_1_scaled_to_num(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'downsample', '--num', '500', sig47)

    # stdout should be new signature
    out = c.last_result.out

    actual_downsample_sig = sourmash.load_one_signature(out)
    actual_mins = actual_downsample_sig.minhash.hashes.keys()
    actual_mins = list(actual_mins)
    actual_mins.sort()

    test_downsample_sig = sourmash.load_one_signature(sig47)
    test_mins = test_downsample_sig.minhash.hashes.keys()
    test_mins = list(test_mins)
    test_mins.sort()
    test_mins = test_mins[:500]           # take 500 smallest

    assert actual_mins == test_mins


def test_sig_downsample_check_num_bounds_negative(runtmp):
    c=runtmp
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '--num', '-5', sig47)

    assert "ERROR: num value must be positive" in c.last_result.err


def test_sig_downsample_check_num_bounds_less_than_minimum(runtmp):
    c=runtmp
    sig47 = utils.get_test_data('47.fa.sig')

    c.run_sourmash('sig', 'downsample', '--num', '25', sig47)

    assert "WARNING: num value should be >= 50. Continuing anyway." in c.last_result.err


def test_sig_downsample_check_num_bounds_more_than_maximum(runtmp):
    c=runtmp
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '--num', '100000', sig47)

    assert "WARNING: num value should be <= 50000. Continuing anyway." in c.last_result.err


@utils.in_tempdir
def test_sig_downsample_1_scaled_to_num_fail(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '--num', '50000', sig47)


@utils.in_tempdir
def test_sig_downsample_1_scaled_empty(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
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
    test_mh = test_downsample_sig.minhash.downsample(num=500)

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

    test_mins = test_downsample_sig.minhash.hashes.keys()
    actual_mins = actual_downsample_sig.minhash.hashes.keys()

    # select those mins that are beneath the new max hash...
    max_hash = actual_downsample_sig.minhash._max_hash
    test_mins_down = { k for k in test_mins if k < max_hash }
    assert test_mins_down == set(actual_mins)


@utils.in_tempdir
def test_sig_downsample_2_num_to_scaled_fail(c):
    # downsample a num signature and FAIL to convert it into a scaled sig
    # because new scaled is too low
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '--scaled', '100',
                       '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_downsample_2_num_and_scaled_both_fail(c):
    # cannot specify both --num and --scaled
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '--scaled', '100', '--num', '50',
                       '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_downsample_2_num_empty(c):
    # downsample a num signature
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sig', 'downsample', '-k', '21', '--dna', sigs11)


def test_sig_describe_1(runtmp):
    c = runtmp

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


def test_sig_describe_1_fromfile_picklist(runtmp):
    c = runtmp

    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')
    from_file = _write_file(runtmp, 'list.txt', [sig47])
    picklist = _write_file(runtmp, 'pl.csv', ['md5short', '09a08691'])

    c.run_sourmash('sig', 'describe',  '--from-file', from_file,
                   '--picklist', f'{picklist}:md5short:md5short')

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


@utils.in_thisdir
def test_sig_describe_protein(c):
    # test describe on a singleton protein signature
    testdata = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    c.run_sourmash('sig', 'describe', testdata)

    assert 'k=19 molecule=protein num=0 scaled=100 seed=42 track_abundance=0' in c.last_result.out


@utils.in_thisdir
def test_sig_describe_hp(c):
    # test describe on a singleton hp signature
    testdata = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    c.run_sourmash('sig', 'describe', testdata)

    assert 'k=19 molecule=hp num=0 scaled=100 seed=42 track_abundance=0' in c.last_result.out


@utils.in_thisdir
def test_sig_describe_dayhoff(c):
    # test describe on a singleton dayhoff signature
    testdata = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    c.run_sourmash('sig', 'describe', testdata)

    assert 'k=19 molecule=dayhoff num=0 scaled=100 seed=42 track_abundance=0' in c.last_result.out


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
signature: ** no name **
source file: short.fa
md5: e45a080101751e044d6df861d3d0f3fd
k=7 molecule=protein num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: c027e96c3379d38942639219daa24fdc
k=7 molecule=dayhoff num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 4b50ae79657d9dd07a1d543ba8b986a0
k=7 molecule=hp num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 1136a8a68420bd93683e45cdaf109b80
k=21 molecule=DNA num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 4244d1612598af044e799587132f007e
k=10 molecule=protein num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 396dcb7c1875f48ca31e0759bec72ee1
k=10 molecule=dayhoff num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 4c43878296459783dbd6a4a071ab7e9d
k=10 molecule=hp num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

---
signature filename: short.fa.sig
signature: ** no name **
source file: short.fa
md5: 71f7c111c01785e5f38efad45b00a0e1
k=30 molecule=DNA num=500 scaled=0 seed=42 track_abundance=0
size: 500
signature license: CC0

""".splitlines()
    for line in out.splitlines():
        cleaned_line = line.strip().replace(
            testdata_dirname, '').replace(location, '')
        assert cleaned_line in expected_output, cleaned_line


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
def test_sig_describe_1_sbt(c):
    # get basic info on multiple signatures in an SBT
    sigs = utils.get_test_data('prot/protein.sbt.zip')
    c.run_sourmash('sig', 'describe', sigs)

    out = c.last_result.out
    print(c.last_result)

    expected_output = """\
signature: GCA_001593925
signature: GCA_001593935
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@utils.in_tempdir
def test_sig_describe_1_lca(c):
    # get basic info on multiple signatures in an LCA database
    sigs = utils.get_test_data('prot/protein.lca.json.gz')
    c.run_sourmash('sig', 'describe', sigs)

    out = c.last_result.out
    print(c.last_result)

    expected_output = """\
signature: GCA_001593925
signature: GCA_001593935
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@utils.in_tempdir
def test_sig_describe_1_dir(c):
    # get basic info on multiple signatures in a directory
    sigs = utils.get_test_data('prot/protein/')
    c.run_sourmash('sig', 'describe', sigs)

    out = c.last_result.out
    print(c.last_result)

    # make sure signature names, as well as full path to .sig file under
    # directory, show up in output.
    expected_output = """\
signature: GCA_001593925
signature: GCA_001593935
prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig
prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@utils.in_tempdir
def test_sig_describe_1_zipfile(c):
    # get basic info on multiple signatures in a zipfile
    sigs = utils.get_test_data('prot/all.zip')
    c.run_sourmash('sig', 'describe', sigs)

    out = c.last_result.out
    print(c.last_result)

    expected_output = """\
k=19 molecule=dayhoff num=0 scaled=100 seed=42 track_abundance=0
k=19 molecule=dayhoff num=0 scaled=100 seed=42 track_abundance=0
k=19 molecule=hp num=0 scaled=100 seed=42 track_abundance=0
k=19 molecule=hp num=0 scaled=100 seed=42 track_abundance=0
k=19 molecule=protein num=0 scaled=100 seed=42 track_abundance=0
k=19 molecule=protein num=0 scaled=100 seed=42 track_abundance=0
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@utils.in_thisdir
def test_sig_describe_stdin(c):
    sig = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    with open(sig, 'rt') as fp:
        data = fp.read()

    c.run_sourmash('sig', 'describe', '-', stdin_data=data)

    assert 'signature: GCA_001593925' in c.last_result.out


@utils.in_tempdir
def test_sig_describe_empty(c):
    sig = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')

    ss = sourmash.load_file_as_signatures(sig)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]

    ss.name = ''
    ss.filename = ''

    outsig = c.output('xxx.sig')
    with open(outsig, 'wt') as fp:
        sourmash.save_signatures([ss], fp)

    ss = sourmash.load_file_as_signatures(outsig)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert ss.name == ''
    assert ss.filename == ''

    c.run_sourmash('sig', 'describe', outsig)
    print(c.last_result.out)
    assert 'signature: ** no name **' in c.last_result.out
    assert 'source file: ** no name **' in c.last_result.out


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
def test_import_export_1_by_md5(c):
    # check to make sure we can import what we've exported!
    inp = utils.get_test_data('genome-s11.fa.gz.sig')
    outp = c.output('export.json')

    c.run_sourmash('sig', 'export', inp, '-o', outp, '--md5', '1437d8eae6')
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


def test_import_mash_csv_to_sig(runtmp):
    # test copied over from 'sourmash import_csv'.
    testdata1 = utils.get_test_data('short.fa.msh.dump')
    testdata2 = utils.get_test_data('short.fa')

    runtmp.sourmash('sig', 'import', '--csv', testdata1, '-o', 'xxx.sig')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,num=970', testdata2)

    runtmp.sourmash('search', '-k', '31', 'short.fa.sig', 'xxx.sig')

    print("RUNTEMP", runtmp)

    assert '1 matches:' in runtmp.last_result.out
    assert '100.0%       short.fa' in runtmp.last_result.out


def test_sig_manifest_1_zipfile(runtmp):
    # make a manifest from a .zip file
    protzip = utils.get_test_data('prot/protein.zip')
    runtmp.sourmash('sig', 'manifest', protzip, '-o', 'SOURMASH-MANIFEST.csv')

    manifest_fn = runtmp.output('SOURMASH-MANIFEST.csv')
    with open(manifest_fn, newline='') as csvfp:
        manifest = CollectionManifest.load_from_csv(csvfp)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_manifest_2_sigfile(runtmp):
    # make a manifest from a .sig file
    sigfile = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')

    runtmp.sourmash('sig', 'manifest', sigfile, '-o', 'SOURMASH-MANIFEST.csv')

    status = runtmp.last_result.status
    out = runtmp.last_result.out
    err = runtmp.last_result.err

    manifest_fn = runtmp.output('SOURMASH-MANIFEST.csv')
    with open(manifest_fn, newline='') as csvfp:
        manifest = CollectionManifest.load_from_csv(csvfp)

    assert len(manifest) == 1
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list


def test_sig_manifest_3_sbt(runtmp):
    # make a manifest from an SBT
    protzip = utils.get_test_data('prot/protein.sbt.zip')
    runtmp.sourmash('sig', 'manifest', protzip, '-o', 'SOURMASH-MANIFEST.csv')

    manifest_fn = runtmp.output('SOURMASH-MANIFEST.csv')
    with open(manifest_fn, newline='') as csvfp:
        manifest = CollectionManifest.load_from_csv(csvfp)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_manifest_4_lca(runtmp):
    # make a manifest from a .lca.json file
    sigfile = utils.get_test_data('prot/protein.lca.json.gz')
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'manifest', sigfile, '-o',
                        'SOURMASH-MANIFEST.csv')

    status = runtmp.last_result.status
    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert status != 0
    assert "ERROR: manifests cannot be generated for this file." in err


def test_sig_manifest_5_dir(runtmp):
    # make a manifest from a directory
    sigfile = utils.get_test_data('prot/protein/')
    runtmp.sourmash('sig', 'manifest', sigfile, '-o', 'SOURMASH-MANIFEST.csv')

    status = runtmp.last_result.status
    out = runtmp.last_result.out
    err = runtmp.last_result.err

    manifest_fn = runtmp.output('SOURMASH-MANIFEST.csv')
    with open(manifest_fn, newline='') as csvfp:
        manifest = CollectionManifest.load_from_csv(csvfp)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_manifest_6_pathlist(runtmp):
    # make a manifest from a pathlist file
    sigfiles = utils.get_test_data('prot/protein/*.sig')
    sigfiles = glob.glob(sigfiles)

    pathlist = runtmp.output('pathlist.txt')
    with open(pathlist, 'wt') as fp:
        fp.write("\n".join(sigfiles))

    runtmp.sourmash('sig', 'manifest', pathlist, '-o', 'SOURMASH-MANIFEST.csv')

    status = runtmp.last_result.status
    out = runtmp.last_result.out
    err = runtmp.last_result.err

    manifest_fn = runtmp.output('SOURMASH-MANIFEST.csv')
    with open(manifest_fn, newline='') as csvfp:
        manifest = CollectionManifest.load_from_csv(csvfp)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_kmers_1_dna(runtmp):
    # test sig kmers on dna
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'dna', seqfile, '-p', 'scaled=1')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'DNA'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 970' in err
    assert 'found 970 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 970

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_1_dna_more_in_query(runtmp):
    # test sig kmers on dna, where query has more than matches
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'dna', seqfile, '-p', 'scaled=1')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'DNA'

    # make a new sequence for query, with more k-mers
    query_seqfile = runtmp.output('query.fa')
    with open(query_seqfile, 'wt') as fp:
        for record in screed.open(seqfile):
            fp.write(f">{record.name}\n{record.sequence}AGTTACGATC\n")

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', query_seqfile)

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 970' in err
    # should only find 970 overlapping hashes here --
    assert 'found 970 distinct matching hashes (100.0%)' in err


def test_sig_kmers_1_dna_empty_seq(runtmp):
    # test sig kmers with empty query seq
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'dna', seqfile, '-p', 'scaled=1')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'DNA'

    # make a new sequence for query, with more k-mers
    query_seqfile = runtmp.output('query.fa')
    with open(query_seqfile, 'wt') as fp:
        pass

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                        '--seq', query_seqfile)

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert "ERROR: no sequences searched!?" in err


def test_sig_kmers_1_dna_empty_sig(runtmp):
    # test sig kmers with empty query sig
    seqfile = utils.get_test_data('short.fa')

    mh = sourmash.MinHash(ksize=31, n=0, scaled=1)
    ss = sourmash.SourmashSignature(mh, name="empty")
    with open(runtmp.output('empty.sig'), 'wt') as fp:
        sourmash.save_signatures([ss], fp)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'kmers', '--sig', 'empty.sig',
                        '--seq', seqfile)

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert "ERROR: no hashes in query signature!?" in err


def test_sig_kmers_1_dna_single_sig(runtmp):
    # test sig kmers with a fabricated query sig with a single hash
    seqfile = utils.get_test_data('short.fa')

    mh = sourmash.MinHash(ksize=31, n=0, scaled=1)
    mh.add_hash(1070961951490202715)
    ss = sourmash.SourmashSignature(mh, name="small")
    with open(runtmp.output('small.sig'), 'wt') as fp:
        sourmash.save_signatures([ss], fp)

    runtmp.sourmash('sig', 'kmers', '--sig', 'small.sig',
                    '--seq', seqfile)

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1' in err
    assert 'found 1 distinct matching hashes (100.0%)' in err


def test_sig_kmers_1_dna_lowscaled(runtmp):
    # test sig kmers on dna with a scaled of 100, so not all k-mers
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'dna', seqfile, '-p', 'scaled=100')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'DNA'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 5' in err
    assert 'found 5 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 5

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_1_dna_num(runtmp):
    # test sig kmers on dna with a scaled of 100, so not all k-mers
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'dna', seqfile, '-p', 'num=50')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'DNA'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 50' in err
    assert 'found 50 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 50

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_1_dna_translate_protein(runtmp):
    # test sig kmers on dna
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'translate', seqfile, '-p', 'scaled=1')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'protein'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa', '--translate')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1942' in err
    assert 'found 1942 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1942

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_1_dna_translate_dayhoff(runtmp):
    # test sig kmers on dna
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'translate', seqfile, '-p', 'scaled=1,dayhoff')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'dayhoff'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa', '--translate')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1906' in err
    assert 'found 1906 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1906

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_1_dna_translate_hp(runtmp):
    # test sig kmers on dna
    seqfile = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch', 'translate', seqfile, '-p', 'scaled=1,hp')
    ss = sourmash.load_one_signature(runtmp.output('short.fa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'hp'

    runtmp.sourmash('sig', 'kmers', '--sig', 'short.fa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'short.csv',
                    '--save-sequences', 'matched.fa', '--translate')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1750' in err
    assert 'found 1750 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 1
    assert len(records[0].sequence) == 1000, len(records[0].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_sequence(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('short.csv'))
    with open(runtmp.output('short.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1750

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_sequence(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_2_protein(runtmp):
    # test out sig kmers on an faa file
    seqfile = utils.get_test_data('ecoli.faa')

    runtmp.sourmash('sketch', 'protein', seqfile, '-p', 'scaled=1')
    ss = sourmash.load_one_signature(runtmp.output('ecoli.faa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'protein'

    runtmp.sourmash('sig', 'kmers', '--sig', 'ecoli.faa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'ecoli.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1112' in err
    assert 'found 1112 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 2
    assert len(records[0].sequence) == 820, len(records[0].sequence)
    assert len(records[1].sequence) == 310, len(records[1].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_protein(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('ecoli.csv'))
    with open(runtmp.output('ecoli.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1112

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_protein(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_2_dayhoff(runtmp):
    # test out sig kmers on an faa file
    seqfile = utils.get_test_data('ecoli.faa')

    runtmp.sourmash('sketch', 'protein', seqfile, '-p', 'scaled=1,dayhoff')
    ss = sourmash.load_one_signature(runtmp.output('ecoli.faa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'dayhoff'

    runtmp.sourmash('sig', 'kmers', '--sig', 'ecoli.faa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'ecoli.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1100' in err
    assert 'found 1100 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 2
    assert len(records[0].sequence) == 820, len(records[0].sequence)
    assert len(records[1].sequence) == 310, len(records[1].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_protein(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('ecoli.csv'))
    with open(runtmp.output('ecoli.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1100

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_protein(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0


def test_sig_kmers_2_hp(runtmp):
    # test out sig kmers on an faa file
    seqfile = utils.get_test_data('ecoli.faa')

    runtmp.sourmash('sketch', 'protein', seqfile, '-p', 'scaled=1,hp')
    ss = sourmash.load_one_signature(runtmp.output('ecoli.faa.sig'))
    mh = ss.minhash
    assert mh.moltype == 'hp'

    runtmp.sourmash('sig', 'kmers', '--sig', 'ecoli.faa.sig',
                    '--seq', seqfile,
                    '--save-kmers', 'ecoli.csv',
                    '--save-sequences', 'matched.fa')

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)

    assert 'total hashes in merged signature: 1048' in err
    assert 'found 1048 distinct matching hashes (100.0%)' in err

    # check FASTA output
    assert os.path.exists(runtmp.output('matched.fa'))
    records = list(screed.open(runtmp.output('matched.fa')))
    assert len(records) == 2
    assert len(records[0].sequence) == 820, len(records[0].sequence)
    assert len(records[1].sequence) == 310, len(records[1].sequence)

    seq_mh = mh.copy_and_clear()
    for record in records:
        seq_mh.add_protein(record.sequence)
    assert seq_mh.similarity(mh) == 1.0

    # check CSV output w/k-mers and hashes etc
    assert os.path.exists(runtmp.output('ecoli.csv'))
    with open(runtmp.output('ecoli.csv'), newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1048

    check_mh = mh.copy_and_clear()
    check_mh2 = mh.copy_and_clear()
    for row in rows:
        check_mh.add_protein(row['kmer'])
        check_mh2.add_hash(int(row['hashval']))
    assert check_mh.similarity(mh) == 1.0
    assert check_mh2.similarity(mh) == 1.0
