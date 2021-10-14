"""
Tests for functions in sourmash_args module.
"""
import os
import csv
import pytest
import gzip
import zipfile
import io
import contextlib

import sourmash_tst_utils as utils
import sourmash
from sourmash import sourmash_args


def test_save_signatures_api_none():
    # save to sigfile
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    with sourmash_args.SaveSignaturesToLocation(None) as save_sig:
        print(repr(save_sig))
        save_sig.add(ss2)
        save_sig.add(ss47)

    # nothing to test - no output!


def test_save_signatures_to_location_1_sig(runtmp):
    # save to sigfile.sig
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.sig')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_stdout():
    # save to stdout
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    output_capture = io.StringIO()
    with contextlib.redirect_stdout(output_capture):
        with sourmash_args.SaveSignaturesToLocation("-") as save_sig:
            save_sig.add(ss2)
            save_sig.add(ss47)

    output = output_capture.getvalue()

    saved = list(sourmash.signature.load_signatures(output))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_sig_is_default(runtmp):
    # save to sigfile.txt
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.txt')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

    saved = list(sourmash.signature.load_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_sig_gz(runtmp):
    # save to sigfile.gz
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.sig.gz')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

    # can we open as a .gz file?
    with gzip.open(outloc, "r") as fp:
        print(save_sig)
        fp.read()

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_zip(runtmp):
    # save to sigfile.zip
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.zip')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

    # can we open as a .zip file?
    with zipfile.ZipFile(outloc, "r") as zf:
        assert list(zf.infolist())

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_zip_dup(runtmp):
    # save to sigfile.zip
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.zip')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

        # here we have to change the names so the sig content is different;
        # exactly duplicates will not be saved, otherwise.
        ss2.name = 'different name for ss2'
        save_sig.add(ss2)
        ss47.name = 'different name for ss47'
        save_sig.add(ss47)

    # can we open as a .zip file?
    with zipfile.ZipFile(outloc, "r") as zf:
        assert list(zf.infolist())

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 4


def test_save_signatures_to_location_1_dirout(runtmp):
    # save to sigout/ (directory)
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('sigout/')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)

    assert os.path.isdir(outloc)

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 2


def test_save_signatures_to_location_1_dirout_duplicate(runtmp):
    # save to sigout/ (directory)
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('sigout/')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)
        save_sig.add(ss47)
        save_sig.add(ss2)
        save_sig.add(ss47)

    assert os.path.isdir(outloc)

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert ss47 in saved
    assert len(saved) == 4


def test_load_empty_zipfile(runtmp):
    outloc = runtmp.output('empty.zip')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        pass

    sigiter = sourmash.load_file_as_signatures(outloc)
    assert list(sigiter) == []


def test_load_many_sigs_empty_file(runtmp):
    # make sure load_many_signatures behaves properly on empty file
    outloc = runtmp.output("empty.sig")
    with open(outloc, "wt") as fp:
        pass

    progress = sourmash_args.SignatureLoadingProgress()

    with contextlib.redirect_stderr(io.StringIO()) as errfp:
        with pytest.raises(SystemExit) as exc:
            for ss, sigloc in sourmash_args.load_many_signatures([outloc],
                                                                 progress):
                pass

    err = errfp.getvalue()
    print(err)
    assert f"ERROR: Error while reading signatures from '{outloc}'." in err
    assert "(continuing)" not in err


def test_load_many_sigs_empty_file_force(runtmp):
    # make sure load_many_signatures behaves properly on empty file w/force
    outloc = runtmp.output("empty.sig")
    with open(outloc, "wt") as fp:
        pass

    progress = sourmash_args.SignatureLoadingProgress()

    with contextlib.redirect_stderr(io.StringIO()) as errfp:
        for ss, sigloc in sourmash_args.load_many_signatures([outloc],
                                                             progress,
                                                             force=True):
            pass

    err = errfp.getvalue()
    print(err)
    assert f"ERROR: Error while reading signatures from '{outloc}'." in err
    assert "(continuing)" in err
