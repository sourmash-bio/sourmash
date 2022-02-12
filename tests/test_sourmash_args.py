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
from sourmash import sourmash_args, manifest


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
        # exact duplicates will not be saved, otherwise.
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


def test_save_signatures_to_location_2_zip_add(runtmp):
    # create sigfile.zip; then, add a new signature.
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    # add only ss2
    outloc = runtmp.output('foo.zip')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)

    # can we open as a .zip file?
    with zipfile.ZipFile(outloc, "r") as zf:
        assert list(zf.infolist())

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert len(saved) == 1

    # now, re-open and add ss47.
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss47)

    # updated file should contain both.
    saved = list(sourmash.load_file_as_signatures(outloc))
    print(saved)
    assert ss47 in saved
    assert ss2 in saved


def test_save_signatures_to_location_2_zip_add_dup(runtmp):
    # create sigfile.zip; then, add a new signature, plus a ~duplicate.
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    # add only ss2
    outloc = runtmp.output('foo.zip')
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss2)

    # can we open as a .zip file?
    with zipfile.ZipFile(outloc, "r") as zf:
        assert list(zf.infolist())

    saved = list(sourmash.load_file_as_signatures(outloc))
    assert ss2 in saved
    assert len(saved) == 1

    # now, re-open and add ss47, plus a slightly renamed ss2.
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss47)

        # add ss2; here we have to change the names so the sig content is
        # different exact duplicates will not be saved, otherwise.
        import copy
        ss2copy = copy.copy(ss2)
        ss2copy.name = 'different name for ss2'
        save_sig.add(ss2copy)

    # updated file should contain all three.
    saved = list(sourmash.load_file_as_signatures(outloc))
    print(len(saved), saved)
    assert ss47 in saved
    assert ss2 in saved
    assert ss2copy in saved


def test_save_signatures_to_location_3_zip_add_fail(runtmp):
    # create sigfile.zip using zipfile, then try to add to it (& fail)
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    # add only ss2, using zipfile API
    outloc = runtmp.output('foo.zip')
    with zipfile.ZipFile(outloc, 'x') as zf:
        with zf.open('xyz.sig', 'w') as fp:
            sourmash.save_signatures([ss2], fp=fp, compression=1)

    # verify it can be loaded, yada yada
    saved = list(sourmash.load_file_as_signatures(outloc))
    print(len(saved), saved)
    assert ss2 in saved

    # now, try to open existing file with SaveSignaturesToLocation...
    with pytest.raises(ValueError) as exc:
        with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
            pass

    assert 'Cannot add to existing zipfile' in str(exc)


def test_save_signatures_to_location_3_zip_add_with_manifest(runtmp):
    # create sigfile.zip using zipfile, then try to add to it (& fail)
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    # add only ss2, using zipfile API; add manifest manually.
    outloc = runtmp.output('foo.zip')
    with zipfile.ZipFile(outloc, 'x') as zf:
        with zf.open('xyz.sig', 'w') as fp:
            sourmash.save_signatures([ss2], fp=fp, compression=1)

        # make a manifest row...
        row = manifest.CollectionManifest.make_manifest_row(ss2, 'xyz.sig',
                                                   include_signature=False)

        # construct & save manifest
        mf = manifest.CollectionManifest([row])
        mf_name = f"SOURMASH-MANIFEST.csv"

        manifest_fp = io.StringIO()
        mf.write_to_csv(manifest_fp, write_header=True)
        manifest_data = manifest_fp.getvalue().encode("utf-8")

        with zf.open(mf_name, 'w') as fp:
            fp.write(manifest_data)

        # fini! made our artisanal hand-crafted zipfile. Now...

    # verify it can be loaded, yada yada
    saved = list(sourmash.load_file_as_signatures(outloc))
    print(len(saved), saved)
    assert ss2 in saved

    # now, try to open existing file with SaveSignaturesToLocation...
    # ...should succeed!
    with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
        print(save_sig)
        save_sig.add(ss47)

    # updated file should contain both
    saved = list(sourmash.load_file_as_signatures(outloc))
    print(len(saved), saved)
    assert ss47 in saved
    assert ss2 in saved


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
