"""
Tests for functions in sourmash_args module.
"""
import sys
import os
import pytest
import gzip
import zipfile
import io
import contextlib
import csv
import argparse
import shutil

import sourmash_tst_utils as utils
import sourmash
from sourmash import sourmash_args, manifest
from sourmash.index import LinearIndex
from sourmash.cli.utils import add_ksize_arg


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


def test_save_signatures_to_location_1_zip_bad(runtmp):
    # try saving to bad sigfile.zip
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    outloc = runtmp.output('foo.zip')

    # create bad zip:
    with open(outloc, 'wt') as fp:
        pass

    # now check for error
    with pytest.raises(ValueError) as exc:
        with sourmash_args.SaveSignaturesToLocation(outloc) as save_sig:
            pass

    assert 'cannot be opened as a zip file' in str(exc)


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
        ss2 = ss2.to_mutable()
        ss2.name = 'different name for ss2'
        save_sig.add(ss2)

        ss47 = ss47.to_mutable()
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
        ss2copy = ss2.to_mutable()
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
        mf_name = "SOURMASH-MANIFEST.csv"

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


def test_get_manifest_1():
    # basic get_manifest retrieves a manifest
    sig47 = utils.get_test_data('47.fa.sig')
    idx = sourmash.load_file_as_index(sig47)

    manifest = sourmash_args.get_manifest(idx)
    assert len(manifest) == 1


def test_get_manifest_2_cannot_build():
    # test what happens when get_manifest cannot build manifest
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47)

    idx = LinearIndex([ss47])

    with pytest.raises(SystemExit) as exc:
        m = sourmash_args.get_manifest(idx)


def test_get_manifest_2_cannot_buildno_require():
    # test what happens when get_manifest cannot build manifest
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47)

    idx = LinearIndex([ss47])

    m = sourmash_args.get_manifest(idx, require=False)

    assert m is None


def test_get_manifest_3_build():
    # check that manifest is building
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47)

    class FakeIndex(LinearIndex):
        was_called = 0
        def _signatures_with_internal(self):
            self.was_called = 1
            return [(ss47, "fakeiloc")]

    idx = FakeIndex([sig47])

    assert not idx.was_called
    m = sourmash_args.get_manifest(idx)
    assert idx.was_called

    print(m)
    assert len(m) == 1
    assert m.rows[0]['internal_location'] == "fakeiloc"


def test_get_manifest_3_build_2():
    # check that manifest is building, but only when asked
    sig47 = utils.get_test_data('47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47)

    class FakeIndex(LinearIndex):
        manifest = None
        was_called = 0

        def _signatures_with_internal(self):
            self.was_called = 1
            return [(ss47, "fakeiloc")]

    idx = FakeIndex([sig47])

    assert not idx.was_called
    m = sourmash_args.get_manifest(idx)
    assert idx.was_called

    # now set and ask again, should not be called
    idx.manifest = m
    idx.was_called = 0

    m2 = sourmash_args.get_manifest(idx)
    assert not idx.was_called
    assert m == m2

    # now, force rebuild
    m3 = sourmash_args.get_manifest(idx, rebuild=True)
    assert idx.was_called
    assert m == m3


class FakeArgs(object):
    picklist = None
    include_db_pattern = None
    exclude_db_pattern = None


def test_pattern_0():
    # test neither --include nor --exclude
    args = FakeArgs()
    args.picklist = None
    args.include_db_pattern = None
    args.exclude_db_pattern = None

    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    assert pattern_search is None


def test_pattern_1():
    # test just --include-pattern handling
    args = FakeArgs()
    args.picklist = None
    args.include_db_pattern = 'foo'
    args.exclude_db_pattern = None

    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    assert pattern_search(['foo', 'bar', 'baz'])
    assert not pattern_search(['bar', 'bif'])


def test_pattern_2():
    # test just --exclude-pattern handling
    args = FakeArgs()
    args.picklist = None
    args.exclude_db_pattern = 'foo'
    args.include_db_pattern = None

    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)
    assert not pattern_search(['foo', 'bar', 'baz'])
    assert pattern_search(['bar', 'baz', 'bif'])


def test_pattern_3():
    # test with --picklist and --exclude: should fail
    args = FakeArgs()
    args.picklist = True
    args.exclude_db_pattern = 'foo'
    args.include_db_pattern = None

    with pytest.raises(SystemExit):
        pattern_search = sourmash_args.load_include_exclude_db_patterns(args)


def test_pattern_4():
    # test with --picklist and --include: should fail
    args = FakeArgs()
    args.picklist = True
    args.include_db_pattern = 'foo'
    args.exclude_db_pattern = None

    with pytest.raises(SystemExit):
        pattern_search = sourmash_args.load_include_exclude_db_patterns(args)


def test_pattern_5():
    # test with --include and --exclude: should fail
    args = FakeArgs()
    args.picklist = None
    args.exclude_db_pattern = 'foo'
    args.include_db_pattern = 'bar'

    with pytest.raises(SystemExit):
        pattern_search = sourmash_args.load_include_exclude_db_patterns(args)


def test_fileinput_csv_1_plain():
    # test basic CSV input

    testfile = utils.get_test_data('tax/test.taxonomy.csv')

    with sourmash_args.FileInputCSV(testfile) as r:
        rows = list(r)
        assert len(rows) == 6


def test_fileinput_csv_1_no_such_file(runtmp):
    # test fail to load file

    noexistfile = runtmp.output('does-not-exist.csv')

    with pytest.raises(FileNotFoundError):
        with sourmash_args.FileInputCSV(noexistfile) as r:
            pass


def test_fileinput_csv_2_gz(runtmp):
    # test basic CSV input from gz file

    testfile = utils.get_test_data('tax/test.taxonomy.csv')
    gzfile = runtmp.output('test.csv.gz')

    with gzip.open(gzfile, 'wt') as outfp:
        with open(testfile, 'rt', newline='') as infp:
            outfp.write(infp.read())

    with sourmash_args.FileInputCSV(gzfile) as r:
        rows = list(r)
        assert len(rows) == 6


def test_fileinput_csv_2_gz_not_csv(runtmp):
    # test basic CSV input from gz file that's not CSV - works

    gzfile = runtmp.output('test.csv.gz')

    with gzip.open(gzfile, 'wt') as outfp:
        outfp.write("hello world!")

    with sourmash_args.FileInputCSV(gzfile) as r:
        assert r.fieldnames == ['hello world!']


def test_fileinput_csv_2_gz_bad_version_header(runtmp):
    # test basic CSV input from gz file with bad version header
    # currently this works; not clear to me how it should fail :grin:

    gzfile = runtmp.output('test.csv.gz')

    with gzip.open(gzfile, 'wt') as outfp:
        outfp.write("# excelsior\nhello world!")

    with sourmash_args.FileInputCSV(gzfile) as r:
        assert r.fieldnames == ['hello world!']
        print(r.version_info)
        assert r.version_info == ['excelsior']


def test_fileinput_csv_2_zip(runtmp):
    # test CSV input from zip file, with component filename

    testfile = utils.get_test_data('tax/test.taxonomy.csv')
    zf_file = runtmp.output('test.zip')

    with zipfile.ZipFile(zf_file, 'w') as outzip:
        with open(testfile, 'rb') as infp:
            with outzip.open('XYZ.csv', 'w') as outfp:
                outfp.write(infp.read())

    with sourmash_args.FileInputCSV(zf_file, default_csv_name='XYZ.csv') as r:
        rows = list(r)
        assert len(rows) == 6
        print(rows)


def test_fileinput_csv_3_load_manifest():
    # test loading a manifest from a zipfile collection, using
    # FileInputCSV.
    testfile = utils.get_test_data('prot/all.zip')

    with sourmash_args.FileInputCSV(testfile, default_csv_name='SOURMASH-MANIFEST.csv') as r:

        rows = list(r)
        assert len(rows) == 8

        assert r.version_info == ['SOURMASH-MANIFEST-VERSION', '1.0']


def test_fileinput_csv_3_load_manifest_no_default():
    # test loading a manifest from a zipfile collection, using
    # FileInputCSV, but with no default_csv_name - should fail
    testfile = utils.get_test_data('prot/all.zip')

    with pytest.raises(csv.Error):
        with sourmash_args.FileInputCSV(testfile) as r:
            print(r.fieldnames)


def test_fileinput_csv_3_load_manifest_zipfile_obj():
    # test loading a manifest from an open zipfile obj, using
    # FileInputCSV.
    testfile = utils.get_test_data('prot/all.zip')

    with zipfile.ZipFile(testfile, "r") as zf:
        with sourmash_args.FileInputCSV(testfile,
                                     default_csv_name='SOURMASH-MANIFEST.csv',
                                     zipfile_obj=zf) as r:
            rows = list(r)
            assert len(rows) == 8

            assert r.version_info == ['SOURMASH-MANIFEST-VERSION', '1.0']


def test_fileinput_csv_3_load_manifest_zipfile_obj_no_defualt():
    # test loading a manifest from an open zipfile obj, using
    # FileInputCSV, but with no default csv name => should fail.
    testfile = utils.get_test_data('prot/all.zip')

    with zipfile.ZipFile(testfile, "r") as zf:
        with pytest.raises(ValueError):
            with sourmash_args.FileInputCSV(testfile,
                                            zipfile_obj=zf) as r:
                pass


def test_fileoutput_csv_1(runtmp):
    # test basic behavior
    outfile = runtmp.output('xxx.csv')

    with sourmash_args.FileOutputCSV(outfile) as fp:
        w = csv.writer(fp)
        w.writerow(['a', 'b', 'c'])
        w.writerow(['x', 'y', 'z'])

    with open(outfile, newline="") as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1
        row = rows[0]
        assert row['a'] == 'x'
        assert row['b'] == 'y'
        assert row['c'] == 'z'


def test_fileoutput_csv_1_gz(runtmp):
    # test basic behavior => gz
    outfile = runtmp.output('xxx.csv.gz')

    with sourmash_args.FileOutputCSV(outfile) as fp:
        w = csv.writer(fp)
        w.writerow(['a', 'b', 'c'])
        w.writerow(['x', 'y', 'z'])

    with gzip.open(outfile, 'rt') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        assert len(rows) == 1
        row = rows[0]
        assert row['a'] == 'x'
        assert row['b'] == 'y'
        assert row['c'] == 'z'


def test_fileoutput_csv_2_stdout():
    # test '-' and 'None' go to sys.stdout

    with sourmash_args.FileOutputCSV('-') as fp:
        assert fp == sys.stdout

    with sourmash_args.FileOutputCSV(None) as fp:
        assert fp == sys.stdout


def test_add_ksize_arg_no_default():
    # test behavior of cli.utils.add_ksize_arg
    p = argparse.ArgumentParser()
    add_ksize_arg(p)
    args = p.parse_args()
    assert args.ksize == None


def test_add_ksize_arg_no_default_specify():
    # test behavior of cli.utils.add_ksize_arg
    p = argparse.ArgumentParser()
    add_ksize_arg(p)
    args = p.parse_args(['-k', '21'])
    assert args.ksize == 21


def test_add_ksize_arg_default_31():
    # test behavior of cli.utils.add_ksize_arg
    p = argparse.ArgumentParser()
    add_ksize_arg(p, default=31)
    args = p.parse_args()
    assert args.ksize == 31


def test_add_ksize_arg_default_31_specify():
    # test behavior of cli.utils.add_ksize_arg
    p = argparse.ArgumentParser()
    add_ksize_arg(p, default=31)
    args = p.parse_args(['-k', '21'])
    assert args.ksize == 21


def test_bug_2370(runtmp):
    # bug - manifest loading code does not catch gzip.BadGzipFile
    sigfile = utils.get_test_data('63.fa.sig')

    # copy sigfile over to a .gz file without compressing it -
    shutil.copyfile(sigfile, runtmp.output('not_really_gzipped.gz'))

    # try running sourmash_args.load_file_as_index
    #runtmp.sourmash('sig', 'describe', runtmp.output('not_really_gzipped.gz'))
    sourmash_args.load_file_as_index(runtmp.output('not_really_gzipped.gz'))
