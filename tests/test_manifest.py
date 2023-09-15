"""
Tests for manifest code in databases, etc.
"""
import pytest
from io import StringIO

import sourmash
from sourmash import index, sourmash_args

import sourmash_tst_utils as utils


def test_generate_manifest():
    # test basic manifest-generating functionality.
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)

    assert len(manifest) == len(rows)
    assert len(manifest) == 2

    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    for sig in siglist:
        assert sig in manifest


def test_manifest_operations():
    # test basic manifest operations - +=
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)
    manifest2 = index.CollectionManifest(rows)
    manifest += manifest2

    assert len(manifest) == 2*len(rows)
    assert len(manifest) == 4

    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_manifest_operations_fail():
    # should not be able to add a manifest to itself - not only makes
    # no sense, but it means you're modifying a generator in place, sometimes.
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)
    with pytest.raises(Exception):
        manifest += manifest


def test_manifest_to_picklist():
    # test manifest/picklist interaction basics
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)
    picklist = manifest.to_picklist()
    assert len(picklist.pickset) == len(manifest)

    new_manifest = manifest.select_to_manifest(picklist=picklist)
    assert len(new_manifest) == len(manifest)


def test_manifest_compare():
    # test saving and loading manifests
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)
    manifest = loader.manifest

    # equal
    rows = list(manifest.rows)

    equal_mf = index.CollectionManifest(rows)
    assert equal_mf == manifest

    # not equal / shorter
    rows = list(manifest.rows)
    rows = rows[:-1]

    short_mf = index.CollectionManifest(rows)
    assert short_mf != manifest

    # not equal / diff values
    rows = list(manifest.rows)
    rows[0] = dict(rows[0])
    rows[0]['internal_location'] += '.foo'

    short_mf = index.CollectionManifest(rows)
    assert short_mf != manifest


def test_save_load_manifest():
    # test saving and loading manifests
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)

    # now, on to CSV
    fp = StringIO()
    manifest.write_csv_header(fp)
    manifest.write_to_csv(fp)

    rfp = StringIO(fp.getvalue())
    manifest2 = index.CollectionManifest.load_from_csv(rfp)

    assert len(manifest) == len(manifest2)

    pick1 = manifest.to_picklist()
    pick2 = manifest2.to_picklist()

    # manifest 1 in manifest2?
    for row in manifest.rows:
        assert pick2.matches_manifest_row(row)

    # manifest 2 in manifest?
    for row in manifest2.rows:
        assert pick1.matches_manifest_row(row)

    # equal?
    assert manifest == manifest2

    # not equal / shorter
    rows = list(manifest.rows)
    rows = rows[1:]

    short_mf = index.CollectionManifest(rows)
    assert short_mf != manifest

    # not equal / diff values
    rows = list(manifest.rows)
    rows[0] = dict(rows[0])
    rows[0]['internal_location'] += '.foo'

    short_mf = index.CollectionManifest(rows)
    assert short_mf != manifest


def test_manifest_to_picklist_bug(runtmp):
    # this tests a fun combination of things that led to a bug.
    # tl;dr we only want to iterate once across a generator...
    # ref #2762
    c = runtmp
    all_zip = utils.get_test_data('prot/all.zip')

    idx = sourmash_args.load_file_as_index(all_zip)
    assert len(idx) == 8

    manifest = sourmash_args.get_manifest(idx)
    assert len(manifest) == 8

    def filter_fn(row):
        # match?
        keep = False
        if "09a0869" in row['md5']:
            keep = True

        return keep

    sub_manifest = manifest.filter_rows(filter_fn)
    sub_picklist = sub_manifest.to_picklist()
    idx = idx.select(picklist=sub_picklist)

    assert len(idx) == 1
    print(idx)

    x = list(idx.signatures())
    assert len(x)


def test_generate_manifest_iterate_once():
    # we should only iterate across manifest rows once
    protzip = utils.get_test_data('prot/protein.zip')

    loader = sourmash.load_file_as_index(protzip)

    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        siglist.append(sig)

    # build generator function => will not allow iteration twice
    def genfn():
        for (sig, loc) in loader._signatures_with_internal():
            row = index.CollectionManifest.make_manifest_row(sig, loc)
            yield row

    manifest = index.CollectionManifest(genfn())

    assert len(manifest) == 2
    assert len(manifest._md5_set) == 2

    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    for sig in siglist:
        assert sig in manifest
