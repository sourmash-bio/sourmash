"""
Tests for 'sourmash sig collect'
"""
import pytest

import sourmash
from sourmash.manifest import BaseCollectionManifest

import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_sig_collect_1_zipfile(runtmp, manifest_db_format):
    # collect a manifest from a .zip file
    protzip = utils.get_test_data('prot/protein.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', protzip, '-o', f'mf.{ext}',
                    '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_collect_2_exists_fail(runtmp, manifest_db_format):
    # collect a manifest from two .zip files
    protzip = utils.get_test_data('prot/protein.zip')
    allzip = utils.get_test_data('prot/protein.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', protzip, '-o', f'mf.{ext}',
                    '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    # now run with same filename - should fail
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                        '-F', manifest_db_format)


def test_sig_collect_2_exists_force(runtmp, manifest_db_format):
    # collect a manifest from two .zip files
    protzip = utils.get_test_data('prot/protein.zip')
    allzip = utils.get_test_data('prot/all.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', protzip, '-o', f'mf.{ext}',
                    '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    # now run with same filename - should fail
    runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                    '-F', manifest_db_format, '--force-overwrite')


def test_sig_collect_3_multiple(runtmp, manifest_db_format):
    # collect a manifest from two .zip files
    protzip = utils.get_test_data('prot/protein.zip')
    hpzip = utils.get_test_data('prot/hp.zip')
    dayzip = utils.get_test_data('prot/dayhoff.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', protzip, hpzip, dayzip,
                    '-o', f'mf.{ext}', '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 6
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list
    assert 'ea2a1ad233c2908529d124a330bcb672' in md5_list
    assert 'bb0e6d90df01b7bd5d0956a5f9e3ed12' in md5_list
    assert 'fbca5e5211e4d58427997fd5c8343e9a' in md5_list
    assert '1cbd888bf910f83ad8f1715509183223' in md5_list

    locations = set([ row['internal_location'] for row in manifest.rows ])
    assert protzip in locations
    assert hpzip in locations
    assert dayzip in locations
    assert len(locations) == 3, locations


def test_sig_collect_4_multiple_from_sig(runtmp, manifest_db_format):
    # collect a manifest from sig files
    sig43 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', sig43, sig63,
                    '-o', f'mf.{ext}', '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '09a08691ce52952152f0e866a59f6261' in md5_list
    assert '38729c6374925585db28916b82a6f513' in md5_list

    locations = set([ row['internal_location'] for row in manifest.rows ])
    assert sig43 in locations
    assert sig63 in locations
    assert len(locations) == 2, locations


def test_sig_collect_5_no_manifest_sbt_fail(runtmp, manifest_db_format):
    # collect a manifest from files that don't have one
    sbt_zip = utils.get_test_data('v6.sbt.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'collect', sbt_zip,
                        '-o', f'mf.{ext}', '-F', manifest_db_format)


def test_sig_collect_5_no_manifest_sbt_succeed(runtmp, manifest_db_format):
    # generate a manifest from files that don't have one when --no-require
    sbt_zip = utils.get_test_data('v6.sbt.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', sbt_zip, '--no-require-manifest',
                    '-o', f'mf.{ext}', '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 7
    locations = set([ row['internal_location'] for row in manifest.rows ])
    assert len(locations) == 1, locations
    assert sbt_zip in locations
