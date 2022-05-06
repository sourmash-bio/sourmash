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
    runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                    '-F', manifest_db_format, '--force-overwrite')
