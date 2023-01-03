"""
Tests for 'sourmash sig collect'
"""
import pytest
import shutil
import os.path
import gzip

import sourmash
from sourmash.manifest import BaseCollectionManifest

import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_sig_collect_0_nothing(runtmp, manifest_db_format):
    # run with just output
    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'
    if manifest_db_format != 'sql': return

    runtmp.sourmash('sig', 'collect', '-o', f'mf.{ext}',
                    '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 0


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


def test_sig_collect_1_zipfile_csv_gz(runtmp):
    # collect a manifest from a .zip file, save to csv.gz
    protzip = utils.get_test_data('prot/protein.zip')

    runtmp.sourmash('sig', 'collect', protzip, '-o', 'mf.csv.gz',
                    '-F', 'csv')

    manifest_fn = runtmp.output('mf.csv.gz')

    # gzip, yes?
    print('XXX', manifest_fn)
    with gzip.open(manifest_fn, 'rt', newline='') as fp:
        fp.read()

    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list


def test_sig_collect_1_zipfile_csv_gz_roundtrip(runtmp):
    # collect a manifest from a .zip file, save to csv.gz; then load again
    protzip = utils.get_test_data('prot/protein.zip')

    runtmp.sourmash('sig', 'collect', protzip, '-o', 'mf.csv.gz',
                    '-F', 'csv')

    manifest_fn = runtmp.output('mf.csv.gz')

    # gzip, yes?
    print('XXX', manifest_fn)
    with gzip.open(manifest_fn, 'rt', newline='') as fp:
        fp.read()

    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    # can we read a csv.gz?
    runtmp.sourmash('sig', 'collect', 'mf.csv.gz', '-o', 'mf2.csv',
                    '-F', 'csv')

    manifest_fn2 = runtmp.output('mf2.csv')
    manifest2 = BaseCollectionManifest.load_from_filename(manifest_fn2)

    assert len(manifest2) == 2
    md5_list = [ row['md5'] for row in manifest2.rows ]
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


def test_sig_collect_2_exists_merge(runtmp, manifest_db_format):
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

    # now run with same filename - should merge
    runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                    '-F', manifest_db_format, '--merge')

    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)
    assert len(manifest) == 10


def test_sig_collect_2_exists_sql_merge_csv(runtmp, manifest_db_format):
    # try to merge csv into sql
    protzip = utils.get_test_data('prot/protein.zip')
    allzip = utils.get_test_data('prot/all.zip')

    ext = 'sqlmf'

    # save as sql...
    runtmp.sourmash('sig', 'collect', protzip, '-o', f'mf.{ext}',
                    '-F', 'sql')

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                        '-F', 'csv', '--merge')

    assert "ERROR loading" in runtmp.last_result.err


def test_sig_collect_2_exists_csv_merge_sql(runtmp):
    # try to merge sql into csv
    protzip = utils.get_test_data('prot/protein.zip')
    allzip = utils.get_test_data('prot/all.zip')

    ext = 'csv'

    # save as csv...
    runtmp.sourmash('sig', 'collect', protzip, '-o', f'mf.{ext}',
                    '-F', 'csv')

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '16869d2c8a1d29d1c8e56f5c561e585e' in md5_list
    assert '120d311cc785cc9d0df9dc0646b2b857' in md5_list

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                        '-F', 'sql', '--merge')

    assert "ERROR loading" in runtmp.last_result.err


def test_sig_collect_2_no_exists_merge(runtmp, manifest_db_format):
    # test 'merge' when args.output doesn't already exist => warning
    protzip = utils.get_test_data('prot/protein.zip')
    allzip = utils.get_test_data('prot/all.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'
    manifest_fn = runtmp.output(f'mf.{ext}')

    # run with --merge but no previous:
    runtmp.sourmash('sig', 'collect', allzip, '-o', manifest_fn,
                    '-F', manifest_db_format, '--merge')

    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)
    assert len(manifest) == 8

    err = runtmp.last_result.err
    print(err)
    assert "WARNING: --merge-previous specified, but output file" in err


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


def test_sig_collect_3_multiple_use_fromfile(runtmp, manifest_db_format):
    # collect a manifest from two .zip files using --from-file
    protzip = utils.get_test_data('prot/protein.zip')
    hpzip = utils.get_test_data('prot/hp.zip')
    dayzip = utils.get_test_data('prot/dayhoff.zip')

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    fromfile = runtmp.output('fromfile.txt')
    with open(fromfile, 'wt') as fp:
        print(protzip, file=fp)
        print(hpzip, file=fp)
        print(dayzip, file=fp)

    runtmp.sourmash('sig', 'collect', '--from-file', 'fromfile.txt',
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


def test_sig_collect_4_multiple_from_sig_abspath(runtmp, manifest_db_format):
    # collect a manifest from sig files, forcing abspath
    sig43 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    shutil.copyfile(sig43, runtmp.output('47.fa.sig'))
    shutil.copyfile(sig63, runtmp.output('63.fa.sig'))

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', '47.fa.sig', '63.fa.sig', '--abspath',
                    '-o', f'mf.{ext}', '-F', manifest_db_format)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '09a08691ce52952152f0e866a59f6261' in md5_list
    assert '38729c6374925585db28916b82a6f513' in md5_list

    locations = set([ row['internal_location'] for row in manifest.rows ])
    print(locations)
    assert len(locations) == 2, locations

    for xx in locations:
        assert xx.startswith('/')


def test_sig_collect_4_multiple_no_abspath(runtmp, manifest_db_format):
    # collect a manifest from sig files, no abspath
    sig43 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    # copy files to tmp, where they will not have full paths
    shutil.copyfile(sig43, runtmp.output('47.fa.sig'))
    shutil.copyfile(sig63, runtmp.output('63.fa.sig'))

    ext = 'sqlmf' if manifest_db_format == 'sql' else 'csv'

    runtmp.sourmash('sig', 'collect', '47.fa.sig', '63.fa.sig',
                    '-o', f'mf.{ext}', '-F', manifest_db_format)

    manifest_fn = runtmp.output(f'mf.{ext}')
    manifest = BaseCollectionManifest.load_from_filename(manifest_fn)

    assert len(manifest) == 2
    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '09a08691ce52952152f0e866a59f6261' in md5_list
    assert '38729c6374925585db28916b82a6f513' in md5_list

    locations = set([ row['internal_location'] for row in manifest.rows ])
    print(locations)
    assert len(locations) == 2, locations
    assert '47.fa.sig' in locations
    assert '63.fa.sig' in locations


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
