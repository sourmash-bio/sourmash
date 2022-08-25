"""
Tests for the 'sourmash signature fileinfo' command line.
"""
import shutil
import os
import glob
import json

import pytest

import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed

## command line tests


def test_fileinfo_1_sig(runtmp):
    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')

    shutil.copyfile(sig47, runtmp.output('sig47.sig'))
    runtmp.run_sourmash('sig', 'fileinfo', 'sig47.sig')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: sig47.sig
is database? no
has manifest? yes
num signatures: 1
total hashes: 5177
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000             5177
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_1_sig_summarize(runtmp):
    # get basic info on a signature with 'summarize' as alias for fileinfo
    sig47 = utils.get_test_data('47.fa.sig')

    shutil.copyfile(sig47, runtmp.output('sig47.sig'))
    runtmp.run_sourmash('sig', 'summarize', 'sig47.sig')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: sig47.sig
is database? no
has manifest? yes
num signatures: 1
total hashes: 5177
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000             5177
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_1_sig_abund(runtmp):
    # get basic info on a signature with abundance
    sig47 = utils.get_test_data('track_abund/47.fa.sig')

    shutil.copyfile(sig47, runtmp.output('sig47.sig'))
    runtmp.run_sourmash('sig', 'fileinfo', 'sig47.sig')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: sig47.sig
is database? no
has manifest? yes
num signatures: 1
5177 total hashes
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000, abund      5177 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_2_lca(runtmp):
    # get basic info on an LCA database
    prot = utils.get_test_data('prot/protein.lca.json.gz')

    shutil.copyfile(prot, runtmp.output('protein.lca.json.gz'))
    runtmp.run_sourmash('sig', 'fileinfo', 'protein.lca.json.gz')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: LCA_Database
location: protein.lca.json.gz
is database? yes
has manifest? no
num signatures: 2
total hashes: 8214
summary of sketches:
   2 sketches with protein, k=19, scaled=100          8214
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_3_sbt_zip(runtmp):
    # test on an SBT.zip
    prot = utils.get_test_data('prot/protein.sbt.zip')

    shutil.copyfile(prot, runtmp.output('protein.sbt.zip'))
    runtmp.run_sourmash('sig', 'fileinfo', 'protein.sbt.zip')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: SBT
location: protein.sbt.zip
is database? yes
has manifest? yes
num signatures: 2
total hashes: 8214
summary of sketches:
   2 sketches with protein, k=19, scaled=100          8214 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()


def test_fileinfo_4_zip(runtmp):
    # test on a ZipFileLinearIndex
    prot = utils.get_test_data('prot/all.zip')

    shutil.copyfile(prot, runtmp.output('all.zip'))
    runtmp.run_sourmash('sig', 'fileinfo', 'all.zip')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    # 'location' will be fully resolved, ignore it for now
    expected_output = """\
path filetype: ZipFileLinearIndex
is database? yes
has manifest? yes
num signatures: 8
total hashes: 31758
summary of sketches:
   2 sketches with dayhoff, k=19, scaled=100          7945 total hashes
   2 sketches with hp, k=19, scaled=100               5184 total hashes
   2 sketches with protein, k=19, scaled=100          8214 total hashes
   2 sketches with DNA, k=31, scaled=1000             10415 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_4_zip_json_out(runtmp):
    # check --json-out
    prot = utils.get_test_data('prot/all.zip')

    shutil.copyfile(prot, runtmp.output('all.zip'))
    runtmp.run_sourmash('sig', 'fileinfo', 'all.zip', '--json-out')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    # should succeed as loading as JSON, with correct info
    vals = json.loads(out)

    assert vals['has_manifest']
    assert vals['is_database']
    assert vals['num_sketches'] == 8
    assert vals['path_filetype'] == 'ZipFileLinearIndex'
    assert vals['total_hashes'] == 31758

    d1 = {'ksize': 19, 'moltype': 'dayhoff', 'scaled': 100, 'num': 0, 'abund': False, 'count': 2, 'n_hashes': 7945}
    d2 = {'ksize': 19, 'moltype': 'hp', 'scaled': 100, 'num': 0, 'abund': False, 'count': 2, 'n_hashes': 5184}
    d3 = {'ksize': 19, 'moltype': 'protein', 'scaled': 100, 'num': 0, 'abund': False, 'count': 2, 'n_hashes': 8214}
    d4 = {'ksize': 31, 'moltype': 'DNA', 'scaled': 1000, 'num': 0, 'abund': False, 'count': 2, 'n_hashes': 10415}

    assert d1 in vals['sketch_info']
    assert d2 in vals['sketch_info']
    assert d3 in vals['sketch_info']
    assert d4 in vals['sketch_info']
    assert len(vals['sketch_info']) == 4


def test_fileinfo_4_zip_rebuild(runtmp):
    # test --rebuild
    prot = utils.get_test_data('prot/all.zip')

    shutil.copyfile(prot, runtmp.output('all.zip'))
    runtmp.run_sourmash('sig', 'fileinfo', 'all.zip', '--rebuild')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    # 'location' will be fully resolved, ignore it for now
    # CTB: note we're missing one of the 8 in the rebuilt, dna-sig.noext,
    # because it is not automatically included unless you load the zipfile
    # with traverse. This is intentional.
    expected_output = """\
path filetype: ZipFileLinearIndex
is database? yes
has manifest? yes
num signatures: 8
total hashes: 26581
summary of sketches:
   2 sketches with dayhoff, k=19, scaled=100          7945 total hashes
   2 sketches with hp, k=19, scaled=100               5184 total hashes
   2 sketches with protein, k=19, scaled=100          8214 total hashes
   1 sketches with DNA, k=31, scaled=1000             5238 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_5_dir(runtmp):
    # test on a directory
    sig47 = utils.get_test_data('47.fa.sig')

    os.mkdir(runtmp.output('subdir'))

    shutil.copyfile(sig47, runtmp.output('subdir/sig47.sig'))
    runtmp.run_sourmash('sig', 'fileinfo', 'subdir/')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: subdir/
is database? no
has manifest? yes
num signatures: 1
total hashes: 5177
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000             5177 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_6_pathlist(runtmp):
    # test on a pathlist
    sig47 = utils.get_test_data('47.fa.sig')
    shutil.copyfile(sig47, runtmp.output("47.fa.sig"))

    with open(runtmp.output('pathlist.txt'), 'wt') as fp:
        fp.write("47.fa.sig\n")

    runtmp.run_sourmash('sig', 'fileinfo', 'pathlist.txt')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: pathlist.txt
is database? no
has manifest? yes
num signatures: 1
total hashes: 5177
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000             5177 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@pytest.mark.parametrize("db", ['v6.sbt.json', 'v5.sbt.json', 'v4.sbt.json',
                                'v3.sbt.json', 'v2.sbt.json', 'v1.sbt.json'])
def test_fileinfo_7_sbt_json(runtmp, db):
    # test on multiple versions of SBT JSON files
    dbfile = utils.get_test_data(db)

    runtmp.run_sourmash('sig', 'fileinfo', dbfile)

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    expected_output = f"""\
path filetype: SBT
location: {dbfile}
is database? yes
has manifest? no
num signatures: 7
total hashes: 3500
summary of sketches:
   7 sketches with DNA, k=31, num=500                 3500 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()


def test_sig_fileinfo_stdin(runtmp):
    # test on stdin
    sig = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    with open(sig, 'rt') as fp:
        data = fp.read()

    runtmp.run_sourmash('sig', 'fileinfo', '-', stdin_data=data)

    out = runtmp.last_result.out
    print(out)

    expected_output = """\
path filetype: MultiIndex
location: -
is database? no
has manifest? yes
num signatures: 1
total hashes: 3409
summary of sketches:
   1 sketches with protein, k=19, scaled=100          3409 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()


def test_sig_fileinfo_does_not_exist(runtmp):
    # test on file that does not exist
    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'fileinfo', 'does-not-exist')

    assert "Cannot open 'does-not-exist' as a sourmash signature collection" in runtmp.last_result.err


def test_sig_fileinfo_8_manifest_works(runtmp):
    # test on a manifest with relative paths, in proper location
    mf = utils.get_test_data('scaled/mf.csv')
    runtmp.sourmash('sig', 'fileinfo', mf)

    out = runtmp.last_result.out
    print(out)

    assert '15 sketches with DNA, k=31, scaled=10000           717 total hashes' in out
    assert 'num signatures: 15' in out
    assert 'has manifest? yes' in out
    assert 'is database? yes' in out
    assert 'path filetype: StandaloneManifestIndex' in out


def test_sig_fileinfo_8_manifest_works_when_moved(runtmp):
    # test on a manifest with relative paths, when in wrong place
    # note: this works, unlike 'describe', because all the necessary info
    # for 'fileinfo' is in the manifest.
    mf = utils.get_test_data('scaled/mf.csv')
    shutil.copyfile(mf, runtmp.output('mf.csv'))

    runtmp.sourmash('sig', 'fileinfo', 'mf.csv')

    out = runtmp.last_result.out
    print(out)

    assert '15 sketches with DNA, k=31, scaled=10000           717 total hashes' in out
    assert 'num signatures: 15' in out
    assert 'has manifest? yes' in out
    assert 'is database? yes' in out
    assert 'path filetype: StandaloneManifestIndex' in out


def test_sig_fileinfo_9_sqldb_make(runtmp):
    # make a sqldb and run fileinfo on it
    gcf_all = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    sqldb = runtmp.output('some.sqldb')

    runtmp.sourmash('sig', 'cat', '-k', '31', *gcf_all, '-o', sqldb)

    runtmp.sourmash('sig', 'fileinfo', sqldb)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "12 sketches with DNA, k=31, scaled=10000           4540 total hashes" in out


def test_sig_fileinfo_9_sqldb_exists(runtmp):
    # run fileinfo on existing sqldb
    sqldb = utils.get_test_data('sqlite/index.sqldb')
    runtmp.sourmash('sig', 'fileinfo', sqldb)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "path filetype: SqliteIndex" in out
    assert "2 sketches with DNA, k=31, scaled=1000             10415 total hashes" in out


def test_sig_fileinfo_9_sql_manifest(runtmp):
    # run fileinfo on existing sqldb
    sqldb = utils.get_test_data('sqlite/prot.sqlmf')
    runtmp.sourmash('sig', 'fileinfo', sqldb)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "path filetype: StandaloneManifestIndex" in out
    assert "num signatures: 7" in out
    assert "1 sketches with DNA, k=31, scaled=1000             5238 total hashes" in out
    assert "2 sketches with hp, k=19, scaled=100               5184 total hashes" in out
    assert "2 sketches with dayhoff, k=19, scaled=100          7945 total hashes" in out
    assert "2 sketches with protein, k=19, scaled=100          8214 total hashes" in out


def test_sig_fileinfo_9_sql_lca_db(runtmp):
    # run fileinfo on existing sqldb
    sqldb = utils.get_test_data('sqlite/lca.sqldb')
    runtmp.sourmash('sig', 'fileinfo', sqldb)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "path filetype: LCA_SqliteDatabase" in out
    assert "2 sketches with DNA, k=31, scaled=10000            1431 total hashes" in out
