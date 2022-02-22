"""
Tests for the 'sourmash signature fileinfo' command line.
"""
import csv
import shutil
import os
import glob

import pytest
import screed
import yaml, json

import sourmash_tst_utils as utils
import sourmash
from sourmash.signature import load_signatures
from sourmash.manifest import CollectionManifest
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
5177 total hashes
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000             5177 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_1_sig_abund(runtmp):
    # get basic info on a signature with abundance
    sig47 = utils.get_test_data('47.abunds.fa.sig')

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
8214 total hashes
summary of sketches:
   2 sketches with protein, k=19, scaled=100          8214 total hashes
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
num signatures: 3
8214 total hashes
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
    expected_output = f"""\
path filetype: ZipFileLinearIndex
is database? yes
has manifest? yes
num signatures: 8
31758 total hashes
summary of sketches:
   2 sketches with dayhoff, k=19, scaled=100          7945 total hashes
   2 sketches with hp, k=19, scaled=100               5184 total hashes
   2 sketches with protein, k=19, scaled=100          8214 total hashes
   2 sketches with DNA, k=31, scaled=1000             10415 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_4_zip_yaml_out(runtmp):
    # check --yaml-out
    prot = utils.get_test_data('prot/all.zip')

    shutil.copyfile(prot, runtmp.output('all.zip'))
    runtmp.run_sourmash('sig', 'fileinfo', 'all.zip', '--yaml-out')

    out = runtmp.last_result.out
    print(runtmp.last_result.out)

    # 'location' will be fully resolved, ignore it for now
    expected_output = f"""\
has_manifest: true
is_database: true
num_sketches: 8
path_filetype: ZipFileLinearIndex
sketch_info:
- abund: false
  count: 2
  ksize: 19
  moltype: dayhoff
  n_hashes: 7945
  num: 0
  scaled: 100
- abund: false
  count: 2
  ksize: 19
  moltype: hp
  n_hashes: 5184
  num: 0
  scaled: 100
- abund: false
  count: 2
  ksize: 19
  moltype: protein
  n_hashes: 8214
  num: 0
  scaled: 100
- abund: false
  count: 2
  ksize: 31
  moltype: DNA
  n_hashes: 10415
  num: 0
  scaled: 1000
total_hashes: 31758
""".splitlines()
    for line in expected_output:
        assert line.strip() in out

    # should succeed as loading as valid YAML, too
    vals = yaml.safe_load(out)

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
    expected_output = f"""\
path filetype: ZipFileLinearIndex
is database? yes
has manifest? yes
num signatures: 8
26581 total hashes
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
5177 total hashes
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
5177 total hashes
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
num signatures: 13
3500 total hashes
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
3409 total hashes
summary of sketches:
   1 sketches with protein, k=19, scaled=100          3409 total hashes
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()


def test_sig_fileinfo_does_not_exit(runtmp):
    # test on file that does not exist
    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'fileinfo', 'does-not-exist')

    assert "Cannot open 'does-not-exist'." in runtmp.last_result.err
