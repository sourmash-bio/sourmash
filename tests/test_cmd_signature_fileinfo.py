"""
Tests for the 'sourmash signature fileinfo' command line.
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


def test_fileinfo_1_sig(runtmp):
    c = runtmp

    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')

    shutil.copyfile(sig47, runtmp.output('sig47.sig'))
    c.run_sourmash('sig', 'fileinfo', 'sig47.sig')

    out = c.last_result.out
    print(c.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: sig47.sig
is database? no
has manifest? yes
is nonempty? yes
num signatures: 1
5177 total hashes
abundance information available: no
ksizes present: 31
moltypes present: DNA
scaled vals present: 1000
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_2_lca(runtmp):
    c = runtmp

    # get basic info on a signature
    prot = utils.get_test_data('prot/protein.lca.json.gz')

    shutil.copyfile(prot, runtmp.output('protein.lca.json.gz'))
    c.run_sourmash('sig', 'fileinfo', 'protein.lca.json.gz')

    out = c.last_result.out
    print(c.last_result.out)

    expected_output = """\
path filetype: LCA_Database
location: protein.lca.json.gz
is database? yes
has manifest? no
is nonempty? yes
num signatures: 2
8214 total hashes
abundance information available: no
ksizes present: 19
moltypes present: protein
scaled vals present: 100
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_3_sbt_zip(runtmp):
    c = runtmp

    # get basic info on a signature
    prot = utils.get_test_data('prot/protein.sbt.zip')

    shutil.copyfile(prot, runtmp.output('protein.sbt.zip'))
    c.run_sourmash('sig', 'fileinfo', 'protein.sbt.zip')

    out = c.last_result.out
    print(c.last_result.out)

    #abundance information available: no @CTB
    expected_output = """\
path filetype: SBT
location: protein.sbt.zip
is database? yes
has manifest? yes
is nonempty? yes
num signatures: 3
8214 total hashes
ksizes present: 19
moltypes present: protein
scaled vals present: 100
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()


def test_fileinfo_4_zip(runtmp):
    c = runtmp

    # get basic info on a signature
    prot = utils.get_test_data('prot/all.zip')

    shutil.copyfile(prot, runtmp.output('all.zip'))
    c.run_sourmash('sig', 'fileinfo', 'all.zip')

    out = c.last_result.out
    print(c.last_result.out)

    # 'location' will be fully resolved, ignore it for now
    expected_output = f"""\
path filetype: ZipFileLinearIndex
is database? yes
has manifest? yes
is nonempty? yes
num signatures: 8
31758 total hashes
abundance information available: yes
ksizes present: 19, 31
moltypes present: DNA, dayhoff, hp, protein
scaled vals present: 100, 1000
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_5_dir(runtmp):
    c = runtmp

    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')

    os.mkdir(runtmp.output('subdir'))

    shutil.copyfile(sig47, runtmp.output('subdir/sig47.sig'))
    c.run_sourmash('sig', 'fileinfo', 'subdir/')

    out = c.last_result.out
    print(c.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: subdir/
is database? no
has manifest? yes
is nonempty? yes
num signatures: 1
5177 total hashes
abundance information available: no
ksizes present: 31
moltypes present: DNA
scaled vals present: 1000
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


def test_fileinfo_6_pathlist(runtmp):
    c = runtmp

    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')
    shutil.copyfile(sig47, runtmp.output("47.fa.sig"))

    with open(runtmp.output('pathlist.txt'), 'wt') as fp:
        fp.write("47.fa.sig\n")

    c.run_sourmash('sig', 'fileinfo', 'pathlist.txt')

    out = c.last_result.out
    print(c.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: pathlist.txt
is database? no
has manifest? yes
is nonempty? yes
num signatures: 1
5177 total hashes
abundance information available: no
ksizes present: 31
moltypes present: DNA
scaled vals present: 1000
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out


@pytest.mark.parametrize("db", ['v6.sbt.json', 'v5.sbt.json', 'v4.sbt.json',
                                'v3.sbt.json', 'v2.sbt.json', 'v1.sbt.json'])
def test_fileinfo_7_sbt_json(runtmp, db):
    c = runtmp

    # get basic info on an SBT json file
    dbfile = utils.get_test_data(db)

    c.run_sourmash('sig', 'fileinfo', dbfile)

    out = c.last_result.out
    print(c.last_result.out)

    #abundance information available: no @CTB
    expected_output = f"""\
path filetype: SBT
location: {dbfile}
is database? yes
has manifest? no
is nonempty? yes
num signatures: 13
3500 total hashes
abundance information available: no
ksizes present: 31
moltypes present: DNA
no scaled sketches present
num vals present: 500
""".splitlines()
    for line in expected_output:
        assert line.strip() in out, line.strip()
