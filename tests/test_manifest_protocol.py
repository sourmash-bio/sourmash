"""
Tests for the 'CollectionManifest' class and protocol. All subclasses
of BaseCollectionManifest should support this functionality.
"""

import pytest
import sourmash_tst_utils as utils

import sourmash
from sourmash.manifest import BaseCollectionManifest, CollectionManifest
from sourmash.index.sqlite_index import SqliteCollectionManifest


def build_simple_manifest(runtmp):
    # load and return the manifest from prot/all.zip
    filename = utils.get_test_data('prot/all.zip')
    idx = sourmash.load_file_as_index(filename)
    mf = idx.manifest
    assert len(mf) == 8
    return mf


def build_sqlite_manifest(runtmp):
    # return the manifest from prot/all.zip
    filename = utils.get_test_data('prot/all.zip')
    idx = sourmash.load_file_as_index(filename)
    mf = idx.manifest

    # build sqlite manifest from this 'un
    mfdb = runtmp.output('test.sqlmf')
    return SqliteCollectionManifest.load_from_manifest(mf, dbfile=mfdb)
    

def save_load_manifest(runtmp):
    # save/load the manifest from a CSV.
    mf = build_simple_manifest(runtmp)

    mf_csv = runtmp.output('mf.csv')
    mf.write_to_filename(mf_csv)

    load_mf = CollectionManifest.load_from_filename(mf_csv)
    return load_mf
    

@pytest.fixture(params=[build_simple_manifest,
                        save_load_manifest,
                        build_sqlite_manifest])
def manifest_obj(request, runtmp):
    build_fn = request.param

    return build_fn(runtmp)


###
### generic CollectionManifeset tests go here
###

def test_manifest_len(manifest_obj):
    # check that 'len' works
    assert len(manifest_obj) == 8


def test_manifest_rows(manifest_obj):
    # check that '.rows' property works
    rows = list(manifest_obj.rows)
    assert len(rows) == 8

    required_keys = set(BaseCollectionManifest.required_keys)
    for row in rows:
        kk = set(row.keys())
        assert required_keys.issubset(kk)


def test_manifest_bool(manifest_obj):
    # check that 'bool' works
    assert bool(manifest_obj)


def test_make_manifest_row(manifest_obj):
    # build a manifest row from a signature
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sig47)

    row = manifest_obj.make_manifest_row(ss, 'foo', include_signature=False)
    assert not 'signature' in row
    assert row['internal_location'] == 'foo'

    assert row['md5'] == ss.md5sum()
    assert row['md5short'] == ss.md5sum()[:8]
    assert row['ksize'] == 31
    assert row['moltype'] == 'DNA'
    assert row['num'] == 0
    assert row['scaled'] == 1000
    assert row['n_hashes'] == len(ss.minhash)
    assert not row['with_abundance']
    assert row['name'] == ss.name
    assert row['filename'] == ss.filename

    
def test_manifest_create_manifest(manifest_obj):
    # test the 'create_manifest' method
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sig47)

    def yield_sigs():
        yield ss, 'fiz'

    new_mf = manifest_obj.create_manifest(yield_sigs(),
                                          include_signature=False)
    assert len(new_mf) == 1
    new_row = list(new_mf.rows)[0]
    
    row = manifest_obj.make_manifest_row(ss, 'fiz', include_signature=False)

    required_keys = BaseCollectionManifest.required_keys
    for k in required_keys:
        assert new_row[k] == row[k], k


def test_manifest_select_to_manifest(manifest_obj):
    # do some light testing of 'select_to_manifest'
    new_mf = manifest_obj.select_to_manifest(moltype='DNA')
    assert len(new_mf) == 2


def test_manifest_locations(manifest_obj):
    # check the 'locations' method
    locs = set(['dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
                'dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig',
                'hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
                'hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig',
                'protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
                'protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig',
                'dna-sig.noext',
                'dna-sig.sig.gz']
               )
    assert set(manifest_obj.locations()) == locs


def test_manifest___contains__(manifest_obj):
    # check the 'in' operator
    sigfile = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    ss = sourmash.load_one_signature(sigfile)

    assert ss in manifest_obj

    sigfile2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sigfile2, ksize=31)
    assert ss2 not in manifest_obj


def test_manifest_to_picklist(manifest_obj):
    # test 'to_picklist'
    picklist = manifest_obj.to_picklist()
    mf = manifest_obj.select_to_manifest(picklist=picklist)

    assert mf == manifest_obj


def test_manifest_filter_rows(manifest_obj):
    # test filter_rows
    filter_fn = lambda x: 'OS223' in x['name']

    mf = manifest_obj.filter_rows(filter_fn)

    assert len(mf) == 1
    row = list(mf.rows)[0]
    assert row['name'] == 'NC_011663.1 Shewanella baltica OS223, complete genome'


def test_manifest_filter_cols(manifest_obj):
    # test filter_rows
    col_filter_fn = lambda x: 'OS223' in x[0]

    mf = manifest_obj.filter_on_columns(col_filter_fn, ['name'])

    assert len(mf) == 1
    row = list(mf.rows)[0]
    assert row['name'] == 'NC_011663.1 Shewanella baltica OS223, complete genome'


def test_manifest_iadd(manifest_obj):
    # test the 'create_manifest' method
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sig47)

    def yield_sigs():
        yield ss, 'fiz'

    new_mf = manifest_obj.create_manifest(yield_sigs(),
                                          include_signature=False)
    assert len(new_mf) == 1

    new_mf += manifest_obj
    assert len(new_mf) == len(manifest_obj) + 1


def test_manifest_add(manifest_obj):
    # test the 'create_manifest' method
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sig47)

    def yield_sigs():
        yield ss, 'fiz'

    new_mf = manifest_obj.create_manifest(yield_sigs(),
                                          include_signature=False)
    assert len(new_mf) == 1

    new_mf2 = new_mf + manifest_obj
    assert len(new_mf2) == len(manifest_obj) + len(new_mf)
