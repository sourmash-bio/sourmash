"""
Tests for the 'sourmash lca' command line and high level API.
"""
import os
import shutil
import csv
import pytest
import glob

import sourmash_tst_utils as utils
import sourmash
from sourmash import load_one_signature, SourmashSignature, sourmash_args

from sourmash.search import make_jaccard_search_query
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import LineagePair
from sourmash.picklist import SignaturePicklist, PickStyle
from sourmash_tst_utils import SourmashCommandFailed


def test_api_create_search():
    # create a database and then search for result.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    assert len(lca_db) == 0
    assert not lca_db

    count = lca_db.insert(ss)
    assert count == len(ss.minhash)

    assert len(lca_db) == 1
    assert lca_db

    results = lca_db.search(ss, threshold=0.0)
    print(results)
    assert len(results) == 1
    (similarity, match, filename) = results[0]
    assert match.minhash == ss.minhash


def test_api_find_picklist_select():
    # does 'find' respect picklists?

    sig47 = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                        ksize=31)
    sig63 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                        ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(sig47)
    lca_db.insert(sig63)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['09a08691'])

    # run a 'find' with sig63, should find 47 and 63 both.
    search_obj = make_jaccard_search_query(do_containment=True, threshold=0.0)
    results = list(lca_db.find(search_obj, sig63))
    print(results)
    assert len(results) == 2

    # now, select on picklist and do another find...
    lca_db = lca_db.select(picklist=picklist)
    results = list(lca_db.find(search_obj, sig63))
    print(results)
    assert len(results) == 1

    # and check that it is the expected one!
    ss = results[0].signature
    assert ss.minhash.ksize == 31
    assert ss.md5sum().startswith('09a08691c')


def test_api_find_picklist_select_exclude():
    # does 'find' respect picklists?

    sig47 = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                        ksize=31)
    sig63 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                        ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(sig47)
    lca_db.insert(sig63)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8', pickstyle= PickStyle.EXCLUDE)
    picklist.init(['09a08691'])

    # run a 'find' with sig63, should find 47 and 63 both.
    search_obj = make_jaccard_search_query(do_containment=True, threshold=0.0)
    results = list(lca_db.find(search_obj, sig63))
    print(results)
    assert len(results) == 2

    # now, select on picklist and do another find...
    lca_db = lca_db.select(picklist=picklist)
    results = list(lca_db.find(search_obj, sig63))
    print(results)
    assert len(results) == 1

    # and check that it is the expected one!
    ss = results[0].signature
    assert ss.minhash.ksize == 31
    assert ss.md5sum().startswith('38729c637')


def test_api_create_insert():
    # test some internal implementation stuff: create & then insert a sig.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    ident = ss.name
    assert len(lca_db._ident_to_name) == 1
    assert ident in lca_db._ident_to_name
    assert lca_db._ident_to_name[ident] == ident
    assert len(lca_db._ident_to_idx) == 1
    assert lca_db._ident_to_idx[ident] == 0
    assert len(lca_db._hashval_to_idx) == len(ss.minhash)
    assert len(lca_db._idx_to_ident) == 1
    assert lca_db._idx_to_ident[0] == ident

    set_of_values = set()
    for vv in lca_db._hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    assert not lca_db._idx_to_lid          # no lineage added
    assert not lca_db._lid_to_lineage      # no lineage added


def test_api_create_insert_bad_ksize():
    # can we insert a ksize=21 signature into a ksize=31 DB? hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=21, scaled=1000)
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_bad_ident():
    # can we insert a signature with no/empty ident?
    ss1 = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                      ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)
    ss1 = ss1.to_mutable()
    ss2 = ss2.to_mutable()

    ss1.name = ''
    ss1.filename = ''
    ss2.name = ''
    ss2.filename = ''

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss1)
    lca_db.insert(ss2)
    # SUCCESS!
    # would fail, previously :)


def test_api_create_insert_bad_scaled():
    # can we insert a scaled=1000 signature into a scaled=500 DB?
    # hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    assert ss.minhash.scaled == 1000

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=500)
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_bad_moltype():
    # can we insert a DNAsignature into a protein DB?
    # hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    assert ss.minhash.moltype == 'DNA'

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=500, moltype='protein')
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_ident():
    # test some internal implementation stuff: signature inserted with
    # different ident than name.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, ident='foo')

    ident = 'foo'
    assert len(lca_db._ident_to_name) == 1
    assert ident in lca_db._ident_to_name
    assert lca_db._ident_to_name[ident] == ss.name
    assert len(lca_db._ident_to_idx) == 1
    assert lca_db._ident_to_idx[ident] == 0
    assert len(lca_db._hashval_to_idx) == len(ss.minhash)
    assert len(lca_db._idx_to_ident) == 1
    assert lca_db._idx_to_ident[0] == ident

    set_of_values = set()
    for vv in lca_db._hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    assert not lca_db._idx_to_lid          # no lineage added
    assert not lca_db._lid_to_lineage      # no lineage added
    assert not lca_db._lineage_to_lid
    assert not lca_db._lid_to_idx


def test_api_create_insert_two():
    # check internal details if multiple signatures are inserted.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, ident='foo')
    lca_db.insert(ss2, ident='bar')

    ident = 'foo'
    ident2 = 'bar'
    assert len(lca_db._ident_to_name) == 2
    assert ident in lca_db._ident_to_name
    assert ident2 in lca_db._ident_to_name
    assert lca_db._ident_to_name[ident] == ss.name
    assert lca_db._ident_to_name[ident2] == ss2.name

    assert len(lca_db._ident_to_idx) == 2
    assert lca_db._ident_to_idx[ident] == 0
    assert lca_db._ident_to_idx[ident2] == 1

    combined_mins = set(ss.minhash.hashes.keys())
    combined_mins.update(set(ss2.minhash.hashes.keys()))
    assert len(lca_db._hashval_to_idx) == len(combined_mins)

    assert len(lca_db._idx_to_ident) == 2
    assert lca_db._idx_to_ident[0] == ident
    assert lca_db._idx_to_ident[1] == ident2

    set_of_values = set()
    for vv in lca_db._hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 2
    assert set_of_values == { 0, 1 }

    assert not lca_db._idx_to_lid          # no lineage added
    assert not lca_db._lid_to_lineage      # no lineage added
    assert not lca_db._lineage_to_lid
    assert not lca_db._lid_to_idx


def test_api_create_insert_w_lineage():
    # test some internal implementation stuff - insert signature w/lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = ((LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')))

    lca_db.insert(ss, lineage=lineage)

    # basic ident stuff
    ident = ss.name
    assert len(lca_db._ident_to_name) == 1
    assert ident in lca_db._ident_to_name
    assert lca_db._ident_to_name[ident] == ident
    assert len(lca_db._ident_to_idx) == 1
    assert lca_db._ident_to_idx[ident] == 0
    assert len(lca_db._hashval_to_idx) == len(ss.minhash)
    assert len(lca_db._idx_to_ident) == 1
    assert lca_db._idx_to_ident[0] == ident

    # all hash values added
    set_of_values = set()
    for vv in lca_db._hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    # check lineage stuff
    assert len(lca_db._idx_to_lid) == 1
    assert lca_db._idx_to_lid[0] == 0
    assert len(lca_db._lid_to_lineage) == 1
    assert lca_db._lid_to_lineage[0] == lineage
    assert lca_db._lid_to_idx[0] == { 0 }

    assert len(lca_db._lineage_to_lid) == 1
    assert lca_db._lineage_to_lid[lineage] == 0


def test_api_create_insert_w_bad_lineage():
    # test some internal implementation stuff - insert signature w/bad lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = ([LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')],)

    with pytest.raises(ValueError):
        lca_db.insert(ss, lineage=lineage)


def test_api_create_insert_w_bad_lineage_2():
    # test some internal implementation stuff - insert signature w/bad lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = 1 # something non-iterable...

    with pytest.raises(ValueError):
        lca_db.insert(ss, lineage=lineage)


def test_api_create_gather():
    # create a database, and then run gather on it.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    result = lca_db.best_containment(ss, threshold_bp=0)
    print(result)
    assert result
    (similarity, match, filename) = result
    assert match.minhash == ss.minhash


def test_api_add_genome_lineage():
    # LCA_Databases can store/retrieve arbitrary lineages/taxonomies.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    lineage = ((LineagePair('rank1', 'name1'),
               (LineagePair('rank2', 'name2'))))

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, lineage=lineage)

    somehash = next(iter(ss.minhash.hashes.keys()))

    lineages = lca_db.get_lineage_assignments(somehash)
    assert len(lineages) == 1
    assert lineage in lineages


def test_api_insert_update():
    # check that cached parts of LCA_Database are updated when a new
    # signature is inserted.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    all_mh = [ x.minhash for x in lca_db.signatures() ]
    assert ss.minhash in all_mh

    # see decorator @cached_property
    assert hasattr(lca_db, '_cache')
    assert lca_db._cache
    # inserting a signature should delete the cache
    lca_db.insert(ss2)
    assert not hasattr(lca_db, '_cache')

    # check that it's rebuilt etc. etc.
    all_mh = [ x.minhash for x in lca_db.signatures() ]
    assert ss.minhash in all_mh
    assert ss2.minhash in all_mh


def test_api_insert_retrieve_check_name():
    # check that signatures retrieved from LCA_Database objects have the
    # right name.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    sigs = list(lca_db.signatures())
    assert len(sigs) == 1
    retrieved_sig = sigs[0]
    assert retrieved_sig.name == ss.name
    assert retrieved_sig.minhash == ss.minhash


def test_api_create_insert_two_then_scale():
    # construct database, THEN downsample
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)
    lca_db.insert(ss2)

    # downsample everything to 5000
    lca_db.downsample_scaled(5000)

    minhash = ss.minhash.downsample(scaled=5000)
    minhash2 = ss2.minhash.downsample(scaled=5000)

    # & check...
    combined_mins = set(minhash.hashes.keys())
    combined_mins.update(set(minhash2.hashes.keys()))
    assert len(lca_db._hashval_to_idx) == len(combined_mins)


def test_api_create_insert_two_then_scale_then_add():
    # construct database, THEN downsample, then add another
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    # downsample everything to 5000
    lca_db.downsample_scaled(5000)

    # insert another after downsample
    lca_db.insert(ss2)

    # now test -
    ss = ss.to_mutable()
    ss.minhash = ss.minhash.downsample(scaled=5000)

    ss2 = ss2.to_mutable()
    ss2.minhash = ss2.minhash.downsample(scaled=5000)

    # & check...
    combined_mins = set(ss.minhash.hashes.keys())
    combined_mins.update(set(ss2.minhash.hashes.keys()))
    assert len(lca_db._hashval_to_idx) == len(combined_mins)


def test_api_create_insert_scale_two():
    # downsample while constructing database
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    # downsample to 5000 while inserting:
    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=5000)
    count = lca_db.insert(ss)
    assert count == 1037
    assert count == len(ss.minhash.downsample(scaled=5000))
    lca_db.insert(ss2)

    # downsample sigs to 5000
    minhash = ss.minhash.downsample(scaled=5000)
    minhash2 = ss2.minhash.downsample(scaled=5000)

    # & check...
    combined_mins = set(minhash.hashes.keys())
    combined_mins.update(set(minhash2.hashes.keys()))
    assert len(lca_db._hashval_to_idx) == len(combined_mins)


def test_load_single_db():
    filename = utils.get_test_data('lca/delmont-1.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    print(db)

    assert ksize == 31
    assert scaled == 10000


def test_load_single_db_empty(runtmp):
    # test load_single_database on an empty file; should raise ValueError
    empty = runtmp.output('empty.lca.json')

    with open(empty, "wt") as fp:
        pass

    with pytest.raises(ValueError) as exc:
        db, ksize, scaled = lca_utils.load_single_database(empty)

    assert f"'{empty}' is not an LCA database file." in str(exc.value)


def test_databases():
    filename1 = utils.get_test_data('lca/delmont-1.lca.json')
    filename2 = utils.get_test_data('lca/delmont-2.lca.json')
    dblist, ksize, scaled = lca_utils.load_databases([filename1, filename2])

    print(dblist)

    assert len(dblist) == 2
    assert ksize == 31
    assert scaled == 10000


def test_databases_load_fail_on_no_JSON():
    filename1 = utils.get_test_data('prot/protein.zip')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not an LCA database file." in err


def test_databases_load_fail_on_dir():
    filename1 = utils.get_test_data('lca')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not a file and cannot be loaded as an LCA database" in err
    assert not 'found 0 matches total;' in err


def test_databases_load_fail_on_not_exist():
    filename1 = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not a file and cannot be loaded as an LCA database" in err
    assert not 'found 0 matches total;' in err

def test_db_repr():
    filename = utils.get_test_data('lca/delmont-1.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    assert repr(db) == "LCA_Database('{}')".format(filename)


def test_lca_index_signatures_method():
    # test 'signatures' method from base class Index
    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    siglist = list(db.signatures())
    assert len(siglist) == 2


def test_lca_index_select():
    # test 'select' method from Index base class.

    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    xx = db.select(ksize=31)
    assert xx == db

    xx = db.select(moltype='DNA')
    assert xx == db

    xx = db.select(abund=False)
    assert xx == db

    with pytest.raises(ValueError):
        db.select(ksize=21)

    with pytest.raises(ValueError):
        db.select(moltype='protein')

    with pytest.raises(ValueError):
        db.select(abund=True)


def test_lca_index_select_picklist():
    # test 'select' method from Index base class with a picklist.

    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['50a92740'])

    xx = db.select(picklist=picklist)
    assert xx == db

    siglist = list(db.signatures())
    assert len(siglist) == 1
    ss = siglist[0]
    assert ss.md5sum().startswith('50a92740')
    assert ss.minhash.ksize == 31


def test_lca_index_find_picklist_check_overlap():
    # make sure 'find' works for picklists that exclude relevant signatures
    # (bug #1638)

    query_fn = utils.get_test_data('47.fa.sig')
    query_sig = sourmash.load_one_signature(query_fn, ksize=31)
    db_fn = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(db_fn)

    # construct a picklist...
    picklist = SignaturePicklist('ident')
    picklist.init(['NC_009665.1'])

    xx = db.select(picklist=picklist)
    assert xx == db

    results = list(db.search(query_sig, threshold=0.1))
    assert len(results) == 1


def test_lca_index_select_picklist_exclude():
    # test 'select' method from Index base class with a picklist.

    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8', pickstyle=PickStyle.EXCLUDE)
    picklist.init(['50a92740'])

    xx = db.select(picklist=picklist)
    assert xx == db

    siglist = list(db.signatures())
    assert len(siglist) == 1
    ss = siglist[0]
    assert ss.md5sum().startswith('e88dc390')
    assert ss.minhash.ksize == 31


def test_lca_index_select_picklist_twice():
    # test 'select' method from Index base class with a picklist.

    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['50a92740'])

    xx = db.select(picklist=picklist)
    assert xx == db

    with pytest.raises(ValueError) as exc:
        xx = db.select(picklist=picklist)

    assert "we do not (yet) support multiple picklists for LCA databases" in str(exc)



def test_search_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    results = db.search(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    minhash = sig.minhash.downsample(scaled=10000)
    assert minhash == match_sig.minhash


def test_search_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    sig = sig.to_mutable()
    sig.minhash = sig.minhash.downsample(scaled=100000)

    results = db.search(sig, threshold=.01, ignore_abundance=True)
    print(results)
    assert results[0].score == 1.0
    match = results[0].signature

    orig_sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    assert orig_sig.minhash.jaccard(match.minhash, downsample=True) == 1.0


def test_gather_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    result = db.best_containment(sig, threshold=.01, ignore_abundance=True)
    match_sig = result[1]

    minhash = sig.minhash.downsample(scaled=10000)
    assert minhash == match_sig.minhash


def test_gather_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    sig_minhash = sig.minhash.downsample(scaled=100000)

    result = db.best_containment(sig, threshold=.01, ignore_abundance=True)
    match_sig = result[1]

    minhash = match_sig.minhash.downsample(scaled=100000)
    assert sig_minhash == minhash


def test_db_lineage_to_lid():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db._lineage_to_lid
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)

    lin1 = items[0][0][-1]
    assert lin1.rank == 'strain'
    assert lin1.name == 'Shewanella baltica OS185'
    lin1 = items[1][0][-1]
    assert lin1.rank == 'strain'
    assert lin1.name == 'Shewanella baltica OS223'


def test_db_lid_to_idx():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db._lid_to_idx
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, {32}), (48, {48})]


def test_db_idx_to_ident():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db._idx_to_ident
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, 'NC_009665'), (48, 'NC_011663')]


## command line tests


def test_run_sourmash_lca():
    status, out, err = utils.runscript('sourmash', ['lca'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_basic_index(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, 'delmont-1', input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db), lca_db

    assert 'Building LCA database with ksize=31 scaled=10000 moltype=DNA' in runtmp.last_result.err
    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err


def test_basic_index_twice(runtmp, lca_db_format):
    # run 'lca index' twice.
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, 'delmont-1', input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    with pytest.raises(SourmashCommandFailed):
        cmd = ['lca', 'index', taxcsv, 'delmont-1', input_sig, '-F', lca_db_format]
        runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'already exists. Not overwriting.' in runtmp.last_result.err


def test_basic_index_bad_spreadsheet(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/bad-spreadsheet.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db), lca_db

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err


def test_basic_index_broken_spreadsheet(runtmp, lca_db_format):
    # duplicate identifiers in this spreadsheet
    taxcsv = utils.get_test_data('lca/bad-spreadsheet-2.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(*cmd)

    assert runtmp.last_result.status != 0
    assert "multiple lineages for identifier TARA_ASE_MAG_00031" in runtmp.last_result.err


def test_basic_index_too_many_strains_too_few_species(runtmp, lca_db_format):
    # explicit test for #841, where 'n_species' wasn't getting counted
    # if lineage was at strain level resolution.
    taxcsv = utils.get_test_data('lca/podar-lineage.csv')
    input_sig = utils.get_test_data('47.fa.sig')
    lca_db = runtmp.output(f'out.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig,
            '-C', '3', '--split-identifiers', '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    assert not 'error: fewer than 20% of lineages' in runtmp.last_result.err
    assert runtmp.last_result.status == 0


def test_basic_index_too_few_species(runtmp, lca_db_format):
    # spreadsheets with too few species should be flagged, unless -f specified
    taxcsv = utils.get_test_data('lca/tully-genome-sigs.classify.csv')

    # (these don't really matter, should break on load spreadsheet)
    input_sig = utils.get_test_data('47.fa.sig')
    lca_db = runtmp.output(f'out.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-C', '3',
           '-F', lca_db_format]
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(*cmd)

    assert not '"ERROR: fewer than 20% of lineages have species-level resolution' in runtmp.last_result.err
    assert runtmp.last_result.status != 0


def test_basic_index_require_taxonomy(runtmp, lca_db_format):
    # no taxonomy in here
    taxcsv = utils.get_test_data('lca/bad-spreadsheet-3.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', '--require-taxonomy', taxcsv, lca_db, input_sig,
           '-F', lca_db_format]
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(*cmd)

    assert runtmp.last_result.status != 0
    assert "ERROR: no hash values found - are there any signatures?" in runtmp.last_result.err


def test_basic_index_column_start(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-3.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', '-C', '3', taxcsv, lca_db, input_sig,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err


def test_index_empty_sketch_name(runtmp, lca_db_format):
    c = runtmp

    # create two signatures with empty 'name' attributes
    cmd = ['sketch', 'dna', utils.get_test_data('genome-s12.fa.gz'),
           utils.get_test_data('genome-s11.fa.gz')]
    c.run_sourmash(*cmd)

    sig1 = c.output('genome-s11.fa.gz.sig')
    assert os.path.exists(sig1)
    sig2 = c.output('genome-s12.fa.gz.sig')
    assert os.path.exists(sig2)

    outfile = f'zzz.lca.{lca_db_format}'

    # can we insert them both?
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    cmd = ['lca', 'index', taxcsv, outfile, sig1, sig2, '-F', lca_db_format]
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output(outfile))

    print(c.last_result.out)
    print(c.last_result.err)
    assert 'WARNING: no lineage provided for 2 sig' in c.last_result.err


def test_basic_index_and_classify_with_tsv_and_gz(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-1.tsv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    if lca_db_format == 'json':
        lca_db = runtmp.output(f'delmont-1.lca.json.gz')
    else:
        lca_db = runtmp.output(f'delmont-1.lca.sql')

    cmd = ['lca', 'index', '--tabs', '--no-header', taxcsv, lca_db, input_sig,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_basic_index_and_classify(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_basic_index_and_classify_dup_lineage(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00007.sig')
    input_sig2 = utils.get_test_data('lca/TARA_ANW_MAG_00005.sig')
    lca_db = runtmp.output(f'delmont-dup.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format, '-f']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'TARA_ASE_MAG_00007,found,Bacteria,Proteobacteria,Gammaproteobacteria,,,,,' in runtmp.last_result.out

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'TARA_ANW_MAG_00005,found,Bacteria,Proteobacteria,Gammaproteobacteria,,,,,' in runtmp.last_result.out


def test_index_traverse(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    in_dir = runtmp.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'index', taxcsv, lca_db, in_dir, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err
    assert 'WARNING: 1 duplicate signatures.' not in runtmp.last_result.err


def test_index_traverse_force(runtmp, lca_db_format):
    c = runtmp
    # test the use of --force to load all files, not just .sig
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    # name signature .txt instead of .sig:
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.txt'))

    # use --force
    cmd = ['lca', 'index', taxcsv, lca_db, in_dir, '-f', '-F', lca_db_format]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    err = c.last_result.err
    print(out)
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err
    assert 'WARNING: 1 duplicate signatures.' not in err


def test_index_from_file_cmdline_sig(runtmp, lca_db_format):
    c = runtmp
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    file_list = c.output('sigs.list')
    with open(file_list, 'wt') as fp:
        print(input_sig, file=fp)

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '--from-file', file_list,
           '-F', lca_db_format]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err
    assert 'WARNING: 1 duplicate signatures.' in err


def test_index_from_file(runtmp, lca_db_format):
    c = runtmp

    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    file_list = c.output('sigs.list')
    with open(file_list, 'wt') as fp:
        print(input_sig, file=fp)

    cmd = ['lca', 'index', taxcsv, lca_db, '--from-file', file_list,
           '-F', lca_db_format]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_index_fail_on_num(runtmp, lca_db_format):
    c = runtmp
    # lca index should yield a decent error message when attempted on 'num'
    sigfile = utils.get_test_data('num/63.fa.sig')
    taxcsv = utils.get_test_data('lca/podar-lineage.csv')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('lca', 'index', taxcsv, f'xxx.lca.{lca_db_format}', sigfile,
                       '-C', '3', '-F', lca_db_format)

    err = c.last_result.err
    print(err)

    assert 'ERROR: cannot insert signature ' in err
    assert 'ERROR: cannot downsample signature; is it a scaled signature?' in err


def test_index_traverse_real_spreadsheet_no_report(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-f',
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in runtmp.last_result.err
    assert 'WARNING: no signatures for 956 spreadsheet rows.' in runtmp.last_result.err
    assert 'WARNING: 105 unused lineages.' in runtmp.last_result.err
    assert '(You can use --report to generate a detailed report.)' in runtmp.last_result.err


def test_index_traverse_real_spreadsheet_report(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')
    report_loc = runtmp.output('report.txt')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '--report',
            report_loc, '-f', '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in runtmp.last_result.err
    assert 'WARNING: no signatures for 956 spreadsheet rows.' in runtmp.last_result.err
    assert 'WARNING: 105 unused lineages.' in runtmp.last_result.err
    assert '(You can use --report to generate a detailed report.)' not in runtmp.last_result.err
    assert os.path.exists(report_loc)


def test_single_classify(runtmp):
    # run a basic 'classify', check output.
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_single_classify_zip_query(runtmp):
    # run 'classify' with a query in a zipfile
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    query_ss = sourmash.load_one_signature(input_sig, ksize=31)
    query_zipfile = runtmp.output('query.zip')
    with sourmash_args.SaveSignaturesToLocation(query_zipfile) as save_sig:
        save_sig.add(query_ss)

    cmd = ['lca', 'classify', '--db', db1, '--query', query_zipfile]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_single_classify_to_output(runtmp):
    db1 = utils.get_test_data(f'lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    cmd = ['lca', 'classify', '--db', db1, '--query', input_sig,
            '-o', runtmp.output('outfile.txt')]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(runtmp.output('outfile.txt'), 'rt') as fp:
        outdata = fp.read()
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in outdata
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_single_classify_to_output_no_name(runtmp):
    db1 = utils.get_test_data(f'lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    ss = sourmash.load_one_signature(input_sig, ksize=31)

    outsig_filename = runtmp.output('q.sig')
    with open(outsig_filename, 'wt') as fp:
        # remove name from signature here --
        new_sig = sourmash.SourmashSignature(ss.minhash, filename='xyz')
        sourmash.save_signatures([new_sig], fp)

    cmd = ['lca', 'classify', '--db', db1, '--query', outsig_filename,
            '-o', runtmp.output('outfile.txt')]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    with open(runtmp.output('outfile.txt'), 'rt') as fp:
        outdata = fp.read()
    print((outdata,))
    assert 'xyz,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in outdata
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_single_classify_empty(runtmp):
    db1 = utils.get_test_data(f'lca/both.lca.json')
    input_sig = utils.get_test_data('GCF_000005845.2_ASM584v2_genomic.fna.gz.sig')

    cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'GCF_000005845,nomatch,,,,,,,,' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_single_classify_traverse(runtmp):
    db1 = utils.get_test_data(f'lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = runtmp.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_multi_query_classify_traverse(runtmp):
    # both.lca.json is built from both dir and dir2
    db1 = utils.get_test_data(f'lca/both.lca.json')
    dir1 = utils.get_test_data('lca/dir1')
    dir2 = utils.get_test_data('lca/dir2')

    cmd = ['lca', 'classify', '--db', db1, '--query', dir1, dir2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(utils.get_test_data('lca/classify-by-both.csv')) as fp:
        fp_lines = fp.readlines()
        out_lines = runtmp.last_result.out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


@utils.in_tempdir
def test_multi_query_classify_query_from_file(c):
    # both.lca.json is built from both dir and dir2
    db1 = utils.get_test_data('lca/both.lca.json')
    dir1_glob = utils.get_test_data('lca/dir1/*.sig')
    dir1_files = glob.glob(dir1_glob)
    dir2_glob = utils.get_test_data('lca/dir2/*.sig')
    dir2_files = glob.glob(dir2_glob)

    file_list = c.output('file.list')
    with open(file_list, 'wt') as fp:
        print("\n".join(dir1_files), file=fp)
        print("\n".join(dir2_files), file=fp)

    cmd = ['lca', 'classify', '--db', db1, '--query-from-file', file_list]
    c.run_sourmash(*cmd)
    out = c.last_result.out

    with open(utils.get_test_data('lca/classify-by-both.csv')) as fp:
        fp_lines = fp.readlines()
        out_lines = out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


@utils.in_tempdir
def test_multi_query_classify_query_from_file_and_query(c):
    # both.lca.json is built from both dir and dir2
    db1 = utils.get_test_data(f'lca/both.lca.json')
    dir1_glob = utils.get_test_data('lca/dir1/*.sig')
    dir1_files = glob.glob(dir1_glob)
    dir2_glob = utils.get_test_data('lca/dir2/*.sig')
    dir2_files = glob.glob(dir2_glob)

    file_list = c.output('file.list')
    with open(file_list, 'wt') as fp:
        print("\n".join(dir1_files[1:]), file=fp)   # leave off first one
        print("\n".join(dir2_files), file=fp)

    cmd = ['lca', 'classify', '--db', db1, '--query', dir1_files[0],
           '--query-from-file', file_list]
    c.run_sourmash(*cmd)
    out = c.last_result.out

    with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
        fp_lines = fp.readlines()
        out_lines = out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


def test_multi_db_multi_query_classify_traverse(runtmp):
    # two halves of both.lca.json, see above test.
    db1 = utils.get_test_data(f'lca/dir1.lca.json')
    db2 = utils.get_test_data(f'lca/dir2.lca.json')
    dir1 = utils.get_test_data('lca/dir1')
    dir2 = utils.get_test_data('lca/dir2')

    cmd = ['lca', 'classify', '--db', db1, db2, '--query', dir1, dir2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
        fp_lines = fp.readlines()
        out_lines = runtmp.last_result.out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


def test_unassigned_internal_index_and_classify(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-4.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,unassigned,Alteromonadaceae,unassigned,Alteromonas_macleodii' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_unassigned_last_index_and_classify(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-5.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,,,\r\n' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_index_and_classify_internal_unassigned_multi(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    # classify input_sig1
    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,unassigned,unassigned,Alteromonadaceae,,,\r\n' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err

    # classify input_sig2
    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_PSW_MAG_00136,found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus,,\r\n' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 1 LCA databases' in runtmp.last_result.err


def test_classify_majority_vote_1(runtmp, lca_db_format):
    # classify merged signature using lca should yield no results
    c = runtmp

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
                   '-F', lca_db_format)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # merge input_sig1 and input_sig2
    c.run_sourmash('signature', 'merge', input_sig1, input_sig2, '-k', '31', '--flatten', '-o', 'sig1and2.sig')
    sig1and2 = c.output('sig1and2.sig')

    # lca classify should yield no results
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', sig1and2)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'disagree,,,,,,,,' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err



def test_classify_majority_vote_2(runtmp, lca_db_format):
    # classify same signature with same database using --majority
    # should yield results

    c = runtmp

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
                   '-F', lca_db_format)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # merge input_sig1 and input_sig2
    c.run_sourmash('signature', 'merge', input_sig1, input_sig2, '-k', '31', '--flatten', '-o', 'sig1and2.sig')
    sig1and2 = c.output('sig1and2.sig')

    # majority vote classify
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', sig1and2, '--majority')

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err


def test_classify_majority_vote_3(runtmp, lca_db_format):
    # classify signature with nothing in counts
    c = runtmp

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output(f'delmont-1.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
                   '-F', lca_db_format)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # obtain testdata '47.fa.sig'
    testdata1 = utils.get_test_data('47.fa.sig')

    # majority vote classify
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', testdata1, '--majority')

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'nomatch,,,,,,,,' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err


def test_multi_db_classify(runtmp):
    db1 = utils.get_test_data(f'lca/delmont-1.lca.json')
    db2 = utils.get_test_data('lca/delmont-2.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    cmd = ['lca', 'classify', '--db', db1, db2, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in runtmp.last_result.out
    assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,,,,' in runtmp.last_result.out
    assert 'classified 1 signatures total' in runtmp.last_result.err
    assert 'loaded 2 LCA databases' in runtmp.last_result.err


def test_classify_unknown_hashes(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig2, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '(root)' not in runtmp.last_result.out
    assert 'TARA_MED_MAG_00029,found,Archaea,Euryarcheoata,unassigned,unassigned,novelFamily_I' in runtmp.last_result.out


def test_single_summarize(runtmp):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in runtmp.last_result.out


def test_single_summarize_singleton(runtmp):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in runtmp.last_result.out
    assert 'test-data/lca/TARA_ASE_MAG_00031.sig:5b438c6c TARA_ASE_MAG_00031' in runtmp.last_result.out


@utils.in_tempdir
def test_single_summarize_traverse(c):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', in_dir]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'loaded 1 signatures from 1 files total.' in err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out

@utils.in_tempdir
def test_single_summarize_singleton_traverse(c):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', in_dir]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'loaded 1 signatures from 1 files total.' in err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out
    assert 'q.sig:5b438c6c TARA_ASE_MAG_00031' in out


def test_single_summarize_to_output(runtmp):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = runtmp.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,
            '-o', runtmp.output('output.txt')]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(runtmp.output('output.txt'), 'rt') as fp:
        outdata = fp.read()
    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert '200,Bacteria,Proteobacteria,Gammaproteobacteria' in outdata



def test_single_summarize_to_output_check_filename(runtmp):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = runtmp.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', os.path.join(in_dir, 'q.sig'),
            '-o', runtmp.output('output.txt')]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    outdata = open(runtmp.output('output.txt'), 'rt').read()

    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert 'count,superkingdom,phylum,class,order,family,genus,species,strain,filename,sig_name,sig_md5,total_counts\n' in outdata
    assert '200,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii,,'+os.path.join(in_dir, 'q.sig')+',TARA_ASE_MAG_00031,5b438c6c858cdaf9e9b05a207fa3f9f0,200.0\n' in outdata
    print(outdata)


def test_summarize_unknown_hashes_to_output_check_total_counts(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig2, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1,
           '-o', 'out.csv']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '(root)' not in runtmp.last_result.out
    assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in runtmp.last_result.out

    with open(runtmp.output('out.csv'), newline="") as fp:
        r = csv.DictReader(fp)
        rows = list(r)
        pairs = [ (row['count'], row['total_counts']) for row in rows ]
        pairs = [ (float(x), float(y)) for x, y in pairs ]
        pairs = set(pairs)

        assert pairs == { (27.0, 234.0) }


def test_single_summarize_scaled(runtmp):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = runtmp.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,
            '--scaled', '100000']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert '100.0%    27   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales'


def test_single_summarize_scaled_zip_query(runtmp):
    # check zipfile as query
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    query_ss = sourmash.load_one_signature(input_sig, ksize=31)
    query_zipfile = runtmp.output('query.zip')
    with sourmash_args.SaveSignaturesToLocation(query_zipfile) as save_sig:
        save_sig.add(query_ss)

    cmd = ['lca', 'summarize', '--db', db1, '--query', query_zipfile,
            '--scaled', '100000']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 1 signatures from 1 files total.' in runtmp.last_result.err
    assert '100.0%    27   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales'


def test_multi_summarize_with_unassigned_singleton(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1,
            input_sig2, '--ignore-abundance']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 2 signatures from 2 files total.' in runtmp.last_result.err

    out_lines = runtmp.last_result.out.splitlines()
    def remove_line_startswith(x, check=None):
        for line in out_lines:
            if line.startswith(x):
                out_lines.remove(line)
                if check:
                    # make sure the check value is in there
                    assert check in line
                return line
        assert 0, "couldn't find {}".format(x)

    # note, proportions/percentages are now per-file
    remove_line_startswith('100.0%   200   Bacteria ', 'TARA_ASE_MAG_00031.sig:5b438c6c')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta ')
    remove_line_startswith('100.0%  1231   Eukaryota ', 'TARA_PSW_MAG_00136.sig:db50b713')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria ')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae ')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned;Alteromonadaceae ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned;Ostreococcus ')
    assert not out_lines


def test_multi_summarize_with_zip_unassigned_singleton(runtmp, lca_db_format):
    # test summarize on multiple queries, in a zipfile.
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    query_zipfile = runtmp.output('query.zip')
    with sourmash_args.SaveSignaturesToLocation(query_zipfile) as save_sig:
        input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        sig1 = sourmash.load_one_signature(input_sig1, ksize=31)
        input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        sig2 = sourmash.load_one_signature(input_sig2, ksize=31)

        save_sig.add(sig1)
        save_sig.add(sig2)

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', 'query.zip',
           '--ignore-abundance']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 2 signatures from 1 files total.' in runtmp.last_result.err

    out_lines = runtmp.last_result.out.splitlines()
    def remove_line_startswith(x, check=None):
        for line in out_lines:
            if line.startswith(x):
                out_lines.remove(line)
                if check:
                    # make sure the check value is in there
                    assert check in line
                return line
        assert 0, "couldn't find {}".format(x)

    # note, proportions/percentages are now per-file
    remove_line_startswith('100.0%   200   Bacteria ', ':5b438c6c')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta ')
    remove_line_startswith('100.0%  1231   Eukaryota ', ':db50b713')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria ')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae ')
    remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned;Alteromonadaceae ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned ')
    remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned;Ostreococcus ')
    assert not out_lines


def test_summarize_to_root(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig2,
            '--ignore-abundance']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '78.6%    99   Archaea' in runtmp.last_result.out
    assert '21.4%    27   (root)' in runtmp.last_result.out


def test_summarize_unknown_hashes(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig2, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '(root)' not in runtmp.last_result.out
    assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in runtmp.last_result.out


def test_summarize_to_root_abund(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2,
           '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '78.9%   101   Archaea' in runtmp.last_result.out
    assert '21.1%    27   (root)' in runtmp.last_result.out


def test_summarize_unknown_hashes_abund(runtmp, lca_db_format):
    taxcsv = utils.get_test_data('lca-root/tax.csv')
    input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
    lca_db = runtmp.output(f'lca-root.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig2, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '(root)' not in runtmp.last_result.out
    assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in runtmp.last_result.out


@utils.in_thisdir
def test_summarize_abund_hmp(c):
    # test lca summarize --with-abundance on some real data
    queryfile = utils.get_test_data('hmp-sigs/G36354.sig.gz')
    dbname = utils.get_test_data('hmp-sigs/G36354-matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile)

    assert '32.1%  1080   p__Firmicutes;c__Bacilli;o__Lactobacillales' in c.last_result.out


@utils.in_thisdir
def test_summarize_abund_fake_no_abund(c):
    # test lca summarize on some known/fake data; see docs for explanation.
    queryfile = utils.get_test_data('fake-abund/query.sig.gz')
    dbname = utils.get_test_data('fake-abund/matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile,
                   '--ignore-abundance')

    assert 'NOTE: discarding abundances in query, since --ignore-abundance' in c.last_result.err
    assert '79.6%   550   Bacteria' in c.last_result.out
    assert '20.4%   141   Archaea' in c.last_result.out


@utils.in_thisdir
def test_summarize_abund_fake_yes_abund(c):
    # test lca summarize abundance weighting on some known/fake data
    queryfile = utils.get_test_data('fake-abund/query.sig.gz')
    dbname = utils.get_test_data('fake-abund/matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile)

    assert '43.2%   563   Bacteria' in c.last_result.out
    assert '56.8%   740   Archaea' in c.last_result.out


def test_rankinfo_on_multi(runtmp):
    db1 = utils.get_test_data('lca/dir1.lca.json')
    db2 = utils.get_test_data('lca/dir2.lca.json')

    cmd = ['lca', 'rankinfo', db1, db2]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    lines = runtmp.last_result.out.splitlines()
    lines.remove('superkingdom: 0 (0.0%)')
    lines.remove('phylum: 464 (12.8%)')
    lines.remove('class: 533 (14.7%)')
    lines.remove('order: 1050 (29.0%)')
    lines.remove('family: 695 (19.2%)')
    lines.remove('genus: 681 (18.8%)')
    lines.remove('species: 200 (5.5%)')
    lines.remove('strain: 0 (0.0%)')

    assert not lines


def test_rankinfo_on_single(runtmp):
    db1 = utils.get_test_data('lca/both.lca.json')

    cmd = ['lca', 'rankinfo', db1]
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    lines = runtmp.last_result.out.splitlines()
    lines.remove('superkingdom: 0 (0.0%)')
    lines.remove('phylum: 464 (12.8%)')
    lines.remove('class: 533 (14.7%)')
    lines.remove('order: 1050 (29.0%)')
    lines.remove('family: 695 (19.2%)')
    lines.remove('genus: 681 (18.8%)')
    lines.remove('species: 200 (5.5%)')
    lines.remove('strain: 0 (0.0%)')

    assert not lines


def test_rankinfo_no_tax(runtmp, lca_db_format):
    # note: TARA_PSW_MAG_00136 is _not_ in delmont-1.csv.
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = runtmp.output(f'delmont-1.lca.{lca_db_format}')

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-F', lca_db_format]
    runtmp.sourmash(*cmd)

    print('cmd:', cmd)
    print('out:', runtmp.last_result.out)
    print('err:', runtmp.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in runtmp.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in runtmp.last_result.err
    assert '0 identifiers used out of 1 distinct identifiers in spreadsheet.' in runtmp.last_result.err

    cmd = ['lca', 'rankinfo', lca_db]
    runtmp.sourmash(*cmd)


def test_rankinfo_with_min(runtmp):
    db1 = utils.get_test_data('lca/dir1.lca.json')
    db2 = utils.get_test_data('lca/dir2.lca.json')

    cmd = ['lca', 'rankinfo', db1, db2, '--minimum-num', '1']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    lines = runtmp.last_result.out.splitlines()
    lines.remove('superkingdom: 0 (0.0%)')
    lines.remove('phylum: 464 (12.8%)')
    lines.remove('class: 533 (14.7%)')
    lines.remove('order: 1050 (29.0%)')
    lines.remove('family: 695 (19.2%)')
    lines.remove('genus: 681 (18.8%)')
    lines.remove('species: 200 (5.5%)')
    lines.remove('strain: 0 (0.0%)')

    assert not lines


def test_rankinfo_with_min_2(runtmp):
    db1 = utils.get_test_data('lca/dir1.lca.json')
    db2 = utils.get_test_data('lca/dir2.lca.json')

    cmd = ['lca', 'rankinfo', db1, db2, '--minimum-num', '2']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "(no hashvals with lineages found)" in runtmp.last_result.err


def test_compare_csv(runtmp):
    a = utils.get_test_data('lca/classify-by-both.csv')
    b = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')

    cmd = ['lca', 'compare_csv', a, b, '-f']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 106 distinct lineages, 957 rows' in runtmp.last_result.err
    assert 'missing 937 assignments in classify spreadsheet.' in runtmp.last_result.err
    assert '20 total assignments, 0 differ between spreadsheets.' in runtmp.last_result.err


def test_compare_csv_real(runtmp):
    a = utils.get_test_data('lca/tully-genome-sigs.classify.csv')
    b = utils.get_test_data('lca/tully-query.delmont-db.sigs.classify.csv')

    cmd = ['lca', 'compare_csv', a, b, '--start-column=3', '-f']
    runtmp.sourmash(*cmd)

    print(cmd)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 87 distinct lineages, 2631 rows' in runtmp.last_result.err
    assert 'missing 71 assignments in classify spreadsheet.' in runtmp.last_result.err
    assert 'missing 1380 assignments in custom spreadsheet.' in runtmp.last_result.err
    assert '(these will not be evaluated any further)' in runtmp.last_result.err
    assert '987 total assignments, 889 differ between spreadsheets.' in runtmp.last_result.err
    assert '296 are compatible (one lineage is ancestor of another.' in runtmp.last_result.err
    assert '593 are incompatible (there is a disagreement in the trees).' in runtmp.last_result.err
    assert '164 incompatible at rank superkingdom' in runtmp.last_result.err
    assert '255 incompatible at rank phylum' in runtmp.last_result.err
    assert '107 incompatible at rank class' in runtmp.last_result.err
    assert '54 incompatible at rank order' in runtmp.last_result.err
    assert '13 incompatible at rank family' in runtmp.last_result.err
    assert '0 incompatible at rank genus' in runtmp.last_result.err
    assert '0 incompatible at rank species' in runtmp.last_result.err


def test_incompat_lca_db_ksize_2_fail(runtmp, lca_db_format):
    # test on gather - create a database with ksize of 25 => fail
    # because of incompatibility.
    c = runtmp
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('sketch', 'dna', '-p', 'k=25,scaled=1000', testdata1,
                   '-o', 'test_db.sig')
    print(c)

    c.run_sourmash('lca', 'index', utils.get_test_data('lca/delmont-1.csv',),
                   f'test.lca.{lca_db_format}', 'test_db.sig',
                    '-k', '25', '--scaled', '10000',
                   '-F', lca_db_format)
    print(c)

    # this should fail: the LCA database has ksize 25, and the query sig has
    # no compatible ksizes.
    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('gather', utils.get_test_data('lca/TARA_ASE_MAG_00031.sig'), f'test.lca.{lca_db_format}')

    err = c.last_result.err
    print(err)

    if lca_db_format == 'sql':
        assert "no compatible signatures found in 'test.lca.sql'" in err
    else:
        assert "ERROR: cannot use 'test.lca.json' for this query." in err
        assert "ksize on this database is 25; this is different from requested ksize of 31"


def test_incompat_lca_db_ksize_2_nofail(runtmp, lca_db_format):
    # test on gather - create a database with ksize of 25, no fail
    # because of --no-fail-on-empty-databases
    c = runtmp
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('sketch', 'dna', '-p', 'k=25,scaled=1000', testdata1,
                   '-o', 'test_db.sig')
    print(c)

    c.run_sourmash('lca', 'index', utils.get_test_data('lca/delmont-1.csv',),
                   f'test.lca.{lca_db_format}', 'test_db.sig',
                    '-k', '25', '--scaled', '10000',
                   '-F', lca_db_format)
    print(c)

    # this should not fail despite mismatched ksize, b/c of --no-fail flag.
    c.run_sourmash('gather', utils.get_test_data('lca/TARA_ASE_MAG_00031.sig'), f'test.lca.{lca_db_format}', '--no-fail-on-empty-database')

    err = c.last_result.err
    print(err)

    if lca_db_format == 'sql':
        assert "no compatible signatures found in 'test.lca.sql'" in err
    else:
        assert "ERROR: cannot use 'test.lca.json' for this query." in err
        assert "ksize on this database is 25; this is different from requested ksize of 31"


def test_lca_index_empty(runtmp, lca_db_format):
    c = runtmp
    # test lca index with an empty taxonomy CSV, followed by a load & gather.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig63 = load_one_signature(sig63file, ksize=31)

    # create an empty spreadsheet
    with open(c.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    # index!
    c.run_sourmash('lca', 'index', 'empty.csv', 'xxx',
                   sig2file, sig47file, sig63file, '--scaled', '1000',
                   '-F', lca_db_format)

    # can we load and search?
    lca_db_filename = c.output(f'xxx.lca.{lca_db_format}')
    db, ksize, scaled = lca_utils.load_single_database(lca_db_filename)

    result = db.best_containment(sig63)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == sig63.minhash
    assert name == lca_db_filename


def test_lca_gather_threshold_1():
    # test gather() method, in some detail; see same tests for sbt.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        db.best_containment(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    result = db.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        db.best_containment(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    result = db.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        db.best_containment(SourmashSignature(new_mh), threshold_bp=5000)


def test_lca_gather_threshold_5():
    # test gather() method, in some detail; see same tests for sbt.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures both have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # add five hashes
    for i in range(5):
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())

    # should get a result with no threshold (any match at all is returned)
    result = db.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # now, check with a threshold_bp that should be meet-able.
    result = db.best_containment(SourmashSignature(new_mh), threshold_bp=5000)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None


def test_gather_multiple_return():
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now, run gather. how many results do we get, and are they in the
    # right order?
    result = db.best_containment(sig63)
    print(result)
    assert result
    assert result.score == 1.0


def test_lca_db_protein_build():
    # test programmatic creation of LCA database with protein sigs in it
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='protein')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db.best_containment(sig2)
    assert result.score == 1.0


@utils.in_tempdir
def test_lca_db_protein_save_load(c):
    # test save/load of programmatically created db with protein sigs in it
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='protein')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'protein'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    print('XXX', mh_list[0].ksize)
    print('YYY', sig1.minhash.ksize)
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


def test_lca_db_protein_command_index(runtmp, lca_db_format):
    # test command-line creation of LCA database with protein sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output(f'protein.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '2', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--protein',
                   '-F', lca_db_format)

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'protein'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


@utils.in_thisdir
def test_lca_db_protein_command_search(c):
    # test command-line search/gather of LCA database with protein sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_lca_db_hp_build():
    # test programmatic creation of LCA database with hp sigs in it
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='hp')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db.best_containment(sig2)
    assert result.score == 1.0


@utils.in_tempdir
def test_lca_db_hp_save_load(c):
    # test save/load of programmatically created db with hp sigs in it
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='hp')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'hp'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


def test_lca_db_hp_command_index(runtmp, lca_db_format):
    # test command-line creation of LCA database with hp sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output(f'hp.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '2', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--hp',
                   '-F', lca_db_format)

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'hp'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


@utils.in_thisdir
def test_lca_db_hp_command_search(c):
    # test command-line search/gather of LCA database with hp sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/hp.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_lca_db_dayhoff_build():
    # test programmatic creation of LCA database with dayhoff sigs in it
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='dayhoff')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db.best_containment(sig2)
    assert result.score == 1.0


@utils.in_tempdir
def test_lca_db_dayhoff_save_load(c):
    # test save/load of programmatically created db with dayhoff sigs in it
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='dayhoff')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'dayhoff'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


def test_lca_db_dayhoff_command_index(runtmp, lca_db_format):
    # test command-line creation of LCA database with dayhoff sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output(f'dayhoff.lca.{lca_db_format}')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '2', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--dayhoff',
                   '-F', lca_db_format)

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'dayhoff'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0


@utils.in_thisdir
def test_lca_db_dayhoff_command_search(c):
    # test command-line search/gather of LCA database with dayhoff sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/dayhoff.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_lca_index_with_picklist(runtmp, lca_db_format):
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    outdb = runtmp.output(f'gcf.lca.{lca_db_format}')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    # create an empty spreadsheet
    with open(runtmp.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    runtmp.sourmash('lca', 'index', 'empty.csv', outdb, *gcf_sigs,
                    '-k', '21', '--picklist', f"{picklist}:md5:md5",
                    '-F', lca_db_format)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    assert "for given picklist, found 3 matches to 9 distinct values" in err
    assert "WARNING: 6 missing picklist values."
    assert "WARNING: no lineage provided for 3 signatures" in err

    siglist = list(sourmash.load_file_as_signatures(outdb))
    assert len(siglist) == 3
    for ss in siglist:
        assert 'Thermotoga' in ss.name


def test_lca_index_with_picklist_exclude(runtmp, lca_db_format):
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    outdb = runtmp.output(f'gcf.lca.{lca_db_format}')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    # create an empty spreadsheet
    with open(runtmp.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    runtmp.sourmash('lca', 'index', 'empty.csv', outdb, *gcf_sigs,
                    '-k', '21', '--picklist', f"{picklist}:md5:md5:exclude",
                    '-F', lca_db_format)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    siglist = list(sourmash.load_file_as_signatures(outdb))
    assert len(siglist) == 9
    for ss in siglist:
        assert 'Thermotoga' not in ss.name


def test_lca_index_select_with_picklist(runtmp, lca_db_format):
    # check what happens with picklists after index
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    outdb = runtmp.output(f'gcf.lca.{lca_db_format}')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    # create an empty spreadsheet
    with open(runtmp.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    runtmp.sourmash('lca', 'index', 'empty.csv', outdb, *gcf_sigs,
                    '-k', '21', '-F', lca_db_format)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    idx = sourmash.load_file_as_index(outdb)
    picklist_obj = SignaturePicklist.from_picklist_args(f"{picklist}:md5:md5")
    picklist_obj.load(picklist_obj.pickfile, picklist_obj.column_name)

    idx = idx.select(picklist=picklist_obj)

    siglist = list(idx.signatures())
    assert len(siglist) == 3
    for ss in siglist:
        assert 'Thermotoga' in ss.name


def test_lca_index_select_with_picklist_exclude(runtmp, lca_db_format):
    # check what happens with picklists after index
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    outdb = runtmp.output(f'gcf.lca.{lca_db_format}')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    # create an empty spreadsheet
    with open(runtmp.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    runtmp.sourmash('lca', 'index', 'empty.csv', outdb, *gcf_sigs,
                    '-k', '21', '-F', lca_db_format)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    idx = sourmash.load_file_as_index(outdb)
    picklist_obj = SignaturePicklist.from_picklist_args(f"{picklist}:md5:md5:exclude")
    picklist_obj.load(picklist_obj.pickfile, picklist_obj.column_name)
    idx = idx.select(picklist=picklist_obj)

    siglist = list(idx.signatures())
    assert len(siglist) == 9
    for ss in siglist:
        assert 'Thermotoga' not in ss.name


def test_lca_jaccard_ordering():
    # this tests a tricky situation where for three sketches A, B, C,
    # |A intersect B| is greater than |A intersect C|
    # _but_
    # |A jaccard B| is less than |A intersect B|
    a = sourmash.MinHash(ksize=31, n=0, scaled=2)
    b = a.copy_and_clear()
    c = a.copy_and_clear()

    a.add_many([1, 2, 3, 4])
    b.add_many([1, 2, 3] + list(range(10, 30)))
    c.add_many([1, 5])

    def _intersect(x, y):
        return x.intersection_and_union_size(y)[0]

    print('a intersect b:', _intersect(a, b))
    print('a intersect c:', _intersect(a, c))
    print('a jaccard b:', a.jaccard(b))
    print('a jaccard c:', a.jaccard(c))
    assert _intersect(a, b) > _intersect(a, c)
    assert a.jaccard(b) < a.jaccard(c)

    # thresholds to use:
    assert a.jaccard(b) < 0.15
    assert a.jaccard(c) > 0.15

    # now - make signatures, try out :)
    ss_a = sourmash.SourmashSignature(a, name='A')
    ss_b = sourmash.SourmashSignature(b, name='B')
    ss_c = sourmash.SourmashSignature(c, name='C')

    db = sourmash.lca.LCA_Database(ksize=31, scaled=2)
    db.insert(ss_a)
    db.insert(ss_b)
    db.insert(ss_c)

    sr = db.search(ss_a, threshold=0.15)
    print(sr)
    assert len(sr) == 2
    assert sr[0].signature == ss_a
    assert sr[0].score == 1.0
    assert sr[1].signature == ss_c
    assert sr[1].score == 0.2


def test_lca_db_protein_save_twice(runtmp, lca_db_format):
    # test save twice
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='protein')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(runtmp.output('xxx'), format=lca_db_format)

    with pytest.raises(ValueError):
        db.save(runtmp.output('xxx'), format=lca_db_format)
