"Tests for SqliteIndex, SqliteCollectionManifest, and LCA_SqliteDatabase"
import os
import pytest
import shutil
import sqlite3

import sourmash
from sourmash.exceptions import IndexNotSupported
from sourmash.index.sqlite_index import (SqliteIndex, load_sqlite_index,
                                         SqliteCollectionManifest,
                                         LCA_SqliteDatabase)

from sourmash.index import StandaloneManifestIndex
from sourmash import load_one_signature, SourmashSignature
from sourmash.picklist import SignaturePicklist, PickStyle
from sourmash.manifest import CollectionManifest
from sourmash.tax.tax_utils import MultiLineageDB

import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed
from sourmash import sqlite_utils


def test_sqlite_index_prefetch_empty():
    # check that an exception is raised upon for an empty database
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    sqlidx = SqliteIndex.create(":memory:")

    # since this is a generator, we need to actually ask for a value to
    # get exception raised.
    g = sqlidx.prefetch(ss2, threshold_bp=0)
    with pytest.raises(ValueError) as e:
        next(g)

    assert "no signatures to search" in str(e.value)


def test_sqlite_index_bad_version(runtmp):
    # create a sqlite database with a bad index version in the
    # sourmash_internal table, see what happens :)

    dbfile = runtmp.output('xyz.sqldb')
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()

    SqliteIndex._create_tables(c)

    # 0.9 doesn't exist/is bad version
    c.execute('UPDATE sourmash_internal SET value=? WHERE key=?',
              ('0.9', 'SqliteIndex'))

    conn.commit()

    with pytest.raises(IndexNotSupported):
        idx = sourmash.load_file_as_index(dbfile)


def test_sqlite_index_bad_version_unique(runtmp):
    # try to insert duplicate sqlite index info into sourmash_internal; fail

    dbfile = runtmp.output('xyz.sqldb')
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()

    SqliteIndex._create_tables(c)

    # can't insert duplicate key
    with pytest.raises(sqlite3.IntegrityError):
        c.execute('INSERT INTO sourmash_internal (value, key) VALUES (?, ?)',
                  ('1.1', 'SqliteIndex'))


def test_index_search_subj_scaled_is_lower():
    # check that subject sketches are appropriately downsampled
    sigfile = utils.get_test_data('scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz')
    ss = sourmash.load_one_signature(sigfile)

    # double check :)
    assert ss.minhash.scaled == 100

    # build a new query that has a scaled of 1000
    qs = SourmashSignature(ss.minhash.downsample(scaled=1000))

    # create Index to search
    sqlidx = SqliteIndex.create(":memory:")
    sqlidx.insert(ss)

    # search!
    results = list(sqlidx.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_sqlite_index_save_load(runtmp):
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    filename = runtmp.output('foo')
    sqlidx = SqliteIndex.create(filename)
    sqlidx.insert(ss2)
    sqlidx.insert(ss47)
    sqlidx.insert(ss63)

    sqlidx.close()

    sqlidx2 = SqliteIndex.load(filename)

    # now, search for sig2
    sr = sqlidx2.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2


def test_sqlite_index_multik_select():
    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex.create(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # select most specifically
    sqlidx2 = sqlidx.select(ksize=31, moltype='DNA')
    assert len(sqlidx2) == 1

    # all are DNA:
    sqlidx2 = sqlidx.select(moltype='DNA')
    assert len(sqlidx2) == 3


def test_sqlite_index_num_select():
    # this will fail on 'num' select, which is not allowed
    sqlidx = SqliteIndex.create(":memory:")
    with pytest.raises(ValueError):
        sqlidx.select(num=100)


def test_sqlite_index_abund_select():
    # this will fail on 'track_abundance' select, which is not allowed
    sqlidx = SqliteIndex.create(":memory:")
    with pytest.raises(ValueError):
        sqlidx.select(track_abundance=True)


def test_sqlite_index_insert_num_fail():
    # cannot insert 'num' signatures
    sqlidx = SqliteIndex.create(":memory:")

    sig47 = utils.get_test_data('num/47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    assert ss47.minhash.num != 0

    with pytest.raises(ValueError) as exc:
        sqlidx.insert(ss47)

    assert "cannot store 'num' signatures in SqliteIndex" in str(exc)


def test_sqlite_index_insert_abund_fail():
    # cannot insert 'num' signatures
    sqlidx = SqliteIndex.create(":memory:")

    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    ss47 = sourmash.load_one_signature(sig47, ksize=31)

    with pytest.raises(ValueError) as exc:
        sqlidx.insert(ss47)

    assert "cannot store signatures with abundance in SqliteIndex" in str(exc)


def test_sqlite_index_moltype_multi_fail():
    # check that we cannot store sigs with multiple scaled values.

    # this loads multiple ksizes (19, 31) and moltypes (DNA, protein, hp, etc)
    filename = utils.get_test_data('prot/all.zip')
    siglist = sourmash.load_file_as_signatures(filename)
    siglist = list(siglist)

    sqlidx = SqliteIndex.create(":memory:")

    sqlidx.insert(siglist[0])
    assert sqlidx.scaled == 100

    with pytest.raises(ValueError) as exc:
        for ss in siglist:
            sqlidx.insert(ss)

    assert "this database can only store scaled values=100" in str(exc)


def test_sqlite_index_picklist_select():
    # test select with a picklist

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex.create(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['f3a90d4e'])

    # select on picklist
    sqlidx2 = sqlidx.select(picklist=picklist)
    assert len(sqlidx2) == 1
    ss = list(sqlidx2.signatures())[0]
    assert ss.minhash.ksize == 31
    assert ss.md5sum().startswith('f3a90d4e55')


def test_sqlite_index_picklist_select_exclude():
    # test select with a picklist, but exclude

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    sqlidx = SqliteIndex.create(":memory:")
    for ss in siglist:
        sqlidx.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8', pickstyle=PickStyle.EXCLUDE)
    picklist.init(['f3a90d4e'])

    # select on picklist
    sqlidx2 = sqlidx.select(picklist=picklist)
    assert len(sqlidx2) == 2
    md5s = set()
    ksizes = set()
    for ss in list(sqlidx2.signatures()):
        md5s.add(ss.md5sum())
        ksizes.add(ss.minhash.ksize)
    assert md5s == set(['f372e47893edd349e5956f8b0d8dcbf7','43f3b48e59443092850964d355a20ac0'])
    assert ksizes == set([21,51])


def test_sqlite_jaccard_ordering():
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

    sqlidx = SqliteIndex.create(":memory:")
    sqlidx.insert(ss_a)
    sqlidx.insert(ss_b)
    sqlidx.insert(ss_c)

    sr = sqlidx.search(ss_a, threshold=0.15)
    print(sr)
    assert len(sr) == 2
    assert sr[0].signature == ss_a
    assert sr[0].score == 1.0
    assert sr[1].signature == ss_c
    assert sr[1].score == 0.2


def test_sqlite_index_scaled1():
    # check on scaled=1 storage.
    sqlidx = SqliteIndex.create(":memory:")

    mh1 = sourmash.MinHash(0, 31, scaled=1)
    mh1.add_hash(2**64 - 1)
    mh1.add_hash(2**64 - 2)
    mh1.add_hash(2**64 - 3)
    ss1 = sourmash.SourmashSignature(mh1, name='ss 1')

    mh2 = sourmash.MinHash(0, 31, scaled=1)
    mh2.add_hash(2**64 - 1)
    mh2.add_hash(2**64 - 2)
    mh2.add_hash(2**64 - 3)
    mh2.add_hash(0)
    mh2.add_hash(1)
    mh2.add_hash(2)
    ss2 = sourmash.SourmashSignature(mh2, name='ss 2')

    sqlidx.insert(ss1)
    sqlidx.insert(ss2)

    # check jaccard search
    results = list(sqlidx.search(ss1, threshold=0))
    print(results)
    assert len(results) == 2
    assert results[0].signature == ss1
    assert results[0].score == 1.0
    assert results[1].signature == ss2
    assert results[1].score == 0.5

    results = list(sqlidx.search(ss1, threshold=0, do_containment=True))
    print(results)
    assert results[0].signature == ss1
    assert results[0].score == 1.0
    assert results[1].signature == ss2
    assert results[1].score == 1.0

    # minhashes retrieved successfully?
    assert len(results[0].signature.minhash) == 3
    assert len(results[1].signature.minhash) == 6


def test_sqlite_index_load_existing():
    # try loading an existing sqlite index
    filename = utils.get_test_data('sqlite/index.sqldb')
    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, SqliteIndex)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2


def test_sqlite_index_create_load_existing(runtmp):
    # try creating then loading an existing sqlite index; create from CLI
    filename = runtmp.output('idx.sqldb')
    sig1 = utils.get_test_data('47.fa.sig')
    sig2 = utils.get_test_data('63.fa.sig')

    runtmp.sourmash('sig', 'cat', sig1, sig2, '-o', filename)

    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, SqliteIndex)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2


def test_sqlite_index_create_load_insert_existing(runtmp):
    # try creating, loading, inserting into an existing sqlite index
    filename = runtmp.output('idx.sqldb')
    sig1 = utils.get_test_data('47.fa.sig')
    sig2 = utils.get_test_data('63.fa.sig')
    sig3 = utils.get_test_data('2.fa.sig')

    runtmp.sourmash('sig', 'cat', sig1, sig2, '-o', filename)

    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, SqliteIndex)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2

    ss3 = sourmash.load_one_signature(sig3, ksize=31)
    sqlidx.insert(ss3)
    sqlidx.commit()

    runtmp.sourmash('sig', 'describe', filename)
    print(runtmp.last_result.out)
    assert "md5: f3a90d4e5528864a5bcc8434b0d0c3b1" in runtmp.last_result.out


def test_sqlite_index_create_load_insert_existing_cli(runtmp):
    # try creating, loading, inserting into an existing sqlite index from cli
    # (aka "append" to existing database)
    filename = runtmp.output('idx.sqldb')
    sig1 = utils.get_test_data('47.fa.sig')
    sig2 = utils.get_test_data('63.fa.sig')
    sig3 = utils.get_test_data('2.fa.sig')

    runtmp.sourmash('sig', 'cat', sig1, sig2, '-o', filename)

    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, SqliteIndex)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2

    # add a third
    runtmp.sourmash('sig', 'cat', sig3, '-o', filename, '-k', '31')

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 3


def test_sqlite_manifest_bad_version(runtmp):
    # create a sqlite database with a bad manifest version in the
    # sourmash_internal table, see what happens :)

    dbfile = runtmp.output('xyz.sqlmf')
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()

    SqliteCollectionManifest._create_tables(c)

    # 0.9 doesn't exist/bad version
    c.execute('UPDATE sourmash_internal SET value=? WHERE key=?',
              ('0.9', 'SqliteManifest'))

    conn.commit()

    with pytest.raises(IndexNotSupported):
        mf = CollectionManifest.load_from_filename(dbfile)


def test_sqlite_manifest_bad_version_unique(runtmp):
    # try to insert duplicate sqlite manifest info into sourmash_internal; fail

    dbfile = runtmp.output('xyz.sqldb')
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()

    SqliteCollectionManifest._create_tables(c)

    # can't insert duplicate key
    with pytest.raises(sqlite3.IntegrityError):
        c.execute('INSERT INTO sourmash_internal (value, key) VALUES (?, ?)',
                  ('1.1', 'SqliteManifest'))


def test_sqlite_manifest_basic():
    # test some features of the SQLite-based manifest.
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    sqlidx = SqliteIndex.create(":memory:")

    # empty manifest tests
    manifest = sqlidx.manifest
    assert not manifest
    assert len(manifest) == 0

    sqlidx.insert(sig47)
    sqlidx.insert(sig63)

    # ok, more full manifest tests!
    assert manifest
    assert len(manifest) == 2

    assert sig47 in manifest
    assert sig2 not in manifest

    # check that we can get a "standard" manifest out
    standard_mf = CollectionManifest.load_from_manifest(sqlidx.manifest)
    assert len(standard_mf) == 2

    picklist = manifest.to_picklist()
    assert sig47 in picklist
    assert sig2 not in picklist


def test_sqlite_manifest_round_trip():
    # check that we can go from regular mf -> sqlite mf -> regular again.
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    rows = []
    rows.append(CollectionManifest.make_manifest_row(sig47, None,
                                                     include_signature=False))
    rows.append(CollectionManifest.make_manifest_row(sig63, None,
                                                     include_signature=False))
    nosql_mf = CollectionManifest(rows)

    sqlite_mf = SqliteCollectionManifest.load_from_manifest(nosql_mf)

    # test roundtrip
    round_mf = CollectionManifest.load_from_manifest(sqlite_mf)

    assert len(round_mf) == 2
    print(round_mf.rows, nosql_mf.rows)
    assert round_mf == nosql_mf

    for mf in (nosql_mf, sqlite_mf, round_mf):
        picklist = mf.to_picklist()
        assert sig47 in picklist
        assert sig2 not in picklist


def test_sqlite_manifest_create(runtmp):
    # test creation and summarization of a manifest of prot.zip
    zipfile = utils.get_test_data('prot/all.zip')

    # create manifest
    runtmp.sourmash('sig', 'manifest', '-F', 'sql', zipfile,
                    '-o', 'mf.sqlmf')

    sqlmf = runtmp.output('mf.sqlmf')
    assert os.path.exists(sqlmf)

    # verify it's loadable as the right type
    idx = load_sqlite_index(sqlmf)
    assert isinstance(idx, StandaloneManifestIndex)

    # summarize
    runtmp.sourmash('sig', 'fileinfo', 'mf.sqlmf')

    out = runtmp.last_result.out
    print(out)

    assert "2 sketches with dayhoff, k=19, scaled=100          7945 total hashes" in out
    assert "2 sketches with hp, k=19, scaled=100               5184 total hashes" in out
    assert "2 sketches with protein, k=19, scaled=100          8214 total hashes" in out
    assert "1 sketches with DNA, k=31, scaled=1000             5238 total hashes" in out

    assert "path filetype: StandaloneManifestIndex" in out
    assert "location: mf.sqlmf" in out
    assert "is database? yes" in out
    assert "has manifest? yes" in out
    assert "num signatures: 7" in out


def test_sqlite_manifest_create_noload_sigs(runtmp):
    # sigs should not be loadable from manifest this way...
    zipfile = utils.get_test_data('prot/all.zip')

    # create manifest
    runtmp.sourmash('sig', 'manifest', '-F', 'sql', zipfile,
                    '-o', 'mf.sqlmf')

    # 'describe' should not be able to load the sqlmf b/c prefix is wrong
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sig', 'describe', 'mf.sqlmf')


def test_sqlite_manifest_create_yesload_sigs(runtmp):
    # should be able to load after copying files
    zipfile = utils.get_test_data('prot/all.zip')
    shutil.copytree(utils.get_test_data('prot'), runtmp.output('prot'))

    # create manifest
    runtmp.sourmash('sig', 'manifest', '-F', 'sql', zipfile,
                    '-o', 'prot/mf.sqlmf')

    # 'describe' should now be able to load the sqlmf, which is cool
    runtmp.sourmash('sig', 'describe', 'prot/mf.sqlmf')
    print(runtmp.last_result.out)


def test_sqlite_manifest_num(runtmp):
    # should be able to produce sql manifests with 'num' sketches in them
    numsig = utils.get_test_data('num/47.fa.sig')

    # create mf
    runtmp.sourmash('sig', 'manifest', '-F', 'sql', numsig,
                    '-o', 'mf.sqlmf')

    # do summarize:
    runtmp.sourmash('sig', 'summarize', 'mf.sqlmf')
    out = runtmp.last_result.out

    print(out)

    assert "1 sketches with DNA, k=21, num=500                 500 total hashes" in out
    assert "1 sketches with DNA, k=31, num=500                 500 total hashes" in out
    assert "1 sketches with DNA, k=51, num=500                 500 total hashes" in out


def test_sqlite_manifest_num_select(runtmp):
    # should be able to _select_ sql manifests with 'num' sketches in them
    numsig = utils.get_test_data('num/47.fa.sig')

    # create mf
    runtmp.sourmash('sig', 'manifest', '-F', 'sql', numsig,
                    '-o', 'mf.sqlmf')

    # load as index
    idx = sourmash.load_file_as_index(runtmp.output('mf.sqlmf'))

    # select
    print(list(idx.manifest.rows))
    idx = idx.select(num=500)
    print(list(idx.manifest.rows))
    assert len(idx) == 3


def test_sqlite_manifest_locations(runtmp):
    # check what locations returns... may return too many, that's ok.
    prot = utils.get_test_data('prot')

    runtmp.sourmash('sig', 'manifest', '-F', 'sql', prot,
                    '-o', 'mf.sqlmf')

    # load as index
    idx = sourmash.load_file_as_index(runtmp.output('mf.sqlmf'))

    picklist = SignaturePicklist('identprefix')
    picklist.pickset = set(['GCA_001593925'])
    idx = idx.select(picklist=picklist)

    sql_locations = set(idx.manifest.locations())
    row_locations = set(row['internal_location'] for row in idx.manifest.rows)

    assert sql_locations.issuperset(row_locations)

    assert 'dna-sig.sig.gz' in sql_locations # this is unnecessary...
    assert 'dna-sig.sig.gz' not in row_locations # ...this is correct :)


def test_sqlite_manifest_create_insert(runtmp):
    # try out creating a sqlite manifest and then running cli on it

    mfname = runtmp.output("some.sqlmf")
    mf = SqliteCollectionManifest.create(mfname)

    sigfile = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sigfile)

    mf._insert_row(mf.conn.cursor(), mf.make_manifest_row(ss, 'some.sig'))
    mf.conn.commit()

    # copy sig in since we want it to resolve...
    shutil.copyfile(sigfile, runtmp.output('some.sig'))

    # 'describe' should work here, to resolve actual sigs.
    runtmp.sourmash('sig', 'describe', mfname)
    print(runtmp.last_result.out)
    assert 'md5: 09a08691ce52952152f0e866a59f6261' in runtmp.last_result.out


def test_sqlite_manifest_create_insert_2(runtmp):
    # try out creating a sqlite manifest from cli and then _insert_row into it

    # copy sig in since we want it to resolve...
    sigfile = utils.get_test_data('47.fa.sig')
    shutil.copyfile(sigfile, runtmp.output('some.sig'))

    runtmp.sourmash('sig', 'manifest', 'some.sig', '-F', 'sql',
                    '-o', 'some.sqlmf')
    mfname = runtmp.output("some.sqlmf")

    mf = CollectionManifest.load_from_filename(mfname)
    ss = sourmash.load_one_signature(runtmp.output('some.sig'))
    mf._insert_row(mf.conn.cursor(), mf.make_manifest_row(ss, 'some.sig'))
    mf.conn.commit()

    # 'describe' should work here, to resolve actual sigs.
    runtmp.sourmash('sig', 'describe', mfname)
    print(runtmp.last_result.out)
    assert 'md5: 09a08691ce52952152f0e866a59f6261' in runtmp.last_result.out


def test_sqlite_manifest_existing(runtmp):
    # try out an existing sqlite manifest

    prefix = runtmp.output('protdir')
    mf = runtmp.output('protdir/prot.sqlmf')
    shutil.copytree(utils.get_test_data('prot'), prefix)
    shutil.copyfile(utils.get_test_data('sqlite/prot.sqlmf'), mf)

    runtmp.sourmash('sig', 'describe', mf)
    print(runtmp.last_result.out)


def test_sqlite_manifest_existing_insert(runtmp):
    # try out an existing sqlite manifest - insert into it

    prefix = runtmp.output('protdir')
    shutil.copytree(utils.get_test_data('prot'), prefix)

    mfname = runtmp.output('protdir/prot.sqlmf')
    shutil.copyfile(utils.get_test_data('sqlite/prot.sqlmf'), mfname)
    mf = CollectionManifest.load_from_filename(mfname)
    assert isinstance(mf, SqliteCollectionManifest)

    sigfile = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sigfile)

    mf._insert_row(mf.conn.cursor(), mf.make_manifest_row(ss, 'some.sig'))
    mf.conn.commit()

    # copy sig in since we want it to resolve...
    shutil.copyfile(sigfile, runtmp.output('protdir/some.sig'))

    # 'describe' should work here.
    runtmp.sourmash('sig', 'describe', mfname)
    print(runtmp.last_result.out)


def test_sqlite_manifest_existing_mf_only(runtmp):
    # try out an existing sqlite manifest, but without underlying files -> fail

    mf = runtmp.output('prot.sqlmf')
    shutil.copyfile(utils.get_test_data('sqlite/prot.sqlmf'), mf)

    # 'fileinfo' should work...
    runtmp.sourmash('sig', 'fileinfo', mf)
    print(runtmp.last_result.out)
    assert 'num signatures: 7' in runtmp.last_result.out

    # ...but 'describe' should fail, since it needs actual sigs.
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash('sig', 'describe', mf)

    print(runtmp.last_result.err)
    assert 'ERROR: Error while reading signatures from' in runtmp.last_result.err


def test_sqlite_manifest_existing_mfonly_insert(runtmp):
    # try out an existing sqlite manifest - insert into it, but fail describe

    mfname = runtmp.output('prot.sqlmf')
    shutil.copyfile(utils.get_test_data('sqlite/prot.sqlmf'), mfname)
    mf = CollectionManifest.load_from_filename(mfname)
    assert isinstance(mf, SqliteCollectionManifest)

    sigfile = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sigfile)

    mf._insert_row(mf.conn.cursor(), mf.make_manifest_row(ss, sigfile))
    mf.conn.commit()

    # 'fileinfo' should work...
    runtmp.sourmash('sig', 'fileinfo', mfname)
    print(runtmp.last_result.out)
    assert 'num signatures: 8' in runtmp.last_result.out

    # ...but 'describe' should fail, since it needs actual sigs.
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash('sig', 'describe', mfname)


def test_sqlite_manifest_load_existing_index():
    # try loading an existing sqlite index as a manifest
    filename = utils.get_test_data('sqlite/index.sqldb')
    mf = CollectionManifest.load_from_filename(filename)
    assert isinstance(mf, SqliteCollectionManifest)

    assert len(mf) == 2


def test_sqlite_manifest_load_existing_index_insert_fail():
    # try loading an existing sqlite index as a manifest; insert should fail
    filename = utils.get_test_data('sqlite/index.sqldb')
    mf = CollectionManifest.load_from_filename(filename)
    assert isinstance(mf, SqliteCollectionManifest)

    assert len(mf) == 2

    # try insert - should fail
    sigfile = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_one_signature(sigfile)

    with pytest.raises(Exception) as exc:
        mf._insert_row(mf.conn.cursor(), mf.make_manifest_row(ss, sigfile))

    assert "must use SqliteIndex.insert to add to this manifest" in str(exc)


def test_sqlite_manifest_create_load_empty(runtmp):
    # try creating an empty manifest, then loading
    mfname = runtmp.output("some.sqlmf")
    mf = SqliteCollectionManifest.create(mfname)
    mf.close()

    mf2 = load_sqlite_index(mfname)
    assert len(mf2) == 0


def test_sqlite_lca_db_load_existing():
    # try loading an existing sqlite index
    filename = utils.get_test_data('sqlite/lca.sqldb')
    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, LCA_SqliteDatabase)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2


def test_sqlite_lca_db_select():
    # try loading an existing sqlite index
    filename = utils.get_test_data('sqlite/lca.sqldb')
    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, LCA_SqliteDatabase)

    sqlidx2 = sqlidx.select(ksize=31)
    x = list(sqlidx2.hashvals)  # only on LCA_SqliteDatabase
    assert isinstance(sqlidx2, LCA_SqliteDatabase)


def test_sqlite_lca_db_create_load_existing(runtmp):
    # try creating (from CLI) then loading (from API) an LCA db
    filename = runtmp.output('lca.sqldb')
    sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')

    runtmp.sourmash('sig', 'flatten', sig1, sig2, '-o', filename, '-k', '31')

    # load tax
    tax_csv = utils.get_test_data('sqlite/delmont-6.csv')
    runtmp.sourmash('tax', 'prepare', '-t', tax_csv,
                    '-o', filename, '-F', 'sql')

    sqlidx = sourmash.load_file_as_index(filename)
    assert isinstance(sqlidx, LCA_SqliteDatabase)

    siglist = list(sqlidx.signatures())
    assert len(siglist) == 2


def test_sqlite_lca_db_load_empty(runtmp):
    # try creating then loading an _empty_ LCA_SqliteDatabase

    dbname = runtmp.output('empty.sqldb')

    # create empty SqliteIndex...
    runtmp.sourmash('sig', 'cat', '-o', dbname)
    assert os.path.exists(dbname)

    # ...and create empty sourmash_taxonomy tables in there...
    empty_tax = utils.get_test_data('scaled/empty-lineage.csv')
    runtmp.sourmash('tax', 'prepare', '-F', 'sql', '-t', empty_tax,
                    '-o', dbname)

    runtmp.sourmash('sig', 'describe', dbname)
    assert 'loaded 0 signatures' in runtmp.last_result.err


def test_sqlite_lca_db_create_readonly(runtmp):
    # try running 'prepare' on a read-only sqlite db, check error message.

    dbname = runtmp.output('empty.sqldb')

    # create empty SqliteIndex...
    runtmp.sourmash('sig', 'cat', '-o', dbname)
    assert os.path.exists(dbname)

    # make it read only...
    from stat import S_IREAD, S_IRGRP, S_IROTH
    os.chmod(dbname, S_IREAD|S_IRGRP|S_IROTH)

    # ...and try creating empty sourmash_taxonomy tables in there...
    empty_tax = utils.get_test_data('scaled/empty-lineage.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash('tax', 'prepare', '-F', 'sql', '-t', empty_tax,
                        '-o', dbname)

    err = runtmp.last_result.err
    print(err)

    assert not "taxonomy table already exists in" in err
    assert "attempt to write a readonly database" in err


def test_sqlite_lca_db_try_load_sqlite_index():
    # try loading a SqliteIndex with no tax tables from .load classmethod
    dbname = utils.get_test_data('sqlite/index.sqldb')

    with pytest.raises(ValueError) as exc:
        db = LCA_SqliteDatabase.load(dbname)

    assert "not a taxonomy database" in str(exc)


def test_sqlite_lca_db_supply_lineage_db():
    # try creating an LCA_SqliteDatabase object with a separate lineage DB.
    dbname = utils.get_test_data('sqlite/index.sqldb')

    tax_csv = utils.get_test_data('sqlite/shewanella-lineage.csv')
    lineage_db = MultiLineageDB.load([tax_csv])

    db = LCA_SqliteDatabase(dbname, lineage_db=lineage_db)

    hashval = next(iter(db.hashvals))
    lineages = db.get_lineage_assignments(hashval)
    print(lineages)
    assert lineages[0][0].rank == 'superkingdom'
    assert lineages[0][0].name == 'd__Bacteria'
    assert lineages[0][-1].rank == 'species'
    assert lineages[0][-1].name == 's__Shewanella baltica'
    assert lineages[1][0].rank == 'superkingdom'
    assert lineages[1][0].name == 'd__Bacteria'
    assert lineages[0][-1].rank == 'species'
    assert lineages[0][-1].name == 's__Shewanella baltica'


def test_bad_sqlite_internal_version():
    # check get_sourmash_internal
    dbname = utils.get_test_data('sqlite/index.sqldb')

    conn = sqlite_utils.open_sqlite_db(dbname)
    c = conn.cursor()
    with pytest.raises(Exception):
        sqlite_utils.add_sourmash_internal(c, 'SqliteIndex', '0.9')
