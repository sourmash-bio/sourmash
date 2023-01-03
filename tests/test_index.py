"""
Tests for Index classes and subclasses.
"""
import pytest
import glob
import os
import zipfile
import shutil

import sourmash
from sourmash import load_one_signature, SourmashSignature
from sourmash.index import (LinearIndex, ZipFileLinearIndex,
                            make_jaccard_search_query, CounterGather,
                            LazyLinearIndex, MultiIndex,
                            StandaloneManifestIndex)
from sourmash.index.revindex import RevIndex
from sourmash.sbt import SBT, GraphFactory
from sourmash import sourmash_args
from sourmash.search import JaccardSearch, SearchType
from sourmash.picklist import SignaturePicklist, PickStyle
from sourmash_tst_utils import SourmashCommandFailed
from sourmash.manifest import CollectionManifest

import sourmash_tst_utils as utils


def test_simple_index(n_children):
    # test basic SBT functionality
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf1_mh.add_sequence("AAAAA")
    leaf1_mh.add_sequence("AAAAT")
    leaf1_mh.add_sequence("AAAAC")
    leaf1_sig = SourmashSignature(leaf1_mh)
    root.insert(leaf1_sig)

    leaf2_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf2_mh.add_sequence("AAAAA")
    leaf2_mh.add_sequence("AAAAT")
    leaf2_mh.add_sequence("AAAAG")
    leaf2_sig = SourmashSignature(leaf2_mh)
    root.insert(leaf2_sig)

    leaf3_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf3_mh.add_sequence("AAAAA")
    leaf3_mh.add_sequence("AAAAT")
    leaf3_mh.add_sequence("CAAAA")
    leaf3_sig = SourmashSignature(leaf3_mh)
    root.insert(leaf3_sig)

    leaf4_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf4_mh.add_sequence("AAAAA")
    leaf4_mh.add_sequence("CAAAA")
    leaf4_mh.add_sequence("GAAAA")
    leaf4_sig = SourmashSignature(leaf4_mh)
    root.insert(leaf4_sig)

    leaf5_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf5_mh.add_sequence("AAAAA")
    leaf5_mh.add_sequence("AAAAT")
    leaf5_mh.add_sequence("GAAAA")
    leaf5_sig = SourmashSignature(leaf5_mh)
    root.insert(leaf5_sig)

    linear = LinearIndex()
    linear.insert(leaf1_sig)
    linear.insert(leaf2_sig)
    linear.insert(leaf3_sig)
    linear.insert(leaf4_sig)
    linear.insert(leaf5_sig)

    search_fn = make_jaccard_search_query(do_containment=True)

    kmers = ["AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA"]
    for kmer in kmers:
        search_mh = sourmash.MinHash(0, 5, scaled=1)
        search_mh.add_sequence(kmer)
        search_sig = sourmash.SourmashSignature(search_mh)

        linear_found = linear.find(search_fn, search_sig)
        linear_found = set(linear_found)

        tree_found = set(root.find(search_fn, search_sig))

        assert tree_found
        assert tree_found == set(linear_found)


def test_linear_index_prefetch_empty():
    # check that an exception is raised upon for an empty LinearIndex
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    lidx = LinearIndex()

    # since this is a generator, we need to actually ask for a value to
    # get exception raised.
    g = lidx.prefetch(ss2, threshold_bp=0)
    with pytest.raises(ValueError) as e:
        next(g)

    assert "no signatures to search" in str(e.value)


def test_linear_index_prefetch_lazy():
    # make sure that prefetch doesn't touch values 'til requested.
    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)
    fake = FakeSignature()

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)
    lidx.insert(fake)

    g = lidx.prefetch(ss47, threshold_bp=0)

    # first value:
    sr = next(g)
    assert sr.signature == ss47

    # second value:
    sr = next(g)
    assert sr.signature == ss63

    # third value: raises exception!
    with pytest.raises(Exception) as e:
        next(g)

    assert "don't touch me!" in str(e.value)


def test_linear_index_search_subj_has_abundance():
    # check that search signatures in the index are flattened appropriately.
    queryfile = utils.get_test_data('47.fa.sig')
    subjfile = utils.get_test_data('track_abund/47.fa.sig')

    qs = sourmash.load_one_signature(queryfile)
    ss = sourmash.load_one_signature(subjfile)

    linear = LinearIndex()
    linear.insert(ss)

    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # note: search returns _original_ signature, not flattened
    assert results[0].signature == ss


def test_linear_index_gather_subj_has_abundance():
    # check that target signatures in the index are flattened appropriately.
    queryfile = utils.get_test_data('47.fa.sig')
    subjfile = utils.get_test_data('track_abund/47.fa.sig')

    qs = sourmash.load_one_signature(queryfile)
    ss = sourmash.load_one_signature(subjfile)

    linear = LinearIndex()
    linear.insert(ss)

    result = linear.best_containment(qs, threshold=0)
    assert result

    # note: gather returns _original_ signature, not flattened
    assert result.signature == ss


def test_index_search_subj_scaled_is_lower():
    # check that subject sketches are appropriately downsampled for scaled
    # sketches.
    sigfile = utils.get_test_data('scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz')
    ss = sourmash.load_one_signature(sigfile)

    # double check :)
    assert ss.minhash.scaled == 100

    # build a new query that has a scaled of 1000
    qs = SourmashSignature(ss.minhash.downsample(scaled=1000))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_index_search_subj_num_is_lower():
    # check that subject sketches are appropriately downsampled for num
    # sketches
    sigfile = utils.get_test_data('num/47.fa.sig')
    ss = sourmash.load_one_signature(sigfile, ksize=31)

    # double check :)
    assert ss.minhash.num == 500

    # build a new query that has a num of 250
    qs = SourmashSignature(ss.minhash.downsample(num=250))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_index_search_query_num_is_lower():
    # check that query sketches are appropriately downsampled for num.
    sigfile = utils.get_test_data('num/47.fa.sig')
    qs = sourmash.load_one_signature(sigfile, ksize=31)

    # double check :)
    assert qs.minhash.num == 500

    # build a new subject that has a num of 250
    ss = SourmashSignature(qs.minhash.downsample(num=250))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    assert results[0].signature == ss


def test_linear_index_search_abund():
    # test Index.search_abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    results = list(lidx.search_abund(ss47, threshold=0))
    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_linear_index_search_abund_downsample_query():
    # test Index.search_abund with query with higher scaled
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    # forcibly downsample ss47 for the purpose of this test :)
    ss47 = ss47.to_mutable()
    ss47.minhash = ss63.minhash.downsample(scaled=2000)
    assert ss63.minhash.scaled != ss47.minhash.scaled

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    results = list(lidx.search_abund(ss47, threshold=0))
    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_linear_index_search_abund_downsample_subj():
    # test Index.search_abund with subj with higher scaled
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    # forcibly downsample ss63 for the purpose of this test :)
    ss63 = ss63.to_mutable()
    ss63.minhash = ss63.minhash.downsample(scaled=2000)
    assert ss63.minhash.scaled != ss47.minhash.scaled

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    results = list(lidx.search_abund(ss47, threshold=0))
    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_linear_index_search_abund_requires_threshold():
    # test that Index.search_abund requires a 'threshold'
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=None))

    assert "'search_abund' requires 'threshold'" in str(exc.value)


def test_linear_index_search_abund_query_flat():
    # test that Index.search_abund requires an abund query sig
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=0))

    assert "'search_abund' requires query signature with abundance information" in str(exc.value)


def test_linear_index_search_abund_subj_flat():
    # test Index.search_abund requires an abund subj
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=0))

    assert "'search_abund' requires subject signatures with abundance information" in str(exc.value)


def test_linear_index_save(runtmp):
    # test save output from LinearIndex => JSON
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    linear = LinearIndex()
    linear.insert(ss2)
    linear.insert(ss47)
    linear.insert(ss63)

    filename = runtmp.output('foo')
    linear.save(filename)

    si = set(sourmash.load_file_as_signatures(filename))

    x = {ss2, ss47, ss63}

    print(len(si))
    print(len(x))

    print('si: ', si)
    print('x: ', x)

    assert si == x, si


def test_linear_index_load(runtmp):
    # test .load class method of LinearIndex
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    filename = runtmp.output('foo')
    with open(filename, 'wt') as fp:
        sourmash.save_signatures([ss2, ss47, ss63], fp)

    linear = LinearIndex.load(filename)

    x = {ss2, ss47, ss63}
    assert set(linear.signatures()) == x, linear.signatures
    assert linear.location == filename


def test_linear_index_save_load(runtmp):
    # LinearIndex save/load round trip
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    linear = LinearIndex()
    linear.insert(ss2)
    linear.insert(ss47)
    linear.insert(ss63)

    filename = runtmp.output('foo')
    linear.save(filename)
    linear2 = LinearIndex.load(filename)

    # now, search for sig2
    sr = linear2.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2


def test_linear_gather_threshold_1():
    # test gather() method, in some detail
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    linear = LinearIndex()

    linear.insert(sig47)
    linear.insert(sig63)
    linear.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        linear.best_containment(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    result = linear.best_containment(SourmashSignature(new_mh))
    assert result

    # it's a namedtuple, so we can unpack like a tuple.
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        linear.best_containment(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    result = linear.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        linear.best_containment(SourmashSignature(new_mh), threshold_bp=5000)


def test_linear_gather_threshold_5():
    # test gather() method above threshold
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    linear = LinearIndex(filename='foo')

    linear.insert(sig47)
    linear.insert(sig63)
    linear.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

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
    result = linear.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name == 'foo'

    # now, check with a threshold_bp that should be meet-able.
    result = linear.best_containment(SourmashSignature(new_mh),
                                     threshold_bp=5000)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name == 'foo'


def test_linear_index_multik_select():
    # test that LinearIndx can load multiple (three) ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # select most specifically
    linear2 = linear.select(ksize=31, moltype='DNA')
    assert len(linear2) == 1

    # all are DNA:
    linear2 = linear.select(moltype='DNA')
    assert len(linear2) == 3


def test_linear_index_moltype_select():
    # this loads two ksizes(21, 10), and two moltypes (DNA and protein)
    filename = utils.get_test_data('genome-s10+s11.sig')
    siglist = sourmash.load_file_as_signatures(filename)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # select most specific DNA
    linear2 = linear.select(ksize=30, moltype='DNA')
    assert len(linear2) == 1

    # select most specific protein
    linear2 = linear.select(ksize=10, moltype='protein')
    assert len(linear2) == 1

    # can leave off ksize, selects all ksizes
    linear2 = linear.select(moltype='DNA')
    assert len(linear2) == 2

    # can leave off ksize, selects all ksizes
    linear2 = linear.select(moltype='protein')
    assert len(linear2) == 2

    # select something impossible
    linear2 = linear.select(ksize=4)
    assert len(linear2) == 0


def test_linear_index_picklist_select():
    # test LinearIndex.select with a picklist

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8')
    picklist.init(['f3a90d4e'])

    # select on picklist
    linear2 = linear.select(picklist=picklist)
    assert len(linear2) == 1
    ss = list(linear2.signatures())[0]
    assert ss.minhash.ksize == 31
    assert ss.md5sum().startswith('f3a90d4e55')


def test_linear_index_picklist_select_exclude():
    # test select with a picklist, but exclude

    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # construct a picklist...
    picklist = SignaturePicklist('md5prefix8', pickstyle=PickStyle.EXCLUDE)
    picklist.init(['f3a90d4e'])

    # select on picklist
    linear2 = linear.select(picklist=picklist)
    assert len(linear2) == 2
    md5s = set()
    ksizes = set()
    for ss in list(linear2.signatures()):
        md5s.add(ss.md5sum())
        ksizes.add(ss.minhash.ksize)
    assert md5s == set(['f372e47893edd349e5956f8b0d8dcbf7','43f3b48e59443092850964d355a20ac0'])
    assert ksizes == set([21,51])


def test_index_same_md5sum_fsstorage(runtmp):
    # check SBT directory 'save' with two signatures that have identical md5
    c = runtmp
    testdata1 = utils.get_test_data('img/2706795855.sig')
    testdata2 = utils.get_test_data('img/638277004.sig')

    c.run_sourmash('index', '-k', '21', 'zzz.sbt.json', testdata1, testdata2)
    assert c.last_result.status == 0

    outfile = c.output('zzz.sbt.json')
    assert os.path.exists(outfile)
    storage = c.output('.sbt.zzz')
    assert len(glob.glob(storage + "/*")) == 4


def test_index_same_md5sum_sbt_zipstorage(runtmp):
    # check SBT zipfile 'save' with two signatures w/identical md5
    c = runtmp
    testdata1 = utils.get_test_data('img/2706795855.sig')
    testdata2 = utils.get_test_data('img/638277004.sig')

    c.run_sourmash('index', '-k', '21', 'zzz.sbt.zip', testdata1, testdata2)
    assert c.last_result.status == 0

    outfile = c.output('zzz.sbt.zip')
    assert os.path.exists(outfile)
    zout = zipfile.ZipFile(outfile, mode='r')
    # should have 3 files, 1 internal and two sigs. We check for 4 because the
    # directory also shows in namelist()
    assert len([f for f in zout.namelist() if f.startswith(".sbt.zzz/")]) == 5


def test_zipfile_does_not_exist(runtmp):
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash('sig', 'describe', 'no-exist.zip')

    # old behavior, pre PR #1777
    assert 'FileNotFoundError: SOURMASH-MANIFEST.csv' not in str(exc)
    assert not os.path.exists(runtmp.output('no-exist.zip'))

    # correct behavior
    assert "ERROR: Error while reading signatures from 'no-exist.zip'." in str(exc)


def test_zipfile_protein_command_search(runtmp):
    # test command-line search/gather of zipfile with protein sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_hp_command_search(runtmp):
    # test command-line search/gather of zipfile with hp sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/hp.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_dayhoff_command_search(runtmp):
    # test command-line search/gather of zipfile with dayhoff sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/dayhoff.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_protein_command_search_combined(runtmp):
    # test command-line search/gather of combined zipfile with protein sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_hp_command_search_combined(runtmp):
    # test command-line search/gather of combined zipfile with hp sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_dayhoff_command_search_combined(runtmp):
    # test command-line search/gather of combined zipfile with dayhoff sigs
    c = runtmp

    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_zipfile_dayhoff_command_search_protein(runtmp):
    # test command-line search/gather of protein sigs in zipfile
    c = runtmp

    # with dayhoff query
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.zip')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')

    print(c.last_result.out)
    print(c.last_result.err)

    assert 'no compatible signatures found in ' in c.last_result.err


def test_zipfile_API_signatures(use_manifest):
    # return all of the .sig and .sig.gz files in all.zip
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    siglist = list(zipidx.signatures())

    if use_manifest:
        assert len(siglist) == 8
        assert len(zipidx) == 8
    else:
        assert len(siglist) == 7
        assert len(zipidx) == 7


def test_zipfile_bool():
    # make sure that zipfile __bool__ doesn't traverse all the signatures
    # by relying on __len__!

    # create fake class that overrides everything useful except for bool -
    class FakeZipFileLinearIndex(ZipFileLinearIndex):
        def __init__(self):
            pass

        def signatures(self):
            yield 'a'
            raise Exception("don't touch me!")

        def __len__(self):
            raise Exception("don't call len!")

    # 'bool' should not touch __len__ or a second signature
    zf = FakeZipFileLinearIndex()
    assert bool(zf)

    # __len__ should raise an exception
    with pytest.raises(Exception) as exc:
        len(zf)
    assert "don't call len!" in str(exc.value)


def test_zipfile_API_signatures_traverse_yield_all(use_manifest):
    # include dna-sig.noext, but not build.sh (cannot be loaded as signature)
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, traverse_yield_all=True,
                                     use_manifest=use_manifest)
    siglist = list(zipidx.signatures())
    assert len(siglist) == 8
    assert len(zipidx) == 8

    # confirm that there are 12 files in there total, incl build.sh and dirs
    allfiles = zipidx.storage._filenames()
    print(allfiles)
    assert len(allfiles) == 13


def test_zipfile_API_signatures_traverse_yield_all_select(use_manifest):
    # include dna-sig.noext
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, traverse_yield_all=True,
                                     use_manifest=use_manifest)
    zipidx = zipidx.select(moltype='DNA')
    siglist = list(zipidx.signatures())
    assert len(siglist) == 2
    assert len(zipidx) == 2


def test_zipfile_API_signatures_traverse_yield_all_manifest():
    # check that manifest len is correct
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, traverse_yield_all=True,
                                     use_manifest=True)
    assert len(zipidx) == 8, len(zipidx)
    assert len(zipidx.manifest) == 8, len(zipidx.manifest)

    zipidx = zipidx.select(moltype='DNA')
    siglist = list(zipidx.signatures())
    assert len(siglist) == 2
    assert len(zipidx) == 2
    assert len(zipidx.manifest) == 2


def test_zipfile_API_signatures_select(use_manifest):
    # include dna-sig.noext
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    ziplist_pre = LinearIndex(zipidx.signatures())
    ziplist_pre = ziplist_pre.select(moltype='DNA')

    zipidx = zipidx.select(moltype='DNA')
    siglist = list(zipidx.signatures())

    if use_manifest:
        assert len(siglist) == 2
        assert len(zipidx) == 2
        assert len(ziplist_pre) == 2
    else:
        assert len(siglist) == 1
        assert len(zipidx) == 1
        assert len(ziplist_pre) == 1


def test_zipfile_API_signatures_select_abund_false(use_manifest):
    # check for abund=False (all signatures match b/c can convert)
    zipfile_db = utils.get_test_data('track_abund/track_abund.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    ziplist_pre = LinearIndex(zipidx.signatures())
    ziplist_pre = ziplist_pre.select(abund=False)

    zipidx = zipidx.select(abund=False)
    siglist = list(zipidx.signatures())

    assert len(siglist) == 2
    assert len(zipidx) == 2
    assert len(ziplist_pre) == 2


def test_zipfile_API_signatures_select_abund_true(use_manifest):
    # find all abund=True (all signatures match, b/c abund)
    zipfile_db = utils.get_test_data('track_abund/track_abund.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    ziplist_pre = LinearIndex(zipidx.signatures())
    ziplist_pre = ziplist_pre.select(abund=True)

    zipidx = zipidx.select(abund=True)
    siglist = list(zipidx.signatures())

    assert len(siglist) == 2
    assert len(zipidx) == 2
    assert len(ziplist_pre) == 2


def test_zipfile_API_signatures_select_abund_none(use_manifest):
    # find all abund=None (all signatures match, b/c no selection criteria)
    zipfile_db = utils.get_test_data('track_abund/track_abund.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    ziplist_pre = LinearIndex(zipidx.signatures())
    ziplist_pre = ziplist_pre.select(abund=None)

    zipidx = zipidx.select(abund=None)
    siglist = list(zipidx.signatures())

    assert len(siglist) == 2
    assert len(zipidx) == 2
    assert len(ziplist_pre) == 2


def test_zipfile_API_signatures_select_twice(use_manifest):
    # include dna-sig.noext
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)
    ziplist_pre = LinearIndex(zipidx.signatures())
    ziplist_pre = ziplist_pre.select(moltype='DNA')
    ziplist_pre = ziplist_pre.select(ksize=31)

    zipidx = zipidx.select(moltype='DNA')
    zipidx = zipidx.select(ksize=31)
    siglist = list(zipidx.signatures())

    if use_manifest:
        assert len(siglist) == 2
        assert len(zipidx) == 2
        assert len(ziplist_pre) == 2
    else:
        assert len(siglist) == 1
        assert len(zipidx) == 1
        assert len(ziplist_pre) == 1


def test_zipfile_API_save():
    # ZipFileLinearIndex.save is not implemented.
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)

    with pytest.raises(NotImplementedError):
        zipidx.save('xxx')


def test_zipfile_API_insert():
    # ZipFileLinearIndex.insert is not implemented.
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)

    with pytest.raises(NotImplementedError):
        # at some point probably want to change this to a real signature :)
        zipidx.insert(None)


def test_zipfile_API_location(use_manifest):
    # test ZipFileLinearIndex.location property
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, use_manifest=use_manifest)

    assert zipidx.location == zipfile_db


def test_zipfile_load_file_as_signatures(use_manifest):
    # make sure that ZipFileLinearIndex.signatures works, and is generator
    from types import GeneratorType

    zipfile_db = utils.get_test_data('prot/all.zip')
    sigs = sourmash_args.load_file_as_signatures(zipfile_db,
                                                 _use_manifest=use_manifest)

    # it's fine if this needs to change, but for now I want to make
    # sure that this is a generator.
    assert isinstance(sigs, GeneratorType)

    sigs = list(sigs)
    if use_manifest:
        assert len(sigs) == 8
    else:
        assert len(sigs) == 7


def test_zipfile_load_file_as_signatures_traverse_yield_all(use_manifest):
    # test with --force, which loads all files
    from types import GeneratorType

    zipfile_db = utils.get_test_data('prot/all.zip')
    sigs = sourmash_args.load_file_as_signatures(zipfile_db,
                                                 yield_all_files=True,
                                                 _use_manifest=use_manifest)

    # it's fine if this needs to change, but for now I want to make
    # sure that this is a generator.
    assert isinstance(sigs, GeneratorType)

    sigs = list(sigs)
    assert len(sigs) == 8


def test_zipfile_load_database_fail_if_not_zip(runtmp):
    # fail _load_database if not .zip
    c = runtmp

    zipfile_db = utils.get_test_data('prot/all.zip')
    badname = c.output('xyz.nada')
    shutil.copyfile(zipfile_db, badname)

    with pytest.raises(ValueError) as exc:
        sigs = sourmash_args.load_file_as_signatures(badname)

    assert 'Error while reading signatures from' in str(exc.value)


def test_multi_index_search():
    # test MultiIndex.search
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIndex with source location override
    lidx = MultiIndex.load([lidx1, lidx2, lidx3], ['A', None, 'C'],
                           None)
    lidx = lidx.select(ksize=31)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2
    assert sr[0][2] == 'A'      # source override

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[0][2] == sig47    # source was set to None, so no override
    assert sr[1][1] == ss63
    assert sr[1][2] == 'C'      # source override

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[0][2] == 'C'      # source override
    assert sr[1][1] == ss47
    assert sr[1][2] == sig47    # source was set to None, so no override

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[0][2] == 'C'      # source override


def test_multi_index_gather():
    # test MultiIndex.best_containment
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIndex with source location override
    lidx = MultiIndex.load([lidx1, lidx2, lidx3], ['A', None, 'C'],
                           None)
    lidx = lidx.select(ksize=31)

    match = lidx.best_containment(ss2)
    assert match
    assert match.score == 1.0
    assert match.location == 'A'

    match = lidx.best_containment(ss47)
    assert match
    assert match.score == 1.0
    assert match.signature == ss47
    assert match.location == sig47     # no source override


def test_multi_index_signatures():
    # test MultiIndex.signatures
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIndex with source location override
    lidx = MultiIndex.load([lidx1, lidx2, lidx3], ['A', None, 'C'],
                           None)
    lidx = lidx.select(ksize=31)

    siglist = list(lidx.signatures())
    assert len(siglist) == 3
    assert ss2 in siglist
    assert ss47 in siglist
    assert ss63 in siglist


def test_multi_index_create():
    # test MultiIndex constructor
    mi = MultiIndex(None, None, prepend_location=False)
    assert len(mi) == 0


def test_multi_index_create_prepend():
    # test MultiIndex constructor - location must be specified if
    # 'prepend_location is True
    with pytest.raises(ValueError):
        mi = MultiIndex(None, None, prepend_location=True)


def test_multi_index_load_from_directory():
    # test MultiIndex loading from a directory. The full paths to the
    # signature files should be available via 'signatures_with_location()'
    dirname = utils.get_test_data('prot/protein')
    mi = MultiIndex.load_from_directory(dirname, force=False)

    assert mi.location == dirname

    sigs = list(mi.signatures())
    assert len(sigs) == 2

    # check to make sure that full paths to expected sig files are returned
    locs = [ x[1] for x in mi.signatures_with_location() ]

    endings = ('GCA_001593925.1_ASM159392v1_protein.faa.gz.sig',
               'GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    for loc in locs:
        found = False
        for end in endings:
            if loc.endswith(end):
                found = True
        assert found, f"could not find full filename in locations for {end}"

    # also check internal locations and parent value --
    assert mi.parent.endswith('prot/protein')

    ilocs = [ x[1] for x in mi._signatures_with_internal() ]
    assert endings[0] in ilocs, ilocs
    assert endings[1] in ilocs, ilocs


def test_multi_index_load_from_directory_2():
    # only load .sig files, currently; not the databases under that directory.
    dirname = utils.get_test_data('prot')
    mi = MultiIndex.load_from_directory(dirname, force=False)

    sigs = list(mi.signatures())
    assert len(sigs) == 7


def test_multi_index_load_from_directory_3_simple_bad_file(runtmp):
    # check that force=False fails properly when confronted with non-JSON
    # files.
    c = runtmp

    with open(runtmp.output('badsig.sig'), 'wt') as fp:
        fp.write('bad content.')

    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_directory(runtmp.location, force=False)


def test_multi_index_load_from_directory_3(runtmp):
    # check that force=False fails properly when confronted with non-JSON
    # files that are legit sourmash files...
    c = runtmp

    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.sig")
            shutil.copyfile(fullname, copyto)
            count += 1

    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_directory(c.location, force=False)


def test_multi_index_load_from_directory_3_yield_all_true(runtmp):
    # check that force works ok on a directory w/force=True
    # Note here that only .sig/.sig.gz files are loaded.
    c = runtmp

    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.something")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_directory(c.location, force=True)

    sigs = list(mi.signatures())
    assert len(sigs) == 8


def test_multi_index_load_from_directory_3_yield_all_true_subdir(runtmp):
    # check that force works ok on subdirectories.
    # Note here that only .sig/.sig.gz files are loaded.
    c = runtmp
    dirname = utils.get_test_data('prot')

    target_dir = c.output("some_subdir")
    os.mkdir(target_dir)

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = os.path.join(target_dir, f"file{count}.something")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_directory(c.location, force=True)

    locations = set([ row['internal_location'] for row in mi.manifest.rows ])
    print(locations)

    sigs = list(mi.signatures())
    assert len(sigs) == 8


def test_multi_index_load_from_directory_3_sig_gz(runtmp):
    # check that we find .sig.gz files, too
    c = runtmp

    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            if not name.endswith('.sig'): # skip non .sig things
                continue
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.sig.gz")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_directory(c.location, force=False)

    assert mi.location == c.location

    sigs = list(mi.signatures())
    assert len(sigs) == 6


def test_multi_index_load_from_directory_3_check_traverse_fn(runtmp):
    # test the actual traverse function... eventually this test can be
    # removed, probably, as we consolidate functionality and test MultiIndex
    # better.
    c = runtmp

    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname]))
    assert len(files) == 7, files

    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 20, files # if this fails, check for extra files!


def test_multi_index_load_from_directory_no_exist():
    # raise ValueError on files that don't exist in load_from_directory
    dirname = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_directory(dirname, force=True)


def test_multi_index_load_from_file_path():
    # test that MultiIndex.load_from_path works fine
    sig2 = utils.get_test_data('2.fa.sig')

    mi = MultiIndex.load_from_path(sig2)
    assert len(mi) == 3
    assert mi.location == sig2


def test_multi_index_load_from_file_path_no_exist():
    # test that load_from_path fails on non-existent files
    filename = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_directory(filename, force=True)


def test_multi_index_load_from_pathlist_no_exist():
    # test that load_from_pathlist fails on non-existent files
    dirname = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_pathlist(dirname)


def test_multi_index_load_from_pathlist_1(runtmp):
    # test functionality of MultiIndex.load_from_pathlist with .sig files
    c = runtmp

    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname]))
    assert len(files) == 7, files

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)
    mi = MultiIndex.load_from_pathlist(file_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 7

    assert mi.location == file_list


def test_multi_index_load_from_pathlist_2(runtmp):
    # create a pathlist file with _all_ files under dir, and try to load it.
    # this will fail on one of several CSV or .sh files in there.

    # CTB note: if you create extra files under this directory,
    # it will fail :)
    c = runtmp
    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 20, files # check there aren't extra files in here!

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)

    with pytest.raises(ValueError) as exc:
        mi = MultiIndex.load_from_pathlist(file_list)

    print(str(exc))
    assert 'Error while reading signatures from' in str(exc)


def test_multi_index_load_from_pathlist_3_zipfile(runtmp):
    # can we load zipfiles in a pathlist? yes please.
    c = runtmp

    zipfile = utils.get_test_data('prot/all.zip')

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print(zipfile, file=fp)

    mi = MultiIndex.load_from_pathlist(file_list)
    assert len(mi) == 8

##
## test a slightly outre version of JaccardSearch - this is a test of the
## JaccardSearch 'collect' protocol, in particular...
##

class JaccardSearchBestOnly_ButIgnore(JaccardSearch):
    "A class that ignores certain results, but still does all the pruning."
    def __init__(self, ignore_list):
        super().__init__(SearchType.JACCARD, threshold=0.1)
        self.ignore_list = ignore_list

    # a collect function that _ignores_ things in the ignore_list
    def collect(self, score, match):
        print('in collect; current threshold:', self.threshold)
        for q in self.ignore_list:
            print('ZZZ', match, match.similarity(q))
            if match.similarity(q) == 1.0:
                print('yes, found.')
                return False

        # update threshold if not perfect match, which could help prune.
        self.threshold = score
        return True


def test_linear_index_gather_ignore():
    # do we properly ignore exact matches in 'search' for LinearIndex?
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    lidx = LinearIndex([ss2, ss47, ss63])

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(lidx.find(search_fn, ss47))
    results = [ sr.signature for sr in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_lca_index_gather_ignore():
    # do we properly ignore exact matches in gather on an LCA DB?
    from sourmash.lca import LCA_Database

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    db = LCA_Database(ksize=31, scaled=1000)
    db.insert(ss2)
    db.insert(ss47)
    db.insert(ss63)

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(db.find(search_fn, ss47))
    results = [ sr.signature for sr in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_sbt_index_gather_ignore():
    # do we properly ignore exact matches in gather on an SBT?
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    factory = GraphFactory(5, 100, 3)
    db = SBT(factory, d=2)

    db.insert(ss2)
    db.insert(ss47)
    db.insert(ss63)

    # ...now search with something that should ignore sig47, the exact match.
    print(f'\n** trying to ignore {ss47}')
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(db.find(search_fn, ss47))
    results = [ sr.signature for sr in results ]

    def is_found(ss, xx):
        for q in xx:
            print('is found?', ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_counter_gather_test_consume():
    # open-box testing of CounterGather.consume(...)
    # (see test_index_protocol.py for generic CounterGather tests.)
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss)
    counter.add(match_ss_1, location='loc a')
    counter.add(match_ss_2, location='loc b')
    counter.add(match_ss_3, location='loc c')

    ### ok, dig into actual counts...
    import pprint
    pprint.pprint(counter.counter)
    pprint.pprint(list(counter.signatures()))
    pprint.pprint(counter.locations)

    assert set(counter.signatures()) == set([match_ss_1, match_ss_2, match_ss_3])
    assert list(sorted(counter.locations.values())) == ['loc a', 'loc b', 'loc c']
    pprint.pprint(counter.counter.most_common())
    assert list(counter.counter.most_common()) == \
        [('26d4943627b33c446f37be1f5baf8d46', 10),
         ('f51cedec90ea666e0ebc11aa274eca61', 8),
         ('f331f8279113d77e42ab8efca8f9cc17', 4)]

    ## round 1

    cur_query = query_ss.minhash.to_mutable()
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_1
    assert len(intersect_mh) == 10
    assert cur_query == query_ss.minhash

    counter.consume(intersect_mh)
    assert set(counter.signatures()) == set([ match_ss_1, match_ss_2, match_ss_3 ])
    assert list(sorted(counter.locations.values())) == ['loc a', 'loc b', 'loc c']
    pprint.pprint(counter.counter.most_common())
    assert list(counter.counter.most_common()) == \
        [('f51cedec90ea666e0ebc11aa274eca61', 5),
         ('f331f8279113d77e42ab8efca8f9cc17', 4)]

    ### round 2

    cur_query.remove_many(intersect_mh.hashes)
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_2
    assert len(intersect_mh) == 5
    assert cur_query != query_ss.minhash

    counter.consume(intersect_mh)
    assert set(counter.signatures()) == set([ match_ss_1, match_ss_2, match_ss_3 ])
    assert list(sorted(counter.locations.values())) == ['loc a', 'loc b', 'loc c']

    pprint.pprint(counter.counter.most_common())
    assert list(counter.counter.most_common()) == \
        [('f331f8279113d77e42ab8efca8f9cc17', 2)]

    ## round 3

    cur_query.remove_many(intersect_mh.hashes)
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_3
    assert len(intersect_mh) == 2
    assert cur_query != query_ss.minhash

    counter.consume(intersect_mh)
    assert set(counter.signatures()) == set([ match_ss_1, match_ss_2, match_ss_3 ])
    assert list(sorted(counter.locations.values())) == ['loc a', 'loc b', 'loc c']
    pprint.pprint(counter.counter.most_common())
    assert list(counter.counter.most_common()) == []

    ## round 4 - nothing left!

    cur_query.remove_many(intersect_mh.hashes)
    results = counter.peek(cur_query)
    assert not results

    counter.consume(intersect_mh)
    assert set(counter.signatures()) == set([ match_ss_1, match_ss_2, match_ss_3 ])
    assert list(sorted(counter.locations.values())) == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.most_common()) == []


def test_counter_gather_identical_md5sum():
    # open-box testing of CounterGather.consume(...)
    # check what happens with identical matches w/different names
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # same as match_mh_1
    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(0, 10))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    # identical md5sum
    assert match_ss_1.md5sum() == match_ss_2.md5sum()

    # load up the counter
    counter = CounterGather(query_ss)
    counter.add(match_ss_1, location='loc a')
    counter.add(match_ss_2, location='loc b')

    assert len(counter.siglist) == 1
    stored_match = list(counter.siglist.values()).pop()
    assert stored_match.name == 'match2'
    # CTB note: this behavior may be changed freely, as the protocol
    # tests simply specify that _one_ of the identical matches is
    # returned. See test_counter_gather_multiple_identical_matches.


def test_lazy_index_1():
    # test some basic features of LazyLinearIndex
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)
    assert len(list(lazy2.signatures())) == 3

    results = lazy2.search(ss2, threshold=1.0)
    assert len(results) == 1
    assert results[0].signature == ss2


def test_lazy_index_2():
    # test laziness by adding a signature that raises an exception when
    # touched.

    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    lidx = LinearIndex()
    lidx.insert(FakeSignature())

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)

    sig_iter = lazy2.signatures()
    with pytest.raises(Exception) as e:
        list(sig_iter)

    assert str(e.value) == "don't touch me!"


def test_lazy_index_3():
    # make sure that you can't do multiple _incompatible_ selects.
    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    lidx = LinearIndex()
    lidx.insert(FakeSignature())

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)
    with pytest.raises(ValueError) as e:
        lazy3 = lazy2.select(ksize=21)

    assert str(e.value) == "cannot select on two different values for ksize"


def test_lazy_index_4_bool():
    # test some basic features of LazyLinearIndex
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    # test bool false/true
    lidx = LinearIndex()
    lazy = LazyLinearIndex(lidx)
    assert not lazy

    lidx.insert(ss2)
    assert lazy


def test_lazy_index_wraps_multi_index_location():
    # check that 'location' works fine when MultiIndex is wrapped by
    # LazyLinearIndex.
    sigdir = utils.get_test_data('prot/protein/')
    sigzip = utils.get_test_data('prot/protein.zip')
    siglca = utils.get_test_data('prot/protein.lca.json.gz')
    sigsbt = utils.get_test_data('prot/protein.sbt.zip')

    db_paths = (sigdir, sigzip, siglca, sigsbt)
    dbs = [ sourmash.load_file_as_index(db_path) for db_path in db_paths ]

    mi = MultiIndex.load(dbs, db_paths, None)
    lazy = LazyLinearIndex(mi)

    mi2 = mi.select(moltype='protein')
    lazy2 = lazy.select(moltype='protein')

    for (ss_tup, ss_lazy_tup) in zip(mi2.signatures_with_location(),
                                     lazy2.signatures_with_location()):
        assert ss_tup == ss_lazy_tup

def test_revindex_index_search():
    # confirm that RevIndex works
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_revindex_gather():
    # check that RevIndex.best_containment works.
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    match = lidx.best_containment(ss2)
    assert match
    assert match.score == 1.0
    assert match.signature == ss2

    match = lidx.best_containment(ss47)
    assert match
    assert match.score == 1.0
    assert match.signature == ss47


def test_revindex_gather_ignore():
    # check that RevIndex gather ignores things properly.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    lidx = RevIndex(template=ss2.minhash, signatures=[ss2, ss47, ss63])

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(lidx.find(search_fn, ss47))
    results = [ ss.signature for ss in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_standalone_manifest_signatures(runtmp):
    # build a StandaloneManifestIndex and test 'signatures' method.

    ## first, build a manifest in memory using MultiIndex
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig47)
    lidx2 = LinearIndex.load(sig63)

    mi = MultiIndex.load([lidx1, lidx2], [sig47, sig63], "")

    ## got a manifest! ok, now test out StandaloneManifestIndex
    mm = StandaloneManifestIndex(mi.manifest, None)

    siglist = [ ss for ss in mm.signatures() ]
    assert len(siglist) == 2
    assert ss47 in siglist
    assert ss63 in siglist


def test_standalone_manifest_signatures_prefix(runtmp):
    # try out 'prefix' for StandaloneManifestIndex

    ## first, build a manifest in memory using MultiIndex
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig47)
    lidx2 = LinearIndex.load(sig63)
    mi = MultiIndex.load([lidx1, lidx2], [sig47, sig63], "")

    # ok, now remove the abspath prefix from iloc
    for row in mi.manifest.rows:
        row['internal_location'] = os.path.basename(row['internal_location'])

    ## this should succeed!
    mm = StandaloneManifestIndex(mi.manifest, None,
                                 prefix=utils.get_test_data(''))

    assert len(list(mm.signatures())) == 2


def test_standalone_manifest_signatures_prefix_fail(runtmp):
    # give StandaloneManifest the wrong prefix

    ## first, build a manifest in memory using MultiIndex
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig47)
    lidx2 = LinearIndex.load(sig63)
    print('XXX', lidx1.location)

    mi = MultiIndex.load([lidx1, lidx2], [sig47, sig63], "")

    # remove prefix from manifest
    for row in mi.manifest.rows:
        row['internal_location'] = os.path.basename(row['internal_location'])

    ## got a manifest! ok, now test out StandaloneManifestIndex
    mm = StandaloneManifestIndex(mi.manifest, None, prefix='foo')

    # should fail
    with pytest.raises(ValueError) as exc:
        list(mm.signatures())

    assert "Error while reading signatures from 'foo/47.fa.sig'" in str(exc)


def test_standalone_manifest_load_from_dir(runtmp):
    # test loading a mf with relative directory paths from test-data
    mf = utils.get_test_data('scaled/mf.csv')
    idx = sourmash.load_file_as_index(mf)

    siglist = list(idx.signatures())
    assert len(siglist) == 15

    assert idx                  # should be 'True'
    assert len(idx) == 15

    with pytest.raises(NotImplementedError):
        idx.insert()

    with pytest.raises(NotImplementedError):
        idx.save('foo')

    assert idx.location == mf


def test_standalone_manifest_lazy_load(runtmp):
    # check that it's actually doing lazy loading
    orig_sig47 = utils.get_test_data('47.fa.sig')
    sig47 = runtmp.output('47.fa.sig')

    # build an external manifest
    shutil.copyfile(orig_sig47, sig47)

    # this is an abspath to sig47
    runtmp.sourmash('sig', 'manifest', sig47, '-o', 'mf.csv')

    # should work to get signatures:
    idx = StandaloneManifestIndex.load(runtmp.output('mf.csv'))

    siglist = list(idx.signatures())
    assert len(siglist) == 1

    # now remove!
    os.unlink(sig47)

    # can still access manifest...
    assert len(idx) == 1

    # ...but we should get an error when we call signatures.
    with pytest.raises(ValueError):
        list(idx.signatures())

    # but put it back, and all is forgiven. yay!
    shutil.copyfile(orig_sig47, sig47)
    x = list(idx.signatures())
    assert len(x) == 1


def test_standalone_manifest_lazy_load_2_prefix(runtmp):
    # check that it's actually doing lazy loading; supply explicit prefix
    orig_sig47 = utils.get_test_data('47.fa.sig')
    sig47 = runtmp.output('47.fa.sig')

    # build an external manifest
    # note, here use a relative path to 47.fa.sig; the manifest will contain
    # just '47.fa.sig' as the location
    shutil.copyfile(orig_sig47, sig47)
    runtmp.sourmash('sig', 'manifest', '47.fa.sig', '-o', 'mf.csv')

    # should work to get signatures:
    idx = StandaloneManifestIndex.load(runtmp.output('mf.csv'),
                                       prefix=runtmp.output(''))

    siglist = list(idx.signatures())
    assert len(siglist) == 1

    # now remove!
    os.unlink(sig47)

    # can still access manifest...
    assert len(idx) == 1

    # ...but we should get an error when we call signatures.
    with pytest.raises(ValueError):
        list(idx.signatures())

    # but put it back, and all is forgiven. yay!
    shutil.copyfile(orig_sig47, sig47)
    x = list(idx.signatures())
    assert len(x) == 1


def test_standalone_manifest_search(runtmp):
    # test a straight up 'search'
    query_sig = utils.get_test_data('scaled/genome-s12.fa.gz.sig')
    mf = utils.get_test_data('scaled/mf.csv')

    runtmp.sourmash('search', query_sig, mf)

    out = runtmp.last_result.out
    print(out)
    assert '100.0%       d84ef28f' in out


def test_standalone_manifest_prefetch_lazy(runtmp):
    # check that prefetch is actually doing lazy loading on manifest index.
    orig_sig47 = utils.get_test_data('47.fa.sig')
    sig47 = runtmp.output('47.fa.sig')
    orig_sig2 = utils.get_test_data('2.fa.sig')
    sig2 = runtmp.output('2.fa.sig')
    orig_sig63 = utils.get_test_data('63.fa.sig')
    sig63 = runtmp.output('63.fa.sig')

    shutil.copyfile(orig_sig47, sig47)
    runtmp.sourmash('sig', 'manifest', sig47, '-o', 'mf1.csv')
    shutil.copyfile(orig_sig2, sig2)
    runtmp.sourmash('sig', 'manifest', sig2, '-o', 'mf2.csv')
    shutil.copyfile(orig_sig63, sig63)
    runtmp.sourmash('sig', 'manifest', sig63, '-o', 'mf3.csv')

    # combine the manifests, manually for now...
    mf1 = CollectionManifest.load_from_filename(runtmp.output('mf1.csv'))
    assert len(mf1) == 1

    mf2 = CollectionManifest.load_from_filename(runtmp.output('mf2.csv'))
    assert len(mf2) == 3

    mf3 = CollectionManifest.load_from_filename(runtmp.output('mf3.csv'))
    assert len(mf3) == 1

    mf = mf1 + mf2 + mf3
    assert len(mf) == 5

    mf.write_to_filename(runtmp.output('mf.csv'))

    # ok! now, remove the last signature, 'sig63'.
    os.unlink(sig63)

    # ...but loading the manifest should still work.
    idx = StandaloneManifestIndex.load(runtmp.output('mf.csv'))

    # double check - third load will fail. this relies on load order :shrug:.
    sig_iter = iter(idx.signatures())
    ss = next(sig_iter)
    print(ss)
    assert '47.fa' in ss.filename

    for i in range(3):
        ss = next(sig_iter)
        print(i, ss)
        assert '2.fa' in ss.filename

    with pytest.raises(ValueError) as exc:
        ss = next(sig_iter)
    assert 'Error while reading signatures from' in str(exc)
    assert '63.fa.sig' in str(exc)

    # ok! now test prefetch... should get one match legit, to 47,
    # and then no matches to 2, and then error.

    ss47 = sourmash.load_one_signature(sig47)
    idx = idx.select(ksize=31)
    g = idx.prefetch(ss47, threshold_bp=0)

    # first value:
    sr = next(g)
    assert sr.signature == ss47

    # second value should raise error.
    with pytest.raises(ValueError) as exc:
        sr = next(g)

    assert 'Error while reading signatures from' in str(exc)
    assert '63.fa.sig' in str(exc)
