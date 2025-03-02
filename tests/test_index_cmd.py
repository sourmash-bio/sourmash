import glob
import os

import pytest
import sourmash_tst_utils as utils

import sourmash
from sourmash import sourmash_args


def test_index_signatures(runtmp):
    # test 'signatures' method from Index base class
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    runtmp.run_sourmash("index", "-k", "31", "zzz.rocksdb", sig47, sig63)

    db = sourmash.load_file_as_index(runtmp.output("zzz.rocksdb"))

    xx = list(db.signatures())
    assert len(xx) == 2

    print(xx)

    ss47 = sourmash_args.load_query_signature(sig47, 31, "DNA")
    assert ss47 in xx
    ss63 = sourmash_args.load_query_signature(sig47, 31, "DNA")
    assert ss63 in xx


def test_search_metagenome(runtmp):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash("search", query_sig, "gcf_all.rocksdb", "-k", "21")

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert (
        " 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T"
        in runtmp.last_result.out
    )
    assert (
        "12 matches above threshold 0.080; showing first 3:" in runtmp.last_result.out
    )


# explanation: you cannot downsample a scaled SBT to match a scaled
# signature, so make sure that when you try such a search, it fails!
# (you *can* downsample a signature to match an SBT.)
def test_search_metagenome_index_downsample_fail(runtmp):
    raise pytest.xfail("mismatch scaled")
    # test downsample on SBT => failure, with --fail-on-empty-databases
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "search", query_sig, "gcf_all.rocksdb", "-k", "21", "--scaled", "100000"
        )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == -1
    assert "ERROR: cannot use 'gcf_all.rocksdb' for this query." in runtmp.last_result.err
    assert (
        "search scaled value 100000 is less than database scaled value of 10000"
        in runtmp.last_result.err
    )


def test_search_metagenome_downsample_containment(runtmp):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash(
        "search",
        query_sig,
        "gcf_all.rocksdb",
        "-k",
        "21",
        "--scaled",
        "100000",
        "--containment",
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert (
        " 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T"
        in runtmp.last_result.out
    )
    assert (
        "12 matches above threshold 0.080; showing first 3:" in runtmp.last_result.out
    )


def test_search_metagenome_downsample_index(runtmp):
    # does same search as search_metagenome_downsample_containment but
    # rescales during indexing

    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    # downscale during indexing, rather than during search.
    runtmp.run_sourmash("index", "gcf_all.rocksdb", *testdata_sigs, "-k", "21", "--scaled", "100000")

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.run_sourmash("search", query_sig, "gcf_all.rocksdb", "-k", "21", "--containment", "--scaled", "100000")
    print(runtmp)

    assert (
        " 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T"
        in str(runtmp)
    )
    assert (
        " 29.7%       NC_003197.2 Salmonella enterica subsp. enterica serovar T"
        in str(runtmp)
    )
    assert "12 matches above threshold 0.080; showing first 3:" in str(runtmp)


def test_gather(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data("short.fa")
    testdata2 = utils.get_test_data("short2.fa")

    runtmp.sourmash("sketch", "dna", "-p", "scaled=10", testdata1, testdata2)

    runtmp.sourmash("sketch", "dna", "-p", "scaled=10", "-o", "query.fa.sig", testdata2)

    runtmp.sourmash("index", "-k", "31", "zzz.rocksdb", "short.fa.sig", "short2.fa.sig")

    assert os.path.exists(runtmp.output("zzz.rocksdb"))

    runtmp.sourmash(
        "gather",
        "query.fa.sig",
        "zzz.rocksdb",
        "-o",
        "foo.csv",
        "--threshold-bp=1",
        linear_gather,
        prefetch_gather,
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "0.9 kbp      100.0%  100.0%" in runtmp.last_result.out


def test_gather_metagenome(runtmp):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash("gather", query_sig, "gcf_all.rocksdb", "-k", "21", "--threshold-bp=0")

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 12 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out
    assert all(
        (
            "4.9 Mbp       33.2%  100.0%" in runtmp.last_result.out,
            "NC_003198.1 Salmonella enterica subsp" in runtmp.last_result.out,
        )
    )
    assert all(
        (
            "4.7 Mbp        0.5%    1.5%" in runtmp.last_result.out,
            "NC_011294.1 Salmonella enterica subs" in runtmp.last_result.out,
        )
    )


def test_gather_metagenome_num_results(runtmp):
    # set a threshold on the number of results to be reported by gather
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.run_sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    cmd = f"gather {query_sig} gcf_all.rocksdb -k 21 --num-results 10"
    cmd = cmd.split(" ")
    runtmp.run_sourmash(*cmd)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    out = runtmp.last_result.out

    assert "found 10 matches total" in out
    assert "(truncated gather because --num-results=10)" in out
    assert "the recovered matches hit 99.4% of the query" in out
    assert all(
        (
            "4.9 Mbp       33.2%  100.0%" in out,
            "NC_003198.1 Salmonella enterica subsp" in out,
        )
    )
    assert "4.3 Mbp        2.1%    7.3%    NC_006511.1 Salmonella enterica subsp" in out


def test_gather_metagenome_threshold_bp(runtmp, linear_gather, prefetch_gather):
    # set a threshold on the gather output
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash(
        "gather",
        query_sig,
        "gcf_all.rocksdb",
        "-k",
        "21",
        "--threshold-bp",
        "2e6",
        linear_gather,
        prefetch_gather,
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 1 matches total" in runtmp.last_result.out
    assert "found less than 2.0 Mbp in common. => exiting" in runtmp.last_result.err
    assert "the recovered matches hit 33.2% of the query" in runtmp.last_result.out
    assert all(
        (
            "4.9 Mbp       33.2%  100.0%" in runtmp.last_result.out,
            "NC_003198.1 Salmonella enterica subsp" in runtmp.last_result.out,
        )
    )


def test_gather_metagenome_threshold_bp_low(runtmp, linear_gather, prefetch_gather):
    # set a threshold on the gather output => too low
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash(
        "gather",
        query_sig,
        "gcf_all.rocksdb",
        "-k",
        "21",
        "--threshold-bp",
        "1",
        linear_gather,
        prefetch_gather,
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 12 matches total" in runtmp.last_result.out
    assert "found less than 1 bp in common. => exiting" in runtmp.last_result.err
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out


def test_gather_metagenome_threshold_bp_too_high(
    runtmp, linear_gather, prefetch_gather
):
    # set a threshold on the gather output => no results
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    cmd = ["index", "gcf_all.rocksdb"]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output("gcf_all.rocksdb"))

    runtmp.sourmash(
        "gather",
        query_sig,
        "gcf_all.rocksdb",
        "-k",
        "21",
        "--threshold-bp",
        "5e6",
        linear_gather,
        prefetch_gather,
    )

    out = runtmp.last_result.out
    err = runtmp.last_result.err
    print(out)
    print(err)

    assert "No matches found for --threshold-bp at 5.0 Mbp." in err


"""
def test_sbt_gather_threshold_1():
    # test gather() method, in some detail
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig2 = load_one_signature(utils.get_test_data("2.fa.sig"), ksize=31)
    sig47 = load_one_signature(utils.get_test_data("47.fa.sig"), ksize=31)
    sig63 = load_one_signature(utils.get_test_data("63.fa.sig"), ksize=31)

    tree.insert(sig47)
    tree.insert(sig63)
    tree.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        tree.best_containment(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    result = tree.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        tree.best_containment(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    result = tree.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a too-high threshold -> should be no results.
    print("len mh", len(new_mh))
    with pytest.raises(ValueError):
        tree.best_containment(SourmashSignature(new_mh), threshold_bp=5000)


def test_sbt_gather_threshold_5():
    # test gather() method above threshold
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig2 = load_one_signature(utils.get_test_data("2.fa.sig"), ksize=31)
    sig47 = load_one_signature(utils.get_test_data("47.fa.sig"), ksize=31)
    sig63 = load_one_signature(utils.get_test_data("63.fa.sig"), ksize=31)

    tree.insert(sig47)
    tree.insert(sig63)
    tree.insert(sig2)

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
    result = tree.best_containment(SourmashSignature(new_mh))
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # now, check with a threshold_bp that should be meet-able.
    tree.best_containment(SourmashSignature(new_mh), threshold_bp=5000)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None


@utils.in_tempdir
def test_gather_single_return(c):
    # test gather() number of returns
    sig2file = utils.get_test_data("2.fa.sig")
    sig47file = utils.get_test_data("47.fa.sig")
    sig63file = utils.get_test_data("63.fa.sig")

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct SBT Database
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    tree.insert(sig2)
    tree.insert(sig47)
    tree.insert(sig63)

    # now, run gather. how many results do we get, and are they in the
    # right order?
    result = tree.best_containment(sig63)
    print(result)
    assert result
    assert result.score == 1.0


def test_sbt_jaccard_ordering(runtmp):
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

    print("a intersect b:", _intersect(a, b))
    print("a intersect c:", _intersect(a, c))
    print("a jaccard b:", a.jaccard(b))
    print("a jaccard c:", a.jaccard(c))
    assert _intersect(a, b) > _intersect(a, c)
    assert a.jaccard(b) < a.jaccard(c)

    # thresholds to use:
    assert a.jaccard(b) < 0.15
    assert a.jaccard(c) > 0.15

    # now - make signatures, try out :)
    ss_a = sourmash.SourmashSignature(a, name="A")
    ss_b = sourmash.SourmashSignature(b, name="B")
    ss_c = sourmash.SourmashSignature(c, name="C")

    factory = GraphFactory(31, 1e5, 4)
    db = SBT(factory, d=2)
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


def test_sbt_protein_command_index(runtmp):
    c = runtmp

    # test command-line creation of SBT database with protein sigs
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = c.output("protein.sbt.zip")

    c.run_sourmash(
        "index", db_out, sigfile1, sigfile2, "--scaled", "100", "-k", "19", "--protein"
    )

    # check to make sure .sbt.protein directory doesn't get created
    assert not os.path.exists(c.output(".sbt.protein"))

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [x.minhash for x in db2.signatures()]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(
        sig1,
        threshold=0.0,
        ignore_abundance=True,
        do_containment=False,
        best_only=False,
    )
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0
    assert result.location == db2._location
    assert result.location == db_out


@utils.in_tempdir
def test_sbt_protein_search_no_threshold(c):
    # test the '.search' method on SBTs w/no threshold
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = c.output("protein.sbt.zip")

    c.run_sourmash(
        "index", db_out, sigfile1, sigfile2, "--scaled", "100", "-k", "19", "--protein"
    )

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)

    # and search, gather
    with pytest.raises(TypeError) as exc:
        db2.search(sig1)
    assert "'search' requires 'threshold'" in str(exc)


@utils.in_thisdir
def test_sbt_protein_command_search(c):
    # test command-line search/gather of SBT database with protein sigs
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    db_out = utils.get_test_data("prot/protein.sbt.zip")

    c.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in c.last_result.out

    c.run_sourmash("gather", sigfile1, db_out)
    assert "found 1 matches total" in c.last_result.out
    assert "the recovered matches hit 100.0% of the query" in c.last_result.out


@utils.in_tempdir
def test_sbt_hp_command_index(c):
    # test command-line creation of SBT database with hp sigs
    sigfile1 = utils.get_test_data(
        "prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = c.output("hp.sbt.zip")

    c.run_sourmash(
        "index", db_out, sigfile1, sigfile2, "--scaled", "100", "-k", "19", "--hp"
    )

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [x.minhash for x in db2.signatures()]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(
        sig1,
        threshold=0.0,
        ignore_abundance=True,
        do_containment=False,
        best_only=False,
    )
    assert results

    result = db2.best_containment(sig2)
    assert result.score == 1.0
    assert result.location == db2._location
    assert result.location == db_out


@utils.in_thisdir
def test_sbt_hp_command_search(c):
    # test command-line search/gather of SBT database with hp sigs
    sigfile1 = utils.get_test_data(
        "prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    db_out = utils.get_test_data("prot/hp.sbt.zip")

    c.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in c.last_result.out

    c.run_sourmash("gather", sigfile1, db_out, "--threshold", "0.0")
    assert "found 1 matches total" in c.last_result.out
    assert "the recovered matches hit 100.0% of the query" in c.last_result.out


@utils.in_tempdir
def test_sbt_dayhoff_command_index(c):
    # test command-line creation of SBT database with dayhoff sigs
    sigfile1 = utils.get_test_data(
        "prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = c.output("dayhoff.sbt.zip")

    c.run_sourmash(
        "index", db_out, sigfile1, sigfile2, "--scaled", "100", "-k", "19", "--dayhoff"
    )

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [x.minhash for x in db2.signatures()]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(
        sig1,
        threshold=0.0,
        ignore_abundance=True,
        do_containment=False,
        best_only=False,
    )
    assert len(results) == 2

    result = db2.best_containment(sig2)
    assert result.score == 1.0
    assert result.location == db2._location
    assert result.location == db_out


@utils.in_thisdir
def test_sbt_dayhoff_command_search(c):
    # test command-line search/gather of SBT database with dayhoff sigs
    sigfile1 = utils.get_test_data(
        "prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    db_out = utils.get_test_data("prot/dayhoff.sbt.zip")

    c.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in c.last_result.out

    c.run_sourmash("gather", sigfile1, db_out, "--threshold", "0.0")
    assert "found 1 matches total" in c.last_result.out
    assert "the recovered matches hit 100.0% of the query" in c.last_result.out


"""
