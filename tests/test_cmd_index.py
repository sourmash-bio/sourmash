"""
Tests for 'sourmash index'.
"""

import glob
import os

import pytest
import sourmash_tst_utils as utils

import sourmash
from sourmash import sourmash_args, SourmashSignature
from sourmash.sourmash_args import load_one_signature


def _index_filename(prefix, index_type):
    if index_type == "SBT":
        return prefix + ".sbt.zip"
    elif index_type == "rocksdb":
        return prefix + ".rocksdb"
    elif index_type == "zip":
        return prefix + ".sig.zip"

    raise Exception(f"unknown index type: {index_type}")


def test_index_signatures(runtmp, disk_index_type):
    # test 'signatures' method from Index base class
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    index_name = _index_filename("zzz", disk_index_type)
    runtmp.run_sourmash(
        "index", "-k", "31", index_name, sig47, sig63, "-F", disk_index_type
    )

    print(runtmp.last_result)
    print("loading from:", runtmp.output(index_name))

    db = sourmash.load_file_as_index(runtmp.output(index_name))

    xx = list(db.signatures())
    assert len(xx) == 2

    print(xx)

    ss47 = sourmash_args.load_query_signature(sig47, 31, "DNA")
    assert ss47 in xx
    ss63 = sourmash_args.load_query_signature(sig47, 31, "DNA")
    assert ss63 in xx


def test_search_metagenome(runtmp, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    index_name = _index_filename("zzz", disk_index_type)
    cmd = ["index", index_name, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output(index_name))

    runtmp.sourmash("search", query_sig, index_name, "-k", "21")

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert (
        " 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T"
        in runtmp.last_result.out
    )
    assert (
        "12 matches above threshold 0.080; showing first 3:" in runtmp.last_result.out
    )


# explanation: you cannot downsample a scaled index to match a scaled
# signature, so make sure that when you try such a search, it fails!
# (you *can* downsample a signature to match an index.)
def test_search_metagenome_index_downsample_fail(runtmp, disk_index_type):
    if disk_index_type in ("zip", "rocksdb"):
        raise pytest.skip("this actually works fine on these files")
    # test downsample on index => failure, with --fail-on-empty-databases
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    db_out = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", "-F", disk_index_type, db_out]
    cmd.extend(testdata_sigs)
    cmd.extend(["-k", "21"])

    print(" ".join(cmd))

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output(db_out))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("search", query_sig, db_out, "-k", "21", "--scaled", "100000")

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == -1
    assert f"ERROR: cannot use '{db_out}' for this query." in runtmp.last_result.err
    assert (
        "search scaled value 100000 is less than database scaled value of 10000"
        in runtmp.last_result.err
    )


def test_search_metagenome_downsample_containment(runtmp, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    db_out = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", db_out, *testdata_sigs, "-k", "21", "-F", disk_index_type]

    runtmp.sourmash(*cmd)

    assert os.path.exists(db_out)

    runtmp.sourmash(
        "search",
        query_sig,
        db_out,
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


def test_search_metagenome_downsample_index(runtmp, disk_index_type):
    # does same search as search_metagenome_downsample_containment but
    # rescales during indexing

    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    db = runtmp.output(_index_filename("gcf_all", disk_index_type))

    # downscale during indexing, rather than during search.
    runtmp.run_sourmash(
        "index",
        db,
        *testdata_sigs,
        "-k",
        "21",
        "--scaled",
        "100000",
        "-F",
        disk_index_type,
    )

    assert os.path.exists(db)

    runtmp.run_sourmash(
        "search",
        query_sig,
        db,
        "-k",
        "21",
        "--containment",
        "--scaled",
        "100000",
    )
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


def test_gather(runtmp, linear_gather, prefetch_gather, disk_index_type):
    testdata1 = utils.get_test_data("short.fa")
    testdata2 = utils.get_test_data("short2.fa")

    runtmp.sourmash("sketch", "dna", "-p", "scaled=10", testdata1, testdata2)

    runtmp.sourmash("sketch", "dna", "-p", "scaled=10", "-o", "query.fa.sig", testdata2)

    dbname = runtmp.output(_index_filename("zzz", disk_index_type))
    runtmp.sourmash(
        "index",
        "-k",
        "31",
        dbname,
        "short.fa.sig",
        "short2.fa.sig",
        "-F",
        disk_index_type,
    )

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        "query.fa.sig",
        dbname,
        "-o",
        "foo.csv",
        "--threshold-bp=1",
        linear_gather,
        prefetch_gather,
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "0.9 kbp      100.0%  100.0%" in runtmp.last_result.out


def test_gather_metagenome(runtmp, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]

    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash("gather", query_sig, dbname, "-k", "21", "--threshold-bp=0")

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


def test_gather_metagenome_num_results(runtmp, disk_index_type):
    # set a threshold on the number of results to be reported by gather
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]

    runtmp.run_sourmash(*cmd)

    assert os.path.exists(dbname)

    cmd = f"gather {query_sig} {dbname} -k 21 --num-results 10"
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


def test_gather_metagenome_threshold_bp(
    runtmp, linear_gather, prefetch_gather, disk_index_type
):
    # set a threshold on the gather output
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
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


def test_gather_metagenome_threshold_bp_low(
    runtmp, linear_gather, prefetch_gather, disk_index_type
):
    # set a threshold on the gather output => too low
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
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
    runtmp, linear_gather, prefetch_gather, disk_index_type
):
    # set a threshold on the gather output => no results
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
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


def test_gather_metagenome_abund(runtmp, disk_index_type):
    testdata_glob = utils.get_test_data("track_abund/*.fa.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("SRR606249.sig.gz")

    dbname = runtmp.output(_index_filename("against", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "31", "-F", disk_index_type]

    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash("gather", query_sig, dbname, "-k", "31", "--threshold-bp=0")

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert (
        "4.4 Mbp        0.6%  100.0%      10.4    NC_011663.1 "
        in runtmp.last_result.out
    )


def test_gather_metagenome_downsample(
    runtmp, prefetch_gather, linear_gather, disk_index_type
):
    # downsample w/scaled of 100,000
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
        "-k",
        "21",
        "--scaled",
        "100000",
        prefetch_gather,
        linear_gather,
        "--threshold-bp",
        "50000",
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 11 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out
    assert all(
        (
            "5.2 Mbp       32.9%  100.0%" in runtmp.last_result.out,
            "NC_003198.1" in runtmp.last_result.out,
        )
    )
    assert all(
        (
            "4.1 Mbp        0.6%    2.4%" in runtmp.last_result.out,
            "4.1 Mbp        4.4%   17.1%" in runtmp.last_result.out,
        )
    )


def test_gather_save_matches(runtmp, linear_gather, prefetch_gather, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
        "-k",
        "21",
        "--save-matches",
        "save.sigs",
        linear_gather,
        prefetch_gather,
        "--threshold-bp",
        "0",
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 12 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out
    assert os.path.exists(runtmp.output("save.sigs"))


def test_gather_save_matches_and_save_prefetch(runtmp, linear_gather, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    runtmp.sourmash(
        "gather",
        query_sig,
        dbname,
        "-k",
        "21",
        "--save-matches",
        "save.sigs",
        "--save-prefetch",
        "save2.sigs",
        linear_gather,
        "--threshold-bp",
        "0",
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 12 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out

    matches_save = runtmp.output("save.sigs")
    prefetch_save = runtmp.output("save2.sigs")
    assert os.path.exists(matches_save)
    assert os.path.exists(prefetch_save)

    matches = list(sourmash.load_file_as_signatures(matches_save))
    prefetch = list(sourmash.load_file_as_signatures(prefetch_save))

    assert set(matches) == set(prefetch)


def test_gather_metagenome_picklist(runtmp, disk_index_type):
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data("gather/combined.sig")

    dbname = runtmp.output(_index_filename("gcf_all", disk_index_type))
    cmd = ["index", dbname, *testdata_sigs, "-k", "21", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    assert os.path.exists(dbname)

    pl = utils.get_test_data("gather/salmonella-picklist.csv")
    runtmp.sourmash(
        "gather", query_sig, dbname, "-k", "21", "--picklist", f"{pl}:name:ident"
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "found 7 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 58.3% of the query" in runtmp.last_result.out
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
    assert (
        "for given picklist, found 8 matches to 8 distinct values"
        in runtmp.last_result.err
    )


def test_index_best_containment_threshold_1(runtmp, disk_index_type):
    # test best_containment() method, in some detail
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = sourmash_args.load_query_signature(sig2, 31, "DNA")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    dbname = runtmp.output(_index_filename("test", disk_index_type))
    cmd = ["index", dbname, sig2, sig47, sig63, "-k", "31", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    db = sourmash.load_file_as_index(dbname)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(ss2.minhash.hashes.keys()))
    new_mh = ss2.minhash.copy_and_clear()

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
    assert match_sig == ss2

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
    assert match_sig == ss2

    # check with a too-high threshold -> should be no results.
    print("len mh", len(new_mh))
    with pytest.raises(ValueError):
        db.best_containment(SourmashSignature(new_mh), threshold_bp=5000)


def test_index_best_containment_threshold_5(runtmp, disk_index_type):
    sig2 = utils.get_test_data("2.fa.sig")
    ss2 = sourmash_args.load_query_signature(sig2, 31, "DNA")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    dbname = runtmp.output(_index_filename("test", disk_index_type))
    cmd = ["index", dbname, sig2, sig47, sig63, "-k", "31", "-F", disk_index_type]
    runtmp.sourmash(*cmd)

    db = sourmash.load_file_as_index(dbname)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(ss2.minhash.hashes.keys()))
    new_mh = ss2.minhash.copy_and_clear()

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
    assert match_sig == ss2

    # now, check with a threshold_bp that should be meet-able.
    db.best_containment(SourmashSignature(new_mh), threshold_bp=5000)
    assert result
    containment, match_sig, name = result
    assert containment == 1.0
    assert match_sig == ss2


def test_gather_single_return(runtmp, disk_index_type):
    # test gather() number of returns
    sig2file = utils.get_test_data("2.fa.sig")
    sig47file = utils.get_test_data("47.fa.sig")
    sig63file = utils.get_test_data("63.fa.sig")

    dbname = runtmp.output(_index_filename("db", disk_index_type))

    runtmp.sourmash(
        "index",
        dbname,
        sig2file,
        sig47file,
        sig63file,
        "-k",
        "31",
        "-F",
        disk_index_type,
    )

    db = sourmash.load_file_as_index(dbname)

    # now, run best_containment. Right match?
    sig63 = sourmash_args.load_query_signature(sig63file, 31, "DNA")
    result = db.best_containment(sig63)
    print(result)
    assert result
    assert result.score == 1.0


def test_sbt_jaccard_ordering(runtmp, disk_index_type):
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

    # now - make signatures, build index, try out.
    ss_a = sourmash.SourmashSignature(a, name="A")
    ss_b = sourmash.SourmashSignature(b, name="B")
    ss_c = sourmash.SourmashSignature(c, name="C")

    sigsfile = runtmp.output("insigs.sig.zip")
    with sourmash_args.SaveSignaturesToLocation(sigsfile) as save_sigs:
        save_sigs.add(ss_a)
        save_sigs.add(ss_b)
        save_sigs.add(ss_c)

    index_name = _index_filename("db", disk_index_type)
    runtmp.sourmash(
        "index", "-F", disk_index_type, index_name, sigsfile, "--scaled", "2"
    )

    db = sourmash.load_file_as_index(runtmp.output(index_name))

    sr = db.search(ss_a, threshold=0.15)
    print(sr)
    assert len(sr) == 2
    assert sr[0].signature == ss_a
    assert sr[0].score == 1.0
    assert sr[1].signature == ss_c
    assert sr[1].score == 0.2


def test_index_protein(runtmp, disk_index_type):
    # test command-line creation of databases with protein sigs
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = runtmp.output(_index_filename("protein", disk_index_type))

    runtmp.run_sourmash(
        "index",
        db_out,
        sigfile1,
        sigfile2,
        "--scaled",
        "100",
        "-k",
        "19",
        "--protein",
        "-F",
        disk_index_type,
    )
    assert os.path.exists(db_out), db_out

    db2 = sourmash.load_file_as_index(db_out)

    sig1 = load_one_signature(sigfile1)
    sig2 = load_one_signature(sigfile2)

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
    assert result.location == db2.location
    assert result.location == db_out


def test_index_protein_search_no_threshold(runtmp, disk_index_type):
    # test the '.search' method on indexes w/no threshold
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = runtmp.output(_index_filename("protein", disk_index_type))

    runtmp.run_sourmash(
        "index",
        db_out,
        sigfile1,
        sigfile2,
        "--scaled",
        "100",
        "-k",
        "19",
        "--protein",
        "-F",
        disk_index_type,
    )

    db2 = sourmash.load_file_as_index(db_out)

    sig1 = load_one_signature(sigfile1)

    # and search, gather
    with pytest.raises(TypeError) as exc:
        db2.search(sig1)
    assert "'search' requires 'threshold'" in str(exc)


def test_index_protein_command_search(runtmp, disk_index_type):
    # test command-line search/gather of on-disk databases with protein sigs
    sigfile1 = utils.get_test_data(
        "prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )

    if disk_index_type == "zip":
        db_out = utils.get_test_data("prot/protein.zip")
    else:
        db_out = utils.get_test_data(_index_filename("prot/protein", disk_index_type))

    runtmp.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in runtmp.last_result.out

    runtmp.run_sourmash("gather", sigfile1, db_out)
    assert "found 1 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out


def test_index_hp_command_index(runtmp, disk_index_type):
    # test command-line creation of on-disk databases with hp sigs
    sigfile1 = utils.get_test_data(
        "prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = runtmp.output(_index_filename("hp", disk_index_type))

    runtmp.run_sourmash(
        "index",
        db_out,
        sigfile1,
        sigfile2,
        "--scaled",
        "100",
        "-k",
        "19",
        "--hp",
        "-F",
        disk_index_type,
    )

    db2 = sourmash.load_file_as_index(db_out)

    sig1 = load_one_signature(sigfile1)
    sig2 = load_one_signature(sigfile2)

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
    assert result.location == db2.location
    assert result.location == db_out


def test_index_hp_command_search(runtmp, disk_index_type):
    # test command-line search/gather of on-disk databases with hp sigs
    sigfile1 = utils.get_test_data(
        "prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    if disk_index_type == "zip":
        db_out = utils.get_test_data("prot/hp.zip")
    else:
        db_out = utils.get_test_data(_index_filename("prot/hp", disk_index_type))

    runtmp.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in runtmp.last_result.out

    runtmp.run_sourmash("gather", sigfile1, db_out, "--threshold", "0.0")
    assert "found 1 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out


def test_index_dayhoff_command_index(runtmp, disk_index_type):
    # test command-line creation of on-disk databases with dayhoff sigs
    sigfile1 = utils.get_test_data(
        "prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    sigfile2 = utils.get_test_data(
        "prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"
    )

    db_out = runtmp.output(_index_filename("dayhoff", disk_index_type))

    runtmp.run_sourmash(
        "index",
        db_out,
        sigfile1,
        sigfile2,
        "--scaled",
        "100",
        "-k",
        "19",
        "--dayhoff",
        "-F",
        disk_index_type,
    )

    db2 = sourmash.load_file_as_index(db_out)

    sig1 = load_one_signature(sigfile1)
    sig2 = load_one_signature(sigfile2)

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
    assert result.location == db2.location
    assert result.location == db_out


def test_index_dayhoff_command_search(runtmp, disk_index_type):
    # test command-line search/gather of on-disk databases with dayhoff sigs
    sigfile1 = utils.get_test_data(
        "prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig"
    )
    if disk_index_type == "zip":
        db_out = utils.get_test_data("prot/dayhoff.zip")
    else:
        db_out = utils.get_test_data(_index_filename("prot/dayhoff", disk_index_type))

    runtmp.run_sourmash("search", sigfile1, db_out, "--threshold", "0.0")
    assert "2 matches" in runtmp.last_result.out

    runtmp.run_sourmash("gather", sigfile1, db_out, "--threshold", "0.0")
    assert "found 1 matches total" in runtmp.last_result.out
    assert "the recovered matches hit 100.0% of the query" in runtmp.last_result.out


def test_index_skipm1n3_command_search(runtmp, disk_index_type):
    # test command-line search/gather of on-disk databases with skipm1n3 sigs
    sig2file = utils.get_test_data("skipmers/2.skip.sig.zip")
    sig47file = utils.get_test_data("skipmers/47.skip.sig.zip")
    sig63file = utils.get_test_data("skipmers/63.skip.sig.zip")

    db_out = runtmp.output(_index_filename("skipm1n3", disk_index_type))

    runtmp.sourmash(
        "index",
        "-F",
        disk_index_type,
        db_out,
        sig2file,
        sig47file,
        sig63file,
        "--skipm1n3",
        "-k",
        "31",
    )
    assert os.path.exists(runtmp.output(db_out))

    db = sourmash.load_file_as_index(db_out)
    sig47 = sourmash_args.load_query_signature(sig47file, 31, "skipm1n3")
    sig63 = sourmash_args.load_query_signature(sig63file, 31, "skipm1n3")

    # check reconstruction --
    mh_list = [x.minhash for x in db.signatures()]
    assert len(mh_list) == 3
    assert sig47.minhash in mh_list
    assert sig63.minhash in mh_list

    # and search, gather
    results = db.search(
        sig47,
        threshold=0.0,
        ignore_abundance=True,
        do_containment=False,
        best_only=False,
    )
    assert len(results) == 3

    result = db.best_containment(sig63)
    assert result.score == 1.0
    assert result.location == db.location
    assert result.location == db_out


def test_index_skipm2n3_command_search(runtmp, disk_index_type):
    # test command-line search/gather of on-disk databases with skipm2n3 sigs
    sig2file = utils.get_test_data("skipmers/2.skip.sig.zip")
    sig47file = utils.get_test_data("skipmers/47.skip.sig.zip")
    sig63file = utils.get_test_data("skipmers/63.skip.sig.zip")

    db_out = runtmp.output(_index_filename("skipm2n3", disk_index_type))

    runtmp.sourmash(
        "index",
        "-F",
        disk_index_type,
        db_out,
        sig2file,
        sig47file,
        sig63file,
        "--skipm2n3",
        "-k",
        "31",
    )
    assert os.path.exists(runtmp.output(db_out))

    db = sourmash.load_file_as_index(db_out)
    sig47 = sourmash_args.load_query_signature(sig47file, 31, "skipm2n3")
    sig63 = sourmash_args.load_query_signature(sig63file, 31, "skipm2n3")

    # check reconstruction --
    mh_list = [x.minhash for x in db.signatures()]
    assert len(mh_list) == 3
    assert sig47.minhash in mh_list
    assert sig63.minhash in mh_list

    # and search, gather
    results = db.search(
        sig47,
        threshold=0.0,
        ignore_abundance=True,
        do_containment=False,
        best_only=False,
    )
    assert len(results) == 3

    result = db.best_containment(sig63)
    assert result.score == 1.0
    assert result.location == db.location
    assert result.location == db_out
