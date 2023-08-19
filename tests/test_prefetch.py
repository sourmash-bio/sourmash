"""
Tests for `sourmash prefetch` command-line and API functionality.
"""
import os
import csv
import gzip
import pytest
import glob
import random
from sourmash.search import PrefetchResult

import sourmash_tst_utils as utils
import sourmash
from sourmash_tst_utils import SourmashCommandFailed
from sourmash import signature, sourmash_args


def approx_eq(val1, val2):
    if val1 == "":
        return val1 == val2
    return round(float(val1), 3) == round(float(val2), 3)


def test_prefetch_basic(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "query sketch has scaled=1000; will be dynamically downsampled as needed" in c.last_result.err

    err = c.last_result.err
    assert "loaded 5 total signatures from 3 locations." in err
    assert "after selecting signatures compatible with search, 3 remain." in err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_select_query_ksize(runtmp, linear_gather):
    # test prefetch where query and subject db both have multiple ksizes
    c = runtmp

    ss = utils.get_test_data('GCF_000005845.2_ASM584v2_genomic.fna.gz.sig')

    c.run_sourmash('prefetch', ss, ss, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'of 4476 distinct query hashes, 4476 were found in matches above threshold.' in c.last_result.err


def test_prefetch_subject_scaled_is_larger(runtmp, linear_gather):
    # test prefetch where subject scaled is larger
    c = runtmp

    # make a query sketch with scaled=1000
    fa = utils.get_test_data('genome-s10.fa.gz')
    c.run_sourmash('sketch', 'dna', fa, '-o', 'query.sig')
    assert os.path.exists(runtmp.output('query.sig'))

    # this has a scaled of 10000, from same genome:
    against1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    against2 = utils.get_test_data('scaled/all.sbt.zip')
    against3 = utils.get_test_data('scaled/all.lca.json')

    # run against large scaled, then small (self)
    c.run_sourmash('prefetch', 'query.sig', against1, against2, against3,
                   'query.sig', linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'total of 8 matching signatures.' in c.last_result.err
    assert 'of 48 distinct query hashes, 48 were found in matches above threshold.' in c.last_result.err
    assert 'final scaled value (max across query and all matches) is 10000' in c.last_result.err


def test_prefetch_subject_scaled_is_larger_outsigs(runtmp, linear_gather):
    # test prefetch where subject scaled is larger -- output sigs
    c = runtmp

    # make a query sketch with scaled=1000
    fa = utils.get_test_data('genome-s10.fa.gz')
    c.run_sourmash('sketch', 'dna', fa, '-o', 'query.sig')
    assert os.path.exists(runtmp.output('query.sig'))

    # this has a scaled of 10000, from same genome:
    against1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    against2 = utils.get_test_data('scaled/all.sbt.zip')
    against3 = utils.get_test_data('scaled/all.lca.json')

    # run against large scaled, then small (self)
    c.run_sourmash('prefetch', 'query.sig', against1, against2, against3,
                   'query.sig', linear_gather, '--save-matches', 'matches.sig')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'total of 8 matching signatures.' in c.last_result.err
    assert 'of 48 distinct query hashes, 48 were found in matches above threshold.' in c.last_result.err
    assert 'final scaled value (max across query and all matches) is 10000' in c.last_result.err

    # make sure non-downsampled sketches were saved.
    matches = sourmash.load_file_as_signatures(runtmp.output('matches.sig'))
    scaled_vals = set([ match.minhash.scaled for match in matches ])
    assert 1000 in scaled_vals
    assert 10000 in scaled_vals
    assert len(scaled_vals) == 2


def test_prefetch_query_abund(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch w/abund query
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "query sketch has scaled=1000; will be dynamically downsampled as needed" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_subj_abund(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch w/abund signature.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "query sketch has scaled=1000; will be dynamically downsampled as needed" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_csv_out(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with CSV output
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    expected_intersect_bp = [2529000, 5177000]
    with open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_intersect_bp):
            print(row)
            assert int(row['intersect_bp']) == expected


def test_prefetch_csv_gz_out(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with CSV output to a .gz file
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv.gz')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    expected_intersect_bp = [2529000, 5177000]
    with gzip.open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_intersect_bp):
            print(row)
            assert int(row['intersect_bp']) == expected


def test_prefetch_matches(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    matches_out = c.output('matches.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)

    sigs = sourmash.load_file_as_index(matches_out)

    expected_matches = [sig63, sig47]
    for (match, expected) in zip(sigs.signatures(), expected_matches):
        ss = sourmash.load_one_signature(expected, ksize=31)
        assert match == ss


def test_prefetch_matches_to_dir(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches to a directory
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches_dir/')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isdir(matches_out)

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matches_to_sig_gz(runtmp, linear_gather):
    c = runtmp

    import gzip

    # test a basic prefetch, with --save-matches to a sig.gz file
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches.sig.gz')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isfile(matches_out)

    with gzip.open(matches_out, "rt") as fp:
        # can we read this as a gz file?
        fp.read()

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matches_to_zip(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches to a zipfile
    import zipfile

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches.zip')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isfile(matches_out)

    with zipfile.ZipFile(matches_out, "r") as fp:
        # can we read this as a .zip file?
        for zi in fp.infolist():
            pass

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matching_hashes(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    matches_out = c.output('matches.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63,
                   '--save-matching-hashes', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)
    matches = set(ss47.minhash.hashes) & set(ss63.minhash.hashes)

    intersect = ss47.minhash.copy_and_clear()
    intersect.add_many(matches)

    ss = sourmash.load_one_signature(matches_out)
    assert ss.name.endswith('-known')
    assert ss.minhash == intersect


def test_prefetch_nomatch_hashes(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    nomatch_out = c.output('unmatched_hashes.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2,
                   '--save-unmatched-hashes', nomatch_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(nomatch_out)

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    remain = ss47.minhash.to_mutable()
    remain.remove_many(ss63.minhash.hashes)

    ss = sourmash.load_one_signature(nomatch_out)
    assert ss.name.endswith('-unknown')
    assert ss.minhash == remain


def test_prefetch_no_num_query(runtmp, linear_gather):
    c = runtmp

    # can't do prefetch with num signatures for query
    sig47 = utils.get_test_data('num/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig47,
                       linear_gather)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0


def test_prefetch_no_num_subj(runtmp, linear_gather):
    c = runtmp

    # can't do prefetch with num signatures for query; no matches!
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('num/63.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, linear_gather)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "ERROR in prefetch: after picklists and patterns, no signatures to search!?" in c.last_result.err


def test_prefetch_db_fromfile(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    from_file = c.output('from-list.txt')

    with open(from_file, 'wt') as fp:
        print(sig63, file=fp)
        print(sig2, file=fp)
        print(sig47, file=fp)

    c.run_sourmash('prefetch', '-k', '31', sig47, linear_gather,
                   '--db-from-file', from_file)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "query sketch has scaled=1000; will be dynamically downsampled as needed" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_no_db(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch with no databases/signatures
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('prefetch', '-k', '31', sig47, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "ERROR: no databases or signatures to search!?" in c.last_result.err


def test_prefetch_check_scaled_bounds_negative(runtmp, linear_gather):
    c = runtmp

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                    '--scaled', '-5', linear_gather)

    assert "ERROR: scaled value must be positive" in str(exc.value)


def test_prefetch_check_scaled_bounds_less_than_minimum(runtmp, linear_gather):
    c = runtmp

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                    '--scaled', '50', linear_gather)

    assert "WARNING: scaled value should be >= 100. Continuing anyway." in str(exc.value)


def test_prefetch_check_scaled_bounds_more_than_maximum(runtmp, linear_gather):
    c = runtmp

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                    '--scaled', '1e9', linear_gather)

    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in str(exc.value)


def test_prefetch_downsample_scaled(runtmp, linear_gather):
    c = runtmp

    # test --scaled
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--scaled', '1e5', linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "downsampling query from scaled=1000 to 10000" in c.last_result.err




def test_prefetch_downsample_multiple(runtmp, linear_gather):
    # test multiple different downsamplings in prefetch code
    query_sig = utils.get_test_data('GCF_000006945.2-s500.sig')

    # load in the hashes and do split them into four bins, randomly.
    ss = sourmash.load_one_signature(query_sig)
    hashes = list(ss.minhash.hashes)

    random.seed(a=1)            # fix seed so test is reproducible
    random.shuffle(hashes)

    # split into 4 bins:
    mh_bins = [ ss.minhash.copy_and_clear() for i in range(4) ]
    for i, hashval in enumerate(hashes):
        mh_bins[i % 4].add_hash(hashval)

    # downsample with different scaleds; initial scaled is 500, note.
    mh_bins[0] = mh_bins[0].downsample(scaled=750)
    mh_bins[1] = mh_bins[1].downsample(scaled=600)
    mh_bins[2] = mh_bins[2].downsample(scaled=1000)
    mh_bins[3] = mh_bins[3].downsample(scaled=650)

    gathersigs = []
    for i in range(4):
        binsig = signature.SourmashSignature(mh_bins[i], name=f"bin{i}")

        with open(runtmp.output(f"bin{i}.sig"), "wb") as fp:
            sourmash.save_signatures([binsig], fp)

        gathersigs.append(f"bin{i}.sig")

    runtmp.sourmash('prefetch', linear_gather, query_sig, *gathersigs)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "final scaled value (max across query and all matches) is 1000" in runtmp.last_result.err


def test_prefetch_empty(runtmp, linear_gather):
    c = runtmp

    # test --scaled
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                       '--scaled', '1e9', linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "no query hashes!? exiting." in c.last_result.err


def test_prefetch_basic_many_sigs(runtmp, linear_gather):
    c = runtmp

    # test what happens with many (and duplicate) signatures
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    manysigs = [sig63, sig2, sig47] * 5

    c.run_sourmash('prefetch', '-k', '31', sig47, *manysigs, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "total of 10 matching signatures so far." in c.last_result.err
    assert "total of 10 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_with_picklist(runtmp):
    # test 'sourmash prefetch' with picklists
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '--picklist', f"{picklist}:md5:md5")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 3 matches to 9 distinct values" in err
    # these are the different ksizes
    assert "WARNING: 6 missing picklist values." in err

    out = runtmp.last_result.out
    print(out)

    assert "total of 3 matching signatures." in err
    assert "of 1466 distinct query hashes, 453 were found in matches above threshold." in err
    assert "a total of 1013 query hashes remain unmatched." in err


def test_prefetch_with_picklist_exclude(runtmp):
    # test 'sourmash prefetch' with picklists, exclude
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '--picklist', f"{picklist}:md5:md5:exclude")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 9 matches by excluding 9 distinct values" in err
    # these are the different ksizes

    out = runtmp.last_result.out
    print(out)

    assert "total of 9 matching signatures." in err
    assert "of 1466 distinct query hashes, 1013 were found in matches above threshold." in err
    assert "a total of 453 query hashes remain unmatched." in err


def test_prefetch_with_pattern_include(runtmp):
    # test 'sourmash prefetch' with --include-db-pattern
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '--include', 'thermotoga')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "total of 3 matching signatures." in err
    assert "of 1466 distinct query hashes, 453 were found in matches above threshold." in err
    assert "a total of 1013 query hashes remain unmatched." in err


def test_prefetch_with_pattern_exclude(runtmp):
    # test 'sourmash prefetch' with --exclude-db-pattern
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '--exclude', 'thermotoga')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "total of 9 matching signatures." in err
    assert "of 1466 distinct query hashes, 1013 were found in matches above threshold." in err
    assert "a total of 453 query hashes remain unmatched." in err


def test_prefetch_output_with_abundance(runtmp, prefetch_gather, linear_gather):
    c = runtmp
    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against = utils.get_test_data('gather-abund/genome-s10.fa.gz.sig')

    c.run_sourmash('prefetch', linear_gather, query, against,
                   '--save-matching-hashes', c.output('match-hash.sig'),
                   '--save-unmatched-hashes', c.output('nomatch-hash.sig'))

    print(c.last_result.out)

    assert os.path.exists(c.output('match-hash.sig'))
    ss = list(sourmash.load_file_as_signatures(c.output('match-hash.sig')))[0]
    assert ss.minhash.track_abundance

    assert os.path.exists(c.output('nomatch-hash.sig'))
    ss = list(sourmash.load_file_as_signatures(c.output('nomatch-hash.sig')))[0]
    assert ss.minhash.track_abundance


def test_prefetch_ani_csv_out(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with CSV output
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    prefetch_result_names = PrefetchResult.prefetch_write_cols
    exp1 = {'q_ani': '0.9771552502238963','m_ani': '0.9767860811200507',
                      'ac_ani': '0.9769706656719734','mc_ani': '0.9771552502238963',
                      'pfn': 'False'}
    exp2 = {'q_ani': '1.0','m_ani': '1.0',
                      'ac_ani': '1.0','mc_ani': '1.0',
                      'pfn': 'False'}
    expected_ani_vals = [exp1, exp2]
    with open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_ani_vals):
            print(row)
            assert prefetch_result_names == list(row.keys())
            assert approx_eq(row['query_containment_ani'], expected['q_ani'])
            assert approx_eq(row['match_containment_ani'], expected['m_ani'])
            assert approx_eq(row['max_containment_ani'], expected['mc_ani'])
            assert approx_eq(row['average_containment_ani'], expected['ac_ani'])
            assert row['potential_false_negative'] == expected['pfn']


def test_prefetch_ani_csv_out_estimate_ci(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with CSV output
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout, linear_gather, '--estimate-ani-ci')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)
    prefetch_result_names_ci = PrefetchResult.prefetch_write_cols_ci

    exp1 = {'q_ani': '0.9771552502238963','m_ani': '0.9767860811200507',
            'q_ani_low': "0.9762537506990911", 'q_ani_high': "0.9780336875157754",
            'm_ani_low': "0.9758801604653301", "m_ani_high": "0.9776692390768575",
            'ac_ani': '0.9769706656719734','mc_ani': '0.9771552502238963',
            'pfn': 'False'}
    exp2 = {'q_ani': '1.0','m_ani': '1.0',
            'q_ani_low': "1.0", 'q_ani_high': "1.0",
            'm_ani_low': "1.0", "m_ani_high": "1.0",
                      'ac_ani': '1.0','mc_ani': '1.0',
                      'pfn': 'False'}

    expected_ani_vals = [exp1, exp2]
    with open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_ani_vals):
            print(row)
            assert prefetch_result_names_ci == list(row.keys())
            assert approx_eq(row['query_containment_ani'],expected['q_ani'])
            assert approx_eq(row['query_containment_ani_low'], expected['q_ani_low'])
            assert approx_eq(row['query_containment_ani_high'], expected['q_ani_high'])
            assert approx_eq(row['match_containment_ani'], expected['m_ani'])
            assert approx_eq(row['match_containment_ani_low'], expected['m_ani_low'])
            assert approx_eq(row['match_containment_ani_high'], expected['m_ani_high'])
            assert approx_eq(row['max_containment_ani'], expected['mc_ani'])
            assert approx_eq(row['average_containment_ani'], expected['ac_ani'])
            assert row['potential_false_negative'] == expected['pfn']


def test_prefetch_ani_containment_asymmetry(runtmp):
    # test contained_by asymmetries, viz #2215
    query_sig = utils.get_test_data('47.fa.sig')
    merged_sig = utils.get_test_data('47-63-merge.sig')

    runtmp.sourmash('prefetch', query_sig, merged_sig, '-o',
                    'query-in-merged.csv')
    runtmp.sourmash('prefetch', merged_sig, query_sig, '-o',
                    'merged-in-query.csv')

    with sourmash_args.FileInputCSV(runtmp.output('query-in-merged.csv')) as r:
        query_in_merged = list(r)[0]

    with sourmash_args.FileInputCSV(runtmp.output('merged-in-query.csv')) as r:
        merged_in_query = list(r)[0]

    assert query_in_merged['query_containment_ani'] == '1.0'
    assert query_in_merged['match_containment_ani'] == '0.9865155060423993'
    assert query_in_merged['average_containment_ani'] == '0.9932577530211997'

    assert merged_in_query['match_containment_ani'] == '1.0'
    assert merged_in_query['query_containment_ani'] == '0.9865155060423993'
    assert merged_in_query['average_containment_ani'] == '0.9932577530211997'
