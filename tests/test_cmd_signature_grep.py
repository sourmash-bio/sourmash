"""
Tests for the 'sourmash signature grep' command line.
"""
import shutil
import os
import csv
import gzip

import pytest

import sourmash_tst_utils as utils
import sourmash
from sourmash_tst_utils import SourmashCommandFailed
from sourmash.signature import load_signatures

## command line tests


def test_grep_1_sig_name(runtmp):
    # search on substring in name
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', 'Shewanella', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella' in ss.name
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_1_sig_name_case_sensitive(runtmp):
    # search on substring in name
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', 'shewanella', sig47)


def test_grep_1_sig_name_case_insensitive(runtmp):
    # search on substring in name, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '-i', 'shewanella', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella' in ss.name
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_1_sig_name_exclude(runtmp):
    # search on substring in name, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    # no matches!
    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', '-v', 'Shewanella', sig47)


def test_grep_2_sig_md5(runtmp):
    # search on substring in md5
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', 'ce52952152f0', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_2_sig_md5_case_sensitive(runtmp):
    # case sensitive no match
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('sig', 'grep', 'CE52952152f0', sig47)


def test_grep_2_sig_md5_case_insensitive(runtmp):
    # search on substring in md5, case insensitive
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '-i', 'CE52952152f0', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_3_filename(runtmp):
    # filename match
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '47.fa', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert '47.fa' in ss.filename
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_3_filename_regexp(runtmp):
    # search for a regexp on filename
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.run_sourmash('sig', 'grep', '^47.fa', sig47)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert '7.fa' in ss.filename
    assert ss.md5sum() == '09a08691ce52952152f0e866a59f6261'


def test_grep_4_no_manifest(runtmp):
    # fail search when no manifest, by default
    sbt = utils.get_test_data('v6.sbt.zip')

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('sig', 'grep', 'e60265', sbt)

    print(runtmp.last_result.err)
    assert 'ERROR on filename' in runtmp.last_result.err
    assert 'sig grep requires a manifest by default, but no manifest present.' in runtmp.last_result.err


def test_grep_4_no_manifest_ok(runtmp):
    # generate manifest if --no-require-manifest
    sbt = utils.get_test_data('v6.sbt.zip')

    runtmp.run_sourmash('sig', 'grep', 'e60265', sbt, '--no-require-manifest')

    ss = load_signatures(runtmp.last_result.out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'e60265' in ss.md5sum()


def test_grep_5_zip_include(runtmp):
    # search zip, include on case sensitive match to name
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_include_picklist(runtmp):
    # search zip, include on case sensitive match to name
    allzip = utils.get_test_data('prot/all.zip')

    pickfile = runtmp.output('pick.csv')
    with open(pickfile, 'w', newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=['md5'])
        w.writeheader()
        w.writerow(dict(md5='09a08691ce52952152f0e866a59f6261'))
        w.writerow(dict(md5='38729c6374925585db28916b82a6f513'))

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip,
                        '--picklist', f"{pickfile}:md5:md5")

    out = runtmp.last_result.out
    print(out)
    err = runtmp.last_result.err
    print(err)
    assert 'for given picklist, found 2 matches to 2 distinct values' in err

    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_include_case_insensitive(runtmp):
    # search zip, include on case insensitive match to name
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-i', 'os223', allzip)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_exclude(runtmp):
    # search zip, exclude on case-sensitive match
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-v', 'OS185', allzip)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_5_zip_exclude_case_insensitive(runtmp):
    # search zip, exclude on case-insensitive match
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', '-vi', 'os185', allzip)

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_6_zip_manifest_csv(runtmp):
    # do --csv and use result as picklist
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip,
                        '--csv', 'match.csv')

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'

    # now run cat with picklist
    runtmp.run_sourmash('sig', 'cat', allzip,
                        '--picklist', 'match.csv::manifest')

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_grep_6_zip_manifest_csv_gz(runtmp):
    # do --csv and use result as picklist
    allzip = utils.get_test_data('prot/all.zip')

    runtmp.run_sourmash('sig', 'grep', '--dna', 'OS223', allzip,
                        '--csv', 'match.csv.gz')

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'

    # check that match.csv.gz is a gzip file
    with gzip.open(runtmp.output('match.csv.gz'), 'rt', newline='') as fp:
        fp.read()

    # now run cat with picklist
    runtmp.run_sourmash('sig', 'cat', allzip,
                        '--picklist', 'match.csv.gz::manifest')

    out = runtmp.last_result.out
    ss = load_signatures(out)
    ss = list(ss)
    assert len(ss) == 1
    ss = ss[0]
    assert 'Shewanella baltica OS223' in ss.name
    assert ss.md5sum() == '38729c6374925585db28916b82a6f513'


def test_sig_grep_7_lca(runtmp):
    # extract 47 from an LCA database, with --no-require-manifest
    allzip = utils.get_test_data('lca/47+63.lca.json')
    sig47 = utils.get_test_data('47.fa.sig')

    runtmp.sourmash('sig', 'grep', "50a9274021e4", allzip,
                    '--no-require-manifest', '-o', 'matches.sig')

    match = sourmash.load_file_as_signatures(runtmp.output('matches.sig'))
    match = list(match)[0]

    ss47 = sourmash.load_file_as_signatures(sig47)
    ss47 = list(ss47)[0]

    ss47 = ss47.to_mutable()
    ss47.minhash = ss47.minhash.downsample(scaled=10000)

    assert ss47.minhash == match.minhash


def test_sig_grep_7_picklist_md5_lca_fail(runtmp):
    # extract 47 from an LCA database, using a picklist w/full md5 => fail
    allzip = utils.get_test_data('lca/47+63.lca.json')

    # select on any of these attributes
    row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
               md5full='50a9274021e43eda8b2e77f8fa60ae8e',
               md5short='50a9274021e43eda8b2e77f8fa60ae8e'[:8],
               fullIdent='NC_009665.1',
               nodotIdent='NC_009665')

    # make picklist
    picklist_csv = runtmp.output('pick.csv')
    with open(picklist_csv, 'w', newline='') as csvfp:
        w = csv.DictWriter(csvfp, fieldnames=row.keys())
        w.writeheader()
        w.writerow(row)

    picklist_arg = f"{picklist_csv}:md5full:md5"
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash('sig', 'grep', '50a92740', allzip,
                        '--picklist', picklist_arg,
                        '--no-require-manifest')

    # this happens b/c the implementation of 'grep' uses picklists, and
    # LCA databases don't support multiple picklists.
    print(runtmp.last_result.err)
    assert "This input collection doesn't support 'grep' with picklists." in runtmp.last_result.err


def test_sig_grep_8_count(runtmp):
    zips = ['prot/all.zip',
            'prot/dayhoff.sbt.zip',
            'prot/dayhoff.zip',
            'prot/hp.sbt.zip',
            'prot/hp.zip',
            'prot/protein.sbt.zip',
            'prot/protein.zip']

    zip_src = [ utils.get_test_data(x) for x in zips ]

    os.mkdir(runtmp.output('prot'))
    for src, dest in zip(zip_src, zips):
        shutil.copyfile(src, runtmp.output(dest))
    
    runtmp.sourmash('sig', 'grep', '-c', '0015939', *zips)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    assert "(no signatures will be saved because of --silent/--count)." in err

    for line in """\
6 matches: prot/all.zip
2 matches: prot/dayhoff.sbt.zip
2 matches: prot/dayhoff.zip
2 matches: prot/hp.sbt.zip
2 matches: prot/hp.zip
2 matches: prot/protein.sbt.zip
2 matches: prot/protein.zip
""".splitlines():
        assert line.strip() in out
