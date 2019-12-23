"""
Tests for the 'sourmash lca' command line.
"""
from __future__ import print_function, unicode_literals
import os
import gzip
import shutil
import time
import screed
import glob
import json
import csv
import pytest

from . import sourmash_tst_utils as utils
import sourmash

from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import *

## lca_utils tests


def test_taxlist_1():
    assert list(taxlist()) == ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']


def test_taxlist_2():
    assert list(taxlist(include_strain=False)) == ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def test_zip_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    assert zip_lineage(x) == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    assert zip_lineage(x, truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    assert zip_lineage(x) == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    with pytest.raises(ValueError) as e:
        zip_lineage(x)

    assert 'incomplete lineage at phylum - is class instead' in str(e.value)


def test_build_tree():
    tree = build_tree([[LineagePair('rank1', 'name1'),
                        LineagePair('rank2', 'name2')]])
    assert tree == { LineagePair('rank1', 'name1'):
                         { LineagePair('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }


def test_build_tree_3():                  # empty 'rank2' name
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', '')]])
    assert tree == { LineagePair('rank1', 'name1'): {} }


def test_build_tree_4():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                      ])

    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ], tree)

    assert tree == { LineagePair('rank1', 'name1'): { LineagePair('rank2', 'name2a') : {},
                                           LineagePair('rank2', 'name2b') : {}} }

def test_build_tree_5():
    with pytest.raises(ValueError):
        tree = build_tree([])


def test_find_lca():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2'),), 0)


def test_find_lca_2():
    tree = build_tree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == ((LineagePair('rank1', 'name1'),), 2)


def test_load_single_db():
    filename = utils.get_test_data('lca/delmont-1.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    print(db)

    assert ksize == 31
    assert scaled == 10000


def test_databases():
    filename1 = utils.get_test_data('lca/delmont-1.lca.json')
    filename2 = utils.get_test_data('lca/delmont-2.lca.json')
    dblist, ksize, scaled = lca_utils.load_databases([filename1, filename2])

    print(dblist)

    assert len(dblist) == 2
    assert ksize == 31
    assert scaled == 10000


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


def test_lca_index_insert_method():
    # test 'signatures' method from base class Index
    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    sig = next(iter(db.signatures()))

    with pytest.raises(NotImplementedError) as e:
        db.insert(sig)


def test_lca_index_find_method():
    # test 'signatures' method from base class Index
    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    sig = next(iter(db.signatures()))

    with pytest.raises(NotImplementedError) as e:
        db.find(None)


def test_search_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    results = db.search(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    sig.minhash = sig.minhash.downsample_scaled(10000)
    assert sig.minhash == match_sig.minhash


def test_search_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    sig.minhash = sig.minhash.downsample_scaled(100000)

    with pytest.raises(ValueError) as e:
        results = db.search(sig, threshold=.01, ignore_abundance=True)


def test_gather_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    results = db.gather(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    sig.minhash = sig.minhash.downsample_scaled(10000)
    assert sig.minhash == match_sig.minhash


def test_gather_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    sig.minhash = sig.minhash.downsample_scaled(100000)

    results = db.gather(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    match_sig.minhash = match_sig.minhash.downsample_scaled(100000)
    assert sig.minhash == match_sig.minhash


def test_db_lineage_to_lids():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db.lineage_to_lids
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

    d = db.lid_to_idx
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, {32}), (48, {48})]


def test_db_idx_to_ident():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db.idx_to_ident
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, 'NC_009665'), (48, 'NC_011663')]


## command line tests


def test_run_sourmash_lca():
    status, out, err = utils.runscript('sourmash', ['lca'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_basic_index():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_basic_index_bad_spreadsheet():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_basic_index_broken_spreadsheet():
    # duplicate identifiers in this spreadsheet
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet-2.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert status != 0
        assert "multiple lineages for identifier TARA_ASE_MAG_00031" in err


def test_basic_index_require_taxonomy():
    # no taxonomy in here
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet-3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', '--require-taxonomy', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert status != 0
        assert "ERROR: no hash values found - are there any signatures?" in err


def test_basic_index_column_start():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', '-C', '3', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_basic_index_and_classify_with_tsv_and_gz():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.tsv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json.gz')

        cmd = ['lca', 'index', '--tabs', '--no-header', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_basic_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_index_traverse():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'index', taxcsv, lca_db, in_dir, '--traverse-directory']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_index_traverse_real_spreadsheet_no_report():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in err
        assert 'WARNING: no signatures for 956 lineage assignments.' in err
        assert 'WARNING: 105 unused lineages.' in err
        assert '(You can use --report to generate a detailed report.)' in err


def test_index_traverse_real_spreadsheet_report():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')
        report_loc = os.path.join(location, 'report.txt')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '--report',
               report_loc, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in err
        assert 'WARNING: no signatures for 956 lineage assignments.' in err
        assert 'WARNING: 105 unused lineages.' in err
        assert '(You can use --report to generate a detailed report.)' not in err
        assert os.path.exists(report_loc)


def test_single_classify():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_empty():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/both.lca.json')
        input_sig = utils.get_test_data('GCF_000005845.2_ASM584v2_genomic.fna.gz.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'data/GCF_000005845.2_ASM584v2_genomic.fna.gz,nomatch,,,,,,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_traverse():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig,
               '--traverse-directory']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_multi_query_classify_traverse():
    with utils.TempDirectory() as location:
        # both.lca.json is built from both dir and dir2
        db1 = utils.get_test_data('lca/both.lca.json')
        dir1 = utils.get_test_data('lca/dir1')
        dir2 = utils.get_test_data('lca/dir2')

        cmd = ['lca', 'classify', '--db', db1, '--query', dir1, dir2,
               '--traverse-directory']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
            fp_lines = fp.readlines()
            out_lines = out.splitlines()

            fp_lines.sort()
            out_lines.sort()

            assert len(fp_lines) == len(out_lines)
            for line1, line2 in zip(fp_lines, out_lines):
                assert line1.strip() == line2.strip(), (line1, line2)


def test_multi_db_multi_query_classify_traverse():
    with utils.TempDirectory() as location:
        # two halves of both.lca.json, see above test.
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')
        dir1 = utils.get_test_data('lca/dir1')
        dir2 = utils.get_test_data('lca/dir2')

        cmd = ['lca', 'classify', '--db', db1, db2, '--query', dir1, dir2,
               '--traverse-directory']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
            fp_lines = fp.readlines()
            out_lines = out.splitlines()

            fp_lines.sort()
            out_lines.sort()

            assert len(fp_lines) == len(out_lines)
            for line1, line2 in zip(fp_lines, out_lines):
                assert line1.strip() == line2.strip(), (line1, line2)


def test_unassigned_internal_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-4.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,unassigned,Alteromonadaceae,unassigned,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_unassigned_last_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-5.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_index_and_classify_internal_unassigned_multi():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-6.csv')
        input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        # classify input_sig1
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,unassigned,unassigned,Alteromonadaceae,,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err

        # classify input_sig2
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_PSW_MAG_00136,found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_multi_db_classify():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        db2 = utils.get_test_data('lca/delmont-2.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'classify', '--db', db1, db2, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 2 LCA databases' in err


def test_classify_unknown_hashes():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '(root)' not in out
        assert 'TARA_MED_MAG_00029,found,Archaea,Euryarcheoata,unassigned,unassigned,novelFamily_I' in out


def test_single_summarize():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out


def test_single_summarize_scaled():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,
               '--scaled', '100000']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '100.0%    27   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales'


def test_multi_summarize_with_unassigned():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-6.csv')
        input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1,
               input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 2 signatures from 2 files total.' in err

        out_lines = out.splitlines()
        out_lines.remove('14.0%   200   Bacteria')
        out_lines.remove('14.0%   200   Bacteria;Proteobacteria;unassigned;unassigned')
        out_lines.remove('86.0%  1231   Eukaryota;Chlorophyta')
        out_lines.remove('86.0%  1231   Eukaryota')
        out_lines.remove('14.0%   200   Bacteria;Proteobacteria')
        out_lines.remove('14.0%   200   Bacteria;Proteobacteria;unassigned')
        out_lines.remove('86.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae')
        out_lines.remove('14.0%   200   Bacteria;Proteobacteria;unassigned;unassigned;Alteromonadaceae')
        out_lines.remove('86.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned')
        out_lines.remove('86.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned')
        out_lines.remove('86.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned;Ostreococcus')
        assert not out_lines


def test_summarize_to_root():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '78.6%    99   Archaea' in out
        assert '21.4%    27   (root)' in out


def test_summarize_unknown_hashes():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '(root)' not in out
        assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in out


def test_rankinfo_on_multi():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')

        cmd = ['lca', 'rankinfo', db1, db2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_rankinfo_on_single():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/both.lca.json')

        cmd = ['lca', 'rankinfo', db1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_rankinfo_no_tax():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'rankinfo', lca_db]
        status, out, err = utils.runscript('sourmash', cmd)


def test_rankinfo_with_min():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')

        cmd = ['lca', 'rankinfo', db1, db2, '--minimum-num', '1']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_compare_csv():
    with utils.TempDirectory() as location:
        a = utils.get_test_data('lca/classify-by-both.csv')
        b = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')

        cmd = ['lca', 'compare_csv', a, b, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 106 distinct lineages, 957 rows' in err
        assert 'missing 937 assignments in classify spreadsheet.' in err
        assert '20 total assignments, 0 differ between spreadsheets.' in err


def test_compare_csv_real():
    with utils.TempDirectory() as location:
        a = utils.get_test_data('lca/tully-genome-sigs.classify.csv')
        b = utils.get_test_data('lca/tully-query.delmont-db.sigs.classify.csv')

        cmd = ['lca', 'compare_csv', a, b, '--start-column=3', '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 87 distinct lineages, 2631 rows' in err
        assert 'missing 71 assignments in classify spreadsheet.' in err
        assert 'missing 1380 assignments in custom spreadsheet.' in err
        assert '(these will not be evaluated any further)' in err
        assert '987 total assignments, 889 differ between spreadsheets.' in err
        assert '296 are compatible (one lineage is ancestor of another.' in err
        assert '593 are incompatible (there is a disagreement in the trees).' in err
        assert '164 incompatible at rank superkingdom' in err
        assert '255 incompatible at rank phylum' in err
        assert '107 incompatible at rank class' in err
        assert '54 incompatible at rank order' in err
        assert '13 incompatible at rank family' in err
        assert '0 incompatible at rank genus' in err
        assert '0 incompatible at rank species' in err


def test_single_gather():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'gather', input_sig, db1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '2.0 Mbp     100.0%  100.0%      Alteromonas_macleodii' in out
        assert 'Query is completely assigned.'


def test_gather_unknown_hashes():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'gather', input_sig1, lca_db]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '270.0 kbp    11.5%   21.4%      Archaea; family novelFamily_I' in out
        assert '88.5% (2.1 Mbp) of hashes have no assignment.' in out


def test_gather_combined_results():
    with utils.TempDirectory() as location:
        query_sig = utils.get_test_data('47+63.fa.sig')
        lca_db = utils.get_test_data('lca/47+63.lca.json')

        cmd = ['lca', 'gather', query_sig, lca_db, '-o', 'matches.csv']
        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        print(cmd)
        print(out)
        print(err)

        assert '5.5 Mbp      69.4%  100.0%      Shewanella baltica OS223' in out
        assert '2.4 Mbp      30.6%   47.1%      Shewanella baltica OS185' in out


def test_gather_equiv_results():
    with utils.TempDirectory() as location:
        query_sig = utils.get_test_data('47+63-intersect.fa.sig')
        lca_db = utils.get_test_data('lca/47+63.lca.json')

        cmd = ['lca', 'gather', query_sig, lca_db, '-o', 'matches.csv']
        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        print(cmd)
        print(out)
        print(err)

        assert '2.7 Mbp     100.0%' in out
        assert 'Shewanella baltica' in out
        assert '(** 1 equal matches)' in out
        assert ('OS223' in out) or ('OS185' in out)

        assert os.path.exists(lca_db)

        r = csv.DictReader(open(os.path.join(location, 'matches.csv')))
        row = next(r)
        assert row['n_equal_matches'] == '1'


def test_gather_old_lca_db():
    with utils.TempDirectory() as location:
        query_sig = utils.get_test_data('47+63.fa.sig')
        lca_db = utils.get_test_data('lca/old-db-format-1.json')

        cmd = ['lca', 'gather', query_sig, lca_db]
        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location,
                                           fail_ok=True)

        print(cmd)
        print(out)
        print(err)
        assert 'Error! This is an old-style LCA DB.' in err
        assert status != 0


@utils.in_tempdir
def test_incompat_lca_db_scaled(c):
    # create a database with scaled of 10000
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('compute', '-k', '25', '--scaled', '10000', testdata1,
                   '-o', 'test_db.sig')
    print(c)

    c.run_sourmash('lca', 'index', utils.get_test_data('lca/delmont-1.csv',),
                   'test.lca.json', 'test_db.sig',
                    '-k', '25', '--scaled', '10000')
    print(c)

    # next, create a query sig with scaled of 100000
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('compute', '-k', '25', '--scaled', '100000', testdata1,
                   '-o', 'test_query.sig')
    print(c)

    with pytest.raises(ValueError) as e:
        c.run_sourmash('lca', 'gather', 'test_query.sig', 'test.lca.json')
        print(c)

    assert 'new scaled 10000 is lower than current sample scaled 10000' in str(e.value)


@utils.in_tempdir
def test_incompat_lca_db_ksize(c):
    # create a database with ksize of 25
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('compute', '-k', '25', '--scaled', '1000', testdata1,
                   '-o', 'test_db.sig')
    print(c)

    c.run_sourmash('lca', 'index', utils.get_test_data('lca/delmont-1.csv',),
                   'test.lca.json', 'test_db.sig',
                    '-k', '25', '--scaled', '10000')
    print(c)

    # this should fail: the LCA database has ksize 25, and the query sig has
    # no compatible ksizes.
    with pytest.raises(ValueError) as e:
        c.run_sourmash('lca', 'gather', utils.get_test_data('lca/TARA_ASE_MAG_00031.sig'), 'test.lca.json')

    assert '0 signatures matching ksize and molecule type;' in str(e.value)
