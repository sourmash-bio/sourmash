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

from . import sourmash_tst_utils as utils
import sourmash_lib

from sourmash_lib.lca.lca_utils import (build_tree, find_lca)

## lca_utils tests

def test_build_tree():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def test_build_tree_3():                  # empty 'rank2' name
    tree = build_tree([[('rank1', 'name1'), ('rank2', '')]])
    assert tree == { ('rank1', 'name1'): {} }


def test_build_tree_4():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                      ])

    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                      ], tree)

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def test_find_lca():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == ((('rank1', 'name1'), ('rank2', 'name2'),), 0)


def test_find_lca_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == ((('rank1', 'name1'),), 2)


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err


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
        assert 'loaded 1 databases for LCA use.' in err


def test_single_classify_empty():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/both.lca.json')
        input_sig = utils.get_test_data('GCF_000005845.2_ASM584v2_genomic.fna.gz.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'data/GCF_000005845.2_ASM584v2_genomic.fna.gz,nomatch,,,,,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err


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
        assert 'loaded 1 databases for LCA use.' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,unassigned,Alteromonadaceae,unassigned,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err


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
        assert '...found 1 genomes with lineage assignments!!' in err
        assert '1 assigned lineages out of 1 distinct lineages in spreadsheet' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err


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
        assert '...found 2 genomes with lineage assignments!!' in err
        assert '2 assigned lineages out of 2 distinct lineages in spreadsheet' in err

        # classify input_sig1
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,unassigned,unassigned,Alteromonadaceae,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err

        # classify input_sig2
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_PSW_MAG_00136,found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 databases for LCA use.' in err


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
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 2 databases for LCA use.' in err


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

        assert 'loaded signatures from 1 files total.' in err
        assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out


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
        assert '...found 2 genomes with lineage assignments!!' in err
        assert '2 assigned lineages out of 2 distinct lineages in spreadsheet' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1,
               input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded signatures from 2 files total.' in err

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
        lines.remove('superkingdom,0')
        lines.remove('phylum,464')
        lines.remove('class,533')
        lines.remove('order,1050')
        lines.remove('family,695')
        lines.remove('genus,681')
        lines.remove('species,200')

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
        lines.remove('superkingdom,0')
        lines.remove('phylum,464')
        lines.remove('class,533')
        lines.remove('order,1050')
        lines.remove('family,695')
        lines.remove('genus,681')
        lines.remove('species,200')

        assert not lines
