"""
Tests for the 'sourmash tax' command line and high level API.
"""
import os
import shutil
import csv
import pytest
import glob

import sourmash_tst_utils as utils
import sourmash
from sourmash import load_one_signature, SourmashSignature

## command line tests
def test_run_sourmash_tax():
    status, out, err = utils.runscript('sourmash', ['tax'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)

def test_summarize_stdout_0(runtmp):
    # test basic summarize
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax,
                   '--split-identifiers')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "rank,fraction,lineage" in c.last_result.out
    assert 'superkingdom,0.131,d__Bacteria' in c.last_result.out
    assert "phylum,0.073,d__Bacteria;p__Bacteroidota" in c.last_result.out
    assert "phylum,0.058,d__Bacteria;p__Proteobacteria" in c.last_result.out
    assert "class,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in c.last_result.out
    assert "class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in c.last_result.out
    assert "order,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in c.last_result.out
    assert "order,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales" in c.last_result.out
    assert "family,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae" in c.last_result.out
    assert "family,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae" in c.last_result.out
    assert "genus,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia" in c.last_result.out
    assert "genus,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella" in c.last_result.out
    assert "genus,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola" in c.last_result.out
    assert "species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out
    assert "species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in c.last_result.out
    assert "species,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in c.last_result.out


def test_summarize_summary_csv_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".summarized.csv"
    csvout = runtmp.output(sum_csv)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '--split-identifiers', '-o', csv_base)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    #expected_intersect_bp = [2529000, 5177000]
    sum_gather_results = [x.rsplit('\n') for x in open(csvout)]
    assert "rank,fraction,lineage" in sum_gather_results[0]
    assert 'superkingdom,0.131,d__Bacteria' in sum_gather_results[1]
    assert "phylum,0.073,d__Bacteria;p__Bacteroidota" in sum_gather_results[2]
    assert "phylum,0.058,d__Bacteria;p__Proteobacteria" in sum_gather_results[3]
    assert "class,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in sum_gather_results[4]
    assert "class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in sum_gather_results[5]
    assert "order,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in sum_gather_results[6]
    assert "order,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales" in sum_gather_results[7]
    assert "family,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae" in sum_gather_results[8]
    assert "family,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae" in sum_gather_results[9]
    assert "genus,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia" in sum_gather_results[10]
    assert "genus,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella" in sum_gather_results[11]
    assert "genus,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola" in sum_gather_results[12]
    assert "species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in sum_gather_results[13]
    assert "species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in sum_gather_results[14]
    assert "species,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in sum_gather_results[15]


## some test ideas to start with -- see test_lca.py for add'l ideas

#def test_summarize_empty_gather_results():
#    pass
#def test_summarize_bad_gather_results():
#    pass
#def test_summarize_empty_lineage_input():
#    pass
#def test_summarize_bad_lineage_input():
#    pass
#def test_summarize_bad_rank():
#    pass
#
#def test_classify_empty_gather_results():
#    pass
#def test_classify_bad_gather_results():
#    pass
#def test_classify_empty_lineage_input():
#    pass
#def test_classify_bad_lineage_input():
#    pass
#def test_single_classify_empty():
#    pass
#def test_mult_classify_empty():
#    pass

