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

    print(c.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "phylum,0.073,d__Bacteria;p__Bacteroidota" in  c.last_result.out
    assert "phylum,0.058,d__Bacteria;p__Proteobacteria" in  c.last_result.out
    assert "class,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in  c.last_result.out
    assert "class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in  c.last_result.out
    assert "order,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in  c.last_result.out
    assert "order,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales" in  c.last_result.out
    assert "family,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae" in  c.last_result.out
    assert "family,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae" in  c.last_result.out
    assert "genus,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia" in  c.last_result.out
    assert "genus,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella" in  c.last_result.out
    assert "genus,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola" in  c.last_result.out
    assert "species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in  c.last_result.out
    assert "species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in  c.last_result.out
    assert "species,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in  c.last_result.out


def test_summarize_csv_out(runtmp):
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = c.output('out.csv')
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '--split-identifiers', '-o', csvout)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    #expected_intersect_bp = [2529000, 5177000]
    #with open(csvout, 'rt', newline="") as fp:
    #    r = csv.DictReader(fp)
    #    for (row, expected) in zip(r, expected_intersect_bp):
    #        assert int(row['intersect_bp']) == expected


## some test ideas to start with -- see test_lca.py for add'l ideas

# test empty gather results

def test_summarize_empty_gather_results(runtmp):
    c = runtmp
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    outcsv = c.output('out.csv')
    g_csv = c.output('g.csv')
    with open(g_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_csv)

    with pytest.raises(ValueError) as exc:
        runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '--split-identifiers', '-o', outcsv)
        assert str(exc.value) == "local variable 'n' referenced before assignment"
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0



#def test_summarize_bad_gather_results():
#    pass

def test_summarize_empty_lineage_input(runtmp):
    c = runtmp
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = c.output('g.csv')
    with open(g_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_csv)

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '--split-identifiers', '-o', g_csv)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(g_csv)
    pass

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

