"""
Tests for the 'sourmash tax' command line and high level API.
"""
import os
import csv
import pytest
import gzip
from collections import Counter

import sourmash
import sourmash_tst_utils as utils
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash_tst_utils import SourmashCommandFailed

from sourmash import sqlite_utils
from sourmash.exceptions import IndexNotSupported
from sourmash import sourmash_args

## command line tests
def test_run_sourmash_tax():
    status, out, err = utils.runscript('sourmash', ['tax'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_metagenome_stdout_0(runtmp):
    # test basic metagenome
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out
    assert 'test1,superkingdom,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,phylum,0.116,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,phylum,0.088,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,phylum,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,class,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,class,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,class,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,order,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,order,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,order,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,family,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,family,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,genus,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella,md5,test1.sig,0.057,444000' in c.last_result.out
    assert 'test1,genus,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,genus,0.028,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola,md5,test1.sig,0.016,138000' in c.last_result.out
    assert 'test1,genus,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000' in c.last_result.out
    assert 'test1,species,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,species,0.028,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus,md5,test1.sig,0.016,138000' in c.last_result.out
    assert 'test1,species,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out


def test_metagenome_stdout_0_db(runtmp):
    # test basic metagenome with sqlite database
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.db')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out
    assert 'test1,superkingdom,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,phylum,0.116,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,phylum,0.088,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,phylum,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,class,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,class,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,class,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,order,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,order,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,order,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,family,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,family,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,genus,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella,md5,test1.sig,0.057,444000' in c.last_result.out
    assert 'test1,genus,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,genus,0.028,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola,md5,test1.sig,0.016,138000' in c.last_result.out
    assert 'test1,genus,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000' in c.last_result.out
    assert 'test1,species,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,species,0.028,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus,md5,test1.sig,0.016,138000' in c.last_result.out
    assert 'test1,species,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out


def test_metagenome_summary_csv_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".summarized.csv"
    csvout = runtmp.output(sum_csv)
    outdir = os.path.dirname(csvout)

    runtmp.run_sourmash('tax', 'metagenome', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-dir', outdir)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    sum_gather_results = [x.rstrip() for x in open(csvout)]
    assert f"saving 'csv_summary' output to '{csvout}'" in runtmp.last_result.err
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in sum_gather_results[0]
    assert 'test1,superkingdom,0.2042281611487834,d__Bacteria,md5,test1.sig,0.13080306238801107,1024000' in  sum_gather_results[1]
    assert 'test1,superkingdom,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[2]
    assert 'test1,phylum,0.11607499002792182,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.07265026877341586,582000' in  sum_gather_results[3]
    assert 'test1,phylum,0.08815317112086159,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[4]
    assert 'test1,phylum,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[5]
    assert 'test1,class,0.11607499002792182,d__Bacteria;p__Bacteroidota;c__Bacteroidia,md5,test1.sig,0.07265026877341586,582000' in sum_gather_results[6]
    assert 'test1,class,0.08815317112086159,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[7]
    assert 'test1,class,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[8]
    assert 'test1,order,0.11607499002792182,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales,md5,test1.sig,0.07265026877341586,582000' in sum_gather_results[9]
    assert 'test1,order,0.08815317112086159,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[10]
    assert 'test1,order,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[11]
    assert 'test1,family,0.11607499002792182,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.07265026877341586,582000' in sum_gather_results[12]
    assert 'test1,family,0.08815317112086159,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[13]
    assert 'test1,family,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[14]
    assert 'test1,genus,0.0885520542481053,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella,md5,test1.sig,0.05701254275940707,444000' in sum_gather_results[15]
    assert 'test1,genus,0.08815317112086159,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[16]
    assert 'test1,genus,0.027522935779816515,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola,md5,test1.sig,0.015637726014008795,138000' in sum_gather_results[17]
    assert 'test1,genus,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[18]
    assert 'test1,species,0.0885520542481053,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.05701254275940707,444000' in sum_gather_results[19]
    assert 'test1,species,0.08815317112086159,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli,md5,test1.sig,0.05815279361459521,442000' in sum_gather_results[20]
    assert 'test1,species,0.027522935779816515,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus,md5,test1.sig,0.015637726014008795,138000' in sum_gather_results[21]
    assert 'test1,species,0.7957718388512166,unclassified,md5,test1.sig,0.8691969376119889,3990000' in sum_gather_results[22]


def test_metagenome_summary_csv_out_empty_gather_force(runtmp):
    # test multiple -g, empty -g file, and --force
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".summarized.csv"
    csvout = runtmp.output(sum_csv)
    outdir = os.path.dirname(csvout)

    gather_empty = runtmp.output('g.csv')
    with open(gather_empty, "w") as fp:
        fp.write("")
    print("g_csv: ", gather_empty)

    runtmp.run_sourmash('tax', 'metagenome', '--gather-csv', g_csv, '-g', gather_empty, '--taxonomy-csv', tax, '-o', csv_base, '--output-dir', outdir, '-f')
    sum_gather_results = [x.rstrip() for x in open(csvout)]
    assert f"saving 'csv_summary' output to '{csvout}'" in runtmp.last_result.err
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in sum_gather_results[0]
    assert 'test1,superkingdom,0.2042281611487834,d__Bacteria,md5,test1.sig,0.13080306238801107,1024000' in  sum_gather_results[1]


def test_metagenome_kreport_out(runtmp):
    # test 'kreport' kraken output format
    g_csv = utils.get_test_data('tax/test1.gather.v450.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".kreport.txt"
    csvout = runtmp.output(sum_csv)
    outdir = os.path.dirname(csvout)

    runtmp.run_sourmash('tax', 'metagenome', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-dir', outdir, '-F', "kreport")

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    kreport_results = [x.rstrip().split('\t') for x in open(csvout)]
    assert f"saving 'kreport' output to '{csvout}'" in runtmp.last_result.err
    print(kreport_results)
    assert ['13.08', '1605999', '0', 'D', '', 'd__Bacteria'] == kreport_results[0]
    assert ['86.92', '10672000', '10672000', 'U', '', 'unclassified'] == kreport_results[1]
    assert ['7.27', '892000', '0', 'P', '', 'p__Bacteroidota'] == kreport_results[2]
    assert ['5.82', '714000', '0', 'P', '', 'p__Proteobacteria'] == kreport_results[3]
    assert ['7.27', '892000', '0', 'C', '', 'c__Bacteroidia'] == kreport_results[4]
    assert ['5.82', '714000', '0', 'C', '', 'c__Gammaproteobacteria'] == kreport_results[5]
    assert ['7.27', '892000', '0', 'O', '', 'o__Bacteroidales'] == kreport_results[6]
    assert ['5.82', '714000', '0', 'O', '', 'o__Enterobacterales'] == kreport_results[7]
    assert ['7.27', '892000', '0', 'F', '', 'f__Bacteroidaceae'] == kreport_results[8]
    assert ['5.82', '714000', '0', 'F', '', 'f__Enterobacteriaceae'] == kreport_results[9]
    assert ['5.70', '700000', '0', 'G', '', 'g__Prevotella']  == kreport_results[10]
    assert ['5.82', '714000', '0', 'G', '', 'g__Escherichia'] == kreport_results[11]
    assert ['1.56', '192000', '0', 'G', '', 'g__Phocaeicola'] == kreport_results[12]
    assert ['5.70', '700000', '700000', 'S', '', 's__Prevotella copri'] == kreport_results[13]
    assert ['5.82', '714000', '714000', 'S', '', 's__Escherichia coli']== kreport_results[14]
    assert ['1.56', '192000', '192000', 'S', '', 's__Phocaeicola vulgatus'] == kreport_results[15]


def test_metagenome_kreport_out_lemonade(runtmp):
    # test 'kreport' kraken output format against lemonade output
    g_csv = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.csv')
    tax = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.tax.csv')
    csv_base = "out"
    sum_csv = csv_base + ".kreport.txt"
    csvout = runtmp.output(sum_csv)
    outdir = os.path.dirname(csvout)

    runtmp.run_sourmash('tax', 'metagenome', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-dir', outdir, '-F', "kreport")

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    kreport_results = [x.rstrip().split('\t') for x in open(csvout)]
    assert f"saving 'kreport' output to '{csvout}'" in runtmp.last_result.err
    print(kreport_results)
    assert ['5.35', '116000', '0', 'D', '', 'd__Bacteria'] == kreport_results[0]
    assert ['94.65', '2054000', '2054000', 'U', '', 'unclassified'] == kreport_results[1]
    assert ['5.35', '116000', '0', 'P', '', 'p__Bacteroidota'] == kreport_results[2]
    assert ['5.35', '116000', '0', 'C', '', 'c__Chlorobia'] == kreport_results[3]
    assert ['5.35', '116000', '0', 'O', '', 'o__Chlorobiales'] == kreport_results[4]
    assert ['5.35', '116000', '0', 'F', '', 'f__Chlorobiaceae'] == kreport_results[5]
    assert ['5.35', '116000', '0', 'G', '', 'g__Prosthecochloris'] == kreport_results[6]
    assert ['5.35', '116000', '116000', 'S', '', 's__Prosthecochloris vibrioformis'] == kreport_results[7]


def test_metagenome_kreport_out_fail(runtmp):
    # kreport cannot be generated with gather results from < v4.5.0
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".kreport.txt"
    csvout = runtmp.output(sum_csv)
    outdir = os.path.dirname(csvout)

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('tax', 'metagenome', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-dir', outdir, '-F', "kreport")

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "ERROR: cannot produce 'kreport' format from gather results before sourmash v4.5.0" in runtmp.last_result.err


def test_metagenome_krona_tsv_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    kr_csv = csv_base + ".krona.tsv"
    csvout = runtmp.output(kr_csv)
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax, '-o', csv_base,
                        '--output-format', 'krona', '--rank', 'genus', '--output-dir', outdir)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)
    assert f"saving 'krona' output to '{csvout}'" in runtmp.last_result.err

    gn_krona_results = [x.rstrip().split('\t') for x in open(csvout)]
    print("species krona results: \n", gn_krona_results)
    assert ['fraction', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus'] == gn_krona_results[0]
    assert ['0.0885520542481053', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Prevotella']  == gn_krona_results[1]
    assert ['0.08815317112086159', 'd__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacterales', 'f__Enterobacteriaceae', 'g__Escherichia'] == gn_krona_results[2]
    assert ['0.027522935779816515', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Phocaeicola'] == gn_krona_results[3]
    assert ['0.7957718388512166', 'unclassified', 'unclassified', 'unclassified', 'unclassified', 'unclassified', 'unclassified'] == gn_krona_results[4]


def test_metagenome_lineage_summary_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    lin_csv = csv_base + ".lineage_summary.tsv"
    csvout = runtmp.output(lin_csv)
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax,
                        '-o', csv_base, '--output-format', 'lineage_summary', '--rank',
                        'genus', '--output-dir', outdir)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)
    assert f"saving 'lineage_summary' output to '{csvout}'" in runtmp.last_result.err

    gn_lineage_summary = [x.rstrip().split('\t') for x in open(csvout)]
    print("species lineage summary results: \n", gn_lineage_summary)
    assert ['lineage', 'test1'] == gn_lineage_summary[0]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola', '0.027522935779816515'] == gn_lineage_summary[1]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella', '0.0885520542481053'] == gn_lineage_summary[2]
    assert ['d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia', '0.08815317112086159'] == gn_lineage_summary[3]
    assert ['unclassified', '0.7957718388512166']  == gn_lineage_summary[4]


def test_metagenome_human_format_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    csvout = runtmp.output(csv_base + '.human.txt')
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax,
                        '-o', csv_base, '--output-format', 'human', '--rank',
                        'genus', '--output-dir', outdir)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)
    assert f"saving 'human' output to '{csvout}'" in runtmp.last_result.err

    with open(csvout) as fp:
        outp = fp.readlines()

    assert len(outp) == 6
    outp = [ x.strip() for x in outp ]

    assert outp[0] == 'sample name    proportion   lineage'
    assert outp[1] == '-----------    ----------   -------'
    assert outp[2] == 'test1             86.9%     unclassified'
    assert outp[3] == 'test1              5.8%     d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia'
    assert outp[4] == 'test1              5.7%     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella'
    assert outp[5] == 'test1              1.6%     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola'


def test_metagenome_no_taxonomy_fail(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome', '-g', g_csv)
    assert "error: the following arguments are required: -t/--taxonomy-csv" in str(exc.value)


def test_metagenome_no_rank_lineage_summary(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'lineage_summary')
    assert "Rank (--rank) is required for krona and lineage_summary output formats." in str(exc.value)


def test_metagenome_no_rank_krona(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'krona')
    assert "Rank (--rank) is required for krona and lineage_summary output formats." in str(exc.value)


def test_genome_no_rank_krona(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'krona')
    assert "Rank (--rank) is required for krona output format." in str(exc.value)


def test_metagenome_rank_not_available(runtmp):
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--rank', 'strain')

    print(str(exc.value))

    assert c.last_result.status == -1
    assert "No taxonomic information provided for rank strain: cannot summarize at this rank" in str(exc.value)


def test_metagenome_duplicated_taxonomy_fail(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1] + 'FOO') # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', duplicated_csv)

    assert "cannot read taxonomy" in str(exc.value)
    assert "multiple lineages for identifier GCF_001881345" in str(exc.value)


def test_metagenome_duplicated_taxonomy_force(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', duplicated_csv, '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    # same as stdout test - just check the first few lines
    assert c.last_result.status == 0
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out
    assert 'test1,superkingdom,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,phylum,0.116,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,phylum,0.088,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,phylum,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out


def test_metagenome_missing_taxonomy(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        subset.write("\n".join(tax[:4]))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', subset_csv)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_003471795" in c.last_result.err

    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.193,d__Bacteria,md5,test1.sig,0.124,970000'in c.last_result.out
    assert 'test1,superkingdom,0.807,unclassified,md5,test1.sig,0.876,4044000' in c.last_result.out
    assert 'test1,phylum,0.105,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.066,528000' in c.last_result.out
    assert 'test1,phylum,0.088,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,phylum,0.807,unclassified,md5,test1.sig,0.876,4044000' in c.last_result.out
    assert 'test1,class,0.105,d__Bacteria;p__Bacteroidota;c__Bacteroidia,md5,test1.sig,0.066,528000' in c.last_result.out


def test_metagenome_missing_fail_taxonomy(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        subset.write("\n".join(tax[:4]))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', subset_csv, '--fail-on-missing-taxonomy')

    print(str(exc.value))

    assert "The following are missing from the taxonomy information: GCF_003471795" in str(exc.value)
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in str(exc.value)
    assert c.last_result.status == -1


def test_metagenome_multiple_taxonomy_files_missing(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # gather against mult databases
    g_csv = utils.get_test_data('tax/test1_x_gtdbrs202_genbank_euks.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', taxonomy_csv, '--force')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert "of 6 gather results, missed 2 lineage assignments." in c.last_result.err
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'multtest,superkingdom,0.204,d__Bacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.131,1024000' in c.last_result.out
    assert 'multtest,superkingdom,0.796,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.869,3990000' in c.last_result.out
    assert 'multtest,phylum,0.116,d__Bacteria;p__Bacteroidota,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out
    assert 'multtest,phylum,0.088,d__Bacteria;p__Proteobacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.058,442000' in c.last_result.out
    assert 'multtest,phylum,0.796,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.869,3990000' in c.last_result.out
    assert 'multtest,class,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out
    assert 'multtest,class,0.088,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.058,442000' in c.last_result.out
    assert 'multtest,class,0.796,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.869,3990000' in c.last_result.out


def test_metagenome_multiple_taxonomy_files(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    protozoa_genbank = utils.get_test_data('tax/protozoa_genbank_lineage.csv')
    bacteria_refseq  = utils.get_test_data('tax/bacteria_refseq_lineage.csv')

    # gather against mult databases
    g_csv = utils.get_test_data('tax/test1_x_gtdbrs202_genbank_euks.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', taxonomy_csv, protozoa_genbank, bacteria_refseq)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'multtest,superkingdom,0.204,Bacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.131,1024000' in c.last_result.out
    assert 'multtest,superkingdom,0.051,Eukaryota,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,superkingdom,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,phylum,0.116,Bacteria;Bacteroidetes,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out
    assert 'multtest,phylum,0.088,Bacteria;Proteobacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.058,442000' in c.last_result.out
    assert 'multtest,phylum,0.051,Eukaryota;Apicomplexa,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,phylum,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,class,0.116,Bacteria;Bacteroidetes;Bacteroidia,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out


def test_metagenome_multiple_taxonomy_files_multiple_taxonomy_args(runtmp):
    c = runtmp
    # pass in mult tax files using mult tax arguments
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    protozoa_genbank = utils.get_test_data('tax/protozoa_genbank_lineage.csv')
    bacteria_refseq  = utils.get_test_data('tax/bacteria_refseq_lineage.csv')

    # gather against mult databases
    g_csv = utils.get_test_data('tax/test1_x_gtdbrs202_genbank_euks.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', taxonomy_csv, '-t', protozoa_genbank, '-t', bacteria_refseq)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'multtest,superkingdom,0.204,Bacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.131,1024000' in c.last_result.out
    assert 'multtest,superkingdom,0.051,Eukaryota,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,superkingdom,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,phylum,0.116,Bacteria;Bacteroidetes,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out
    assert 'multtest,phylum,0.088,Bacteria;Proteobacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.058,442000' in c.last_result.out
    assert 'multtest,phylum,0.051,Eukaryota;Apicomplexa,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,phylum,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,class,0.116,Bacteria;Bacteroidetes;Bacteroidia,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out


def test_metagenome_multiple_taxonomy_files_multiple_taxonomy_args_empty_force(runtmp):
    # pass in mult tax files using mult tax arguments, with one empty,
    # and use --force
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    protozoa_genbank = utils.get_test_data('tax/protozoa_genbank_lineage.csv')
    bacteria_refseq  = utils.get_test_data('tax/bacteria_refseq_lineage.csv')

    tax_empty = runtmp.output('t.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    # gather against mult databases
    g_csv = utils.get_test_data('tax/test1_x_gtdbrs202_genbank_euks.gather.csv')

    c.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', taxonomy_csv, '-t', protozoa_genbank, '-t', bacteria_refseq, '-t', tax_empty, '--force')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'multtest,superkingdom,0.204,Bacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.131,1024000' in c.last_result.out
    assert 'multtest,superkingdom,0.051,Eukaryota,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,superkingdom,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,phylum,0.116,Bacteria;Bacteroidetes,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out
    assert 'multtest,phylum,0.088,Bacteria;Proteobacteria,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.058,442000' in c.last_result.out
    assert 'multtest,phylum,0.051,Eukaryota;Apicomplexa,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.245,258000' in c.last_result.out
    assert 'multtest,phylum,0.744,unclassified,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.624,3732000' in c.last_result.out
    assert 'multtest,class,0.116,Bacteria;Bacteroidetes;Bacteroidia,9687eeed,outputs/abundtrim/HSMA33MX.abundtrim.fq.gz,0.073,582000' in c.last_result.out


def test_metagenome_empty_gather_results(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    #creates empty gather result
    g_csv = runtmp.output('g.csv')
    with open(g_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_csv)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax)

    assert f"Cannot read gather results from '{g_csv}'. Is file empty?" in str(exc.value)
    assert runtmp.last_result.status == -1


def test_metagenome_bad_gather_header(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    bad_g_csv = runtmp.output('g.csv')

    #creates bad gather result
    bad_g = [x.replace("name", "nope") for x in open(g_csv, 'r')]
    with open(bad_g_csv, 'w') as fp:
        for line in bad_g:
            fp.write(line)
    print("bad_gather_results: \n", bad_g)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', bad_g_csv, '--taxonomy-csv', tax)

    assert f"Not all required gather columns are present in '{bad_g_csv}'." in str(exc.value)
    assert runtmp.last_result.status == -1


def test_metagenome_empty_tax_lineage_input(runtmp):
    # test an empty tax CSV
    tax_empty = runtmp.output('t.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)


    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax_empty)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status != 0
    assert "cannot read taxonomy assignments from" in str(exc.value)


def test_metagenome_empty_tax_lineage_input_force(runtmp):
    # test an empty tax CSV with --force
    tax_empty = runtmp.output('t.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)


    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-g', g_csv, '--taxonomy-csv', tax_empty, '--force')

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status != 0
    assert "ERROR: No taxonomic assignments loaded" in str(exc.value)


def test_metagenome_perfect_match_warning(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    perfect_g_csv = runtmp.output('g.csv')

    #create a perfect gather result
    with open(g_csv, 'r') as fp:
        r = csv.DictReader(fp, delimiter=',')
        header = r.fieldnames
        print(header)
        with open(perfect_g_csv, 'w') as out_fp:
            w = csv.DictWriter(out_fp, header)
            w.writeheader()
            for n, row in enumerate(r):
                if n == 0:
                    # make a perfect match
                    row["f_unique_to_query"] = 1.0
                else:
                    # set the rest to 0
                    row["f_unique_to_query"] = 0.0
                w.writerow(row)
                print(row)

    runtmp.run_sourmash('tax', 'metagenome', '-g', perfect_g_csv, '--taxonomy-csv', tax)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert 'WARNING: 100% match! Is query "test1" identical to its database match, GCF_001881345' in runtmp.last_result.err


def test_metagenome_over100percent_error(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    perfect_g_csv = runtmp.output('g.csv')

    #create a perfect gather result
    with open(g_csv, 'r') as fp:
        r = csv.DictReader(fp, delimiter=',')
        header = r.fieldnames
        print(header)
        with open(perfect_g_csv, 'w') as out_fp:
            w = csv.DictWriter(out_fp, header)
            w.writeheader()
            for n, row in enumerate(r):
                if n == 0:
                    row["f_unique_to_query"] = 1.0
                # let the rest stay as they are (should be > 100% match now)
                w.writerow(row)
                print(row)

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('tax', 'metagenome', '-g', perfect_g_csv, '--taxonomy-csv', tax)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == -1
    assert "ERROR: The tax summary of query 'test1' is 1.1160749900279219, which is > 100% of the query!!" in runtmp.last_result.err


def test_metagenome_gather_duplicate_query(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # different filename, contents identical to test1
    g_res2 = runtmp.output("test2.gather.csv")
    with open(g_res2, 'w') as fp:
        for line in open(g_res, 'r'):
            fp.write(line)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome',  '--gather-csv', g_res, g_res2,
                   '--taxonomy-csv', taxonomy_csv)

    assert c.last_result.status == -1
    print(str(exc.value))
    assert "Gather query test1 was found in more than one CSV. Cannot load from " in str(exc.value)


def test_metagenome_gather_duplicate_query_force(runtmp):
    # do not load same query from multiple files.
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # different filename, contents identical to test1
    g_res2 = runtmp.output("test2.gather.csv")
    with open(g_res2, 'w') as fp:
        for line in open(g_res, 'r'):
            fp.write(line)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'metagenome',  '--gather-csv', g_res, g_res2,
                   '--taxonomy-csv', taxonomy_csv, '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1

    assert "Gather query test1 was found in more than one CSV." in c.last_result.err
    assert "Cannot force past duplicated gather query. Exiting." in c.last_result.err


def test_metagenome_gather_duplicate_filename(runtmp):
    # test that a duplicate filename is properly flagged, when passed in
    # twice to a single -g argument.
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'metagenome', '--gather-csv', g_res, g_res, '--taxonomy-csv', taxonomy_csv)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}'
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out


def test_metagenome_gather_duplicate_filename_2(runtmp):
    # test that a duplicate filename is properly flagged, with -g a -g b
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'metagenome', '--gather-csv', g_res, '-g', g_res, '--taxonomy-csv', taxonomy_csv)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}'
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out


def test_metagenome_gather_duplicate_filename_from_file(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'metagenome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}'
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out


def test_genome_empty_gather_results(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    #creates empty gather result
    g_csv = runtmp.output('g.csv')
    with open(g_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_csv)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax)

    assert f"Cannot read gather results from '{g_csv}'. Is file empty?" in str(exc.value)
    assert runtmp.last_result.status == -1


def test_genome_bad_gather_header(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    bad_g_csv = runtmp.output('g.csv')

    #creates bad gather result
    bad_g = [x.replace("f_unique_to_query", "nope") for x in open(g_csv, 'r')]
    with open(bad_g_csv, 'w') as fp:
        for line in bad_g:
            fp.write(line)
    print("bad_gather_results: \n", bad_g)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'genome', '-g', bad_g_csv, '--taxonomy-csv', tax)

    assert f"Not all required gather columns are present in '{bad_g_csv}'." in str(exc.value)
    assert runtmp.last_result.status == -1


def test_genome_empty_tax_lineage_input(runtmp):
    # test an empty tax csv
    tax_empty = runtmp.output('t.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax_empty)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status != 0
    assert "cannot read taxonomy assignments from" in str(exc.value)


def test_genome_rank_stdout_0(runtmp):
    # test basic genome
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'genome', '--gather-csv', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0'  in c.last_result.out


def test_genome_rank_stdout_0_db(runtmp):
    # test basic genome with sqlite database
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.db')

    c.run_sourmash('tax', 'genome', '--gather-csv', g_csv, '--taxonomy-csv',
                   tax, '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0'  in c.last_result.out

    # too stringent of containment threshold:
    c.run_sourmash('tax', 'genome', '--gather-csv', g_csv, '--taxonomy-csv',
                   tax, '--rank', 'species', '--containment-threshold', '1.0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "WARNING: classifying query test1 at desired rank species does not meet containment threshold 1.0" in c.last_result.err
    assert "test1,below_threshold,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0," in c.last_result.out


def test_genome_rank_csv_0(runtmp):
    # test basic genome - output csv
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    cl_csv = csv_base + ".classifications.csv"
    csvout = runtmp.output(cl_csv)
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '-o', csv_base, '--containment-threshold', '0',
                   '--output-dir', outdir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert f"saving 'classification' output to '{csvout}'" in runtmp.last_result.err
    assert c.last_result.status == 0
    cl_results = [x.rstrip() for x in open(csvout)]
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in cl_results[0]
    assert 'test1,match,species,0.0885520542481053,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.05701254275940707,444000.0' in cl_results[1]


def test_genome_rank_krona(runtmp):
    # test basic genome - output csv
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    cl_csv = csv_base + ".krona.tsv"
    csvout = runtmp.output(cl_csv)
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '-o', csv_base, '--containment-threshold', '0',
                   '--output-format', 'krona', '--output-dir', outdir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert f"saving 'krona' output to '{csvout}'" in runtmp.last_result.err
    assert c.last_result.status == 0
    kr_results = [x.rstrip().split('\t') for x in open(csvout)]
    print(kr_results)
    assert ['fraction', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']  == kr_results[0]
    assert ['0.0885520542481053', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Prevotella', 's__Prevotella copri'] == kr_results[1]


def test_genome_rank_human_output(runtmp):
    # test basic genome - output csv
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    csvout = runtmp.output(csv_base + '.human.txt')
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '-o', csv_base, '--containment-threshold', '0',
                   '--output-format', 'human', '--output-dir', outdir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert f"saving 'human' output to '{csvout}'" in runtmp.last_result.err
    assert c.last_result.status == 0

    with open(csvout) as fp:
        outp = fp.readlines()

    assert len(outp) == 3
    outp = [ x.strip() for x in outp ]

    assert outp[0] == 'sample name    proportion   lineage'
    assert outp[1] == '-----------    ----------   -------'
    assert outp[2] == 'test1              5.7%     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri'


def test_genome_rank_lineage_csv_output(runtmp):
    # test basic genome - output csv
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    csvout = runtmp.output(csv_base + '.lineage.csv')
    outdir = os.path.dirname(csvout)
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '-o', csv_base, '--containment-threshold', '0',
                   '--output-format', 'lineage_csv', '--output-dir', outdir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert f"saving 'lineage_csv' output to '{csvout}'" in runtmp.last_result.err
    assert c.last_result.status == 0
    with open(csvout) as fp:
        outp = fp.readlines()

    assert len(outp) == 2
    outp = [ x.strip() for x in outp ]

    assert outp[0] == 'ident,superkingdom,phylum,class,order,family,genus,species'
    assert outp[1] == 'test1,d__Bacteria,p__Bacteroidota,c__Bacteroidia,o__Bacteroidales,f__Bacteroidaceae,g__Prevotella,s__Prevotella copri'


def test_genome_gather_from_file_rank(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_two_files(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # make test2 results (identical to test1 except query_name and filename)
    g_res2 = runtmp.output("test2.gather.csv")
    test2_results = [x.replace("test1", "test2") for x in open(g_res, 'r')]
    with open(g_res2, 'w') as fp:
        for line in test2_results:
            fp.write(line)

    c.run_sourmash('tax', 'genome', '-g', g_res, g_res2, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out
    assert 'test2,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test2.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_two_files_empty_force(runtmp):
    # make test2 results (identical to test1 except query_name and filename)
    # add an empty file too, with --force -> should work
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    g_empty_csv = runtmp.output('g_empty.csv')
    with open(g_empty_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_empty_csv)

    g_res2 = runtmp.output("test2.gather.csv")
    test2_results = [x.replace("test1", "test2") for x in open(g_res, 'r')]
    with open(g_res2, 'w') as fp:
        for line in test2_results:
            fp.write(line)

    c.run_sourmash('tax', 'genome', '-g', g_res, g_res2, '-g', g_empty_csv,
                   '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0',
                   '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out
    assert 'test2,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test2.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_duplicate_filename(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '--gather-csv', g_res, '-g', g_res, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}'
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_from_file_duplicate_filename(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}'
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_from_file_duplicate_query(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # different filename, contents identical to test1
    g_res2 = runtmp.output("test2.gather.csv")
    with open(g_res2, 'w') as fp:
        for line in open(g_res, 'r'):
            fp.write(line)

    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res2}\n")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')
    assert c.last_result.status == -1
    print(str(exc.value))
    assert "Gather query test1 was found in more than one CSV. Cannot load from " in str(exc.value)


def test_genome_gather_from_file_duplicate_query_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # different filename, contents identical to test1
    g_res2 = runtmp.output("test2.gather.csv")
    with open(g_res2, 'w') as fp:
        for line in open(g_res, 'r'):
            fp.write(line)

    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res2}\n")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1

    assert "Gather query test1 was found in more than one CSV." in c.last_result.err
    assert "Cannot force past duplicated gather query. Exiting." in c.last_result.err


def test_genome_gather_cli_and_from_file(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")

    # make test2 results (identical to test1 except query_name)
    g_res2 = runtmp.output("test2.gather.csv")
    test2_results = [x.replace("test1", "test2") for x in open(g_res, 'r')]
    with open(g_res2, 'w') as fp:
        for line in test2_results:
            fp.write(line)

    # write test2 csv to a text file for input
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res2}\n")

    c.run_sourmash('tax', 'genome', '-g', g_res, '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out
    assert 'test2,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test2.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_cli_and_from_file_duplicate_filename(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")

    # also write test1 csv to a text file for input
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'genome', '-g', g_res, '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'ignoring duplicated reference to file: {g_res}' in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_gather_from_file_below_threshold(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--containment-threshold', '1')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,status,rank,fraction,lineage" in c.last_result.out
    assert "test1,below_threshold,,0.000," in c.last_result.out


def test_genome_gather_two_queries(runtmp):
    '''
    This checks for initial bug where classification
    would only happen for one genome per rank when
    doing --containment-threshold classification
    '''
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/47+63_x_gtdb-rs202.gather.csv')

    # split 47+63 into two fake queries: q47, q63
    g_res2 = runtmp.output("two-queries.gather.csv")
    q2_results = [x for x in open(g_res, 'r')]
    # rename queries
    q2_results[1] = q2_results[1].replace('47+63', 'q47')
    q2_results[2] = q2_results[2].replace('47+63', 'q63')
    with open(g_res2, 'w') as fp:
        for line in q2_results:
            print(line)
            fp.write(line)

    c.run_sourmash('tax', 'genome', '-g', g_res2, '--taxonomy-csv', taxonomy_csv,
                   '--containment-threshold', '0')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,status,rank,fraction,lineage" in c.last_result.out
    assert "q63,match,species,0.336,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Shewanellaceae;g__Shewanella;s__Shewanella baltica,491c0a81," in c.last_result.out
    assert "q47,match,species,0.664,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Shewanellaceae;g__Shewanella;s__Shewanella baltica," in c.last_result.out


def test_genome_rank_duplicated_taxonomy_fail(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1] + 'FOO') # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', duplicated_csv,
                       '--rank', 'species')
    assert "cannot read taxonomy assignments" in str(exc.value)
    assert "multiple lineages for identifier GCF_001881345" in str(exc.value)


def test_genome_rank_duplicated_taxonomy_fail_lineages(runtmp):
    # write temp taxonomy with duplicates => lineages-style file
    c = runtmp

    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    taxdb = tax_utils.LineageDB.load(taxonomy_csv)

    for k, v in taxdb.items():
        print(k, v)

    lineage_csv = runtmp.output('lin.csv')
    with open(lineage_csv, 'w', newline="") as fp:
        w = csv.writer(fp)
        w.writerow(['name', 'lineage'])
        for k, v in taxdb.items():
            linstr = lca_utils.display_lineage(v)
            w.writerow([k, linstr])

            # duplicate each row, changing something (truncate species, here)
            v = v[:-1]
            linstr = lca_utils.display_lineage(v)
            w.writerow([k, linstr])

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'summarize', lineage_csv)
        print(c.last_result.out)
        print(c.last_result.err)

    assert "cannot read taxonomy assignments" in str(exc.value)
    assert "multiple lineages for identifier GCF_001881345" in str(exc.value)


def test_genome_rank_duplicated_taxonomy_force(runtmp):
    # write temp taxonomy with duplicates
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', duplicated_csv,
                   '--rank', 'species', '--force', '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_missing_taxonomy_ignore_threshold(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, '--containment-threshold', '0')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_missing_taxonomy_recover_with_second_tax_file(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, '-t', taxonomy_csv, '--containment-threshold', '0')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" not in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_missing_taxonomy_ignore_rank(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, '--rank', 'species')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,below_threshold,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_multiple_taxonomy_files(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    # using mult -t args
    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, '-t', taxonomy_csv)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" not in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,' in c.last_result.out
    # using single -t arg
    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, taxonomy_csv)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" not in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,' in c.last_result.out


def test_genome_multiple_taxonomy_files_empty_force(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry, as well as an empty file,
    # and use force
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    empty_tax = runtmp.output('tax_empty.txt')
    with open(empty_tax, "w") as fp:
        fp.write("")
    
    # using mult -t args
    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv, '-t', taxonomy_csv, '-t', empty_tax, '--force')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" not in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,' in c.last_result.out


def test_genome_missing_taxonomy_fail_threshold(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv,
                       '--fail-on-missing-taxonomy', '--containment-threshold', '0')

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert "The following are missing from the taxonomy information: GCF_001881345" in str(exc.value)
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in str(exc.value)
    assert c.last_result.status == -1


def test_genome_missing_taxonomy_fail_rank(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', subset_csv,
                       '--fail-on-missing-taxonomy', '--rank', 'species')

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert "The following are missing from the taxonomy information: GCF_001881345" in str(exc.value)
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in str(exc.value)
    assert c.last_result.status == -1


def test_genome_rank_not_available(runtmp):
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--rank', 'strain', '--containment-threshold', '0')

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert "No taxonomic information provided for rank strain: cannot classify at this rank" in str(exc.value)


def test_genome_empty_gather_results_with_header_single(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results = [x for x in open(g_csv, 'r')]
    empty_gather_with_header = runtmp.output('g_header.csv')
    # write temp empty gather results (header only)
    with open(empty_gather_with_header, "w") as fp:
        fp.write(gather_results[0])

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', empty_gather_with_header, '--taxonomy-csv', taxonomy_csv)

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f'No gather results loaded from {empty_gather_with_header}.' in str(exc.value)
    assert 'Exiting.' in str(exc.value)


def test_genome_empty_gather_results_single(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results
    empty_tax = runtmp.output('tax_header.csv')
    with open(empty_tax, "w") as fp:
        fp.write("")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', empty_tax, '--taxonomy-csv', taxonomy_csv)


    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f"Cannot read gather results from '{empty_tax}'. Is file empty?" in str(exc.value)
    assert 'Exiting.' in c.last_result.err


def test_genome_empty_gather_results_single_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results (header only)
    empty_tax = runtmp.output('tax_header.csv')
    with open(empty_tax, "w") as fp:
        fp.write("")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', empty_tax, '--taxonomy-csv', taxonomy_csv,
                       '--force')

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert '--force is set. Attempting to continue to next set of gather results.' in str(exc.value)
    assert 'No results for classification. Exiting.' in str(exc.value)


def test_genome_empty_gather_results_with_empty_csv_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results
    empty_tax = runtmp.output('tax_empty.txt')
    with open(empty_tax, "w") as fp:
        fp.write("")

    g_from_file = runtmp.output("tmp-from-csv.csv")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{empty_tax}\n")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', empty_tax, '--from-file', g_from_file,
                       '--taxonomy-csv', taxonomy_csv, '--rank', 'species', '--force')

    print(str(exc.value))
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert '--force is set. Attempting to continue to next set of gather results.' in str(exc.value)
    assert 'No results for classification. Exiting.' in str(exc.value)


def test_genome_empty_gather_results_with_csv_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    # write temp empty gather results
    empty_tax = runtmp.output('tax_empty.csv')
    with open(empty_tax, "w") as fp:
        fp.write("")

    c.run_sourmash('tax', 'genome', '-g', empty_tax, '--from-file', g_from_file,
                   '--taxonomy-csv', taxonomy_csv, '--rank', 'species',
                   '--containment-threshold', '0', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert '--force is set. Attempting to continue to next set of gather results.' in c.last_result.err
    assert 'loaded 4 results total from 1 gather CSVs' in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out


def test_genome_containment_threshold_bounds(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    below_threshold = "-1"

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', tax, '--taxonomy-csv', tax,
                       '--containment-threshold', below_threshold)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Argument must be >0 and <1" in str(exc.value)

    above_threshold = "1.1"
    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--containment-threshold', above_threshold)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Argument must be >0 and <1" in str(exc.value)


def test_genome_containment_threshold_type(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    not_a_float = "str"

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--containment-threshold', not_a_float)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Must be a floating point number" in str(exc.value)


def test_genome_over100percent_error(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    perfect_g_csv = runtmp.output('g.csv')

    #create an impossible gather result
    with open(g_csv, 'r') as fp:
        r = csv.DictReader(fp, delimiter=',')
        header = r.fieldnames
        print(header)
        with open(perfect_g_csv, 'w') as out_fp:
            w = csv.DictWriter(out_fp, header)
            w.writeheader()
            for n, row in enumerate(r):
                if n == 0:
                    row["f_unique_to_query"] = 1.1
                w.writerow(row)
                print(row)

    with pytest.raises(SourmashCommandFailed):
        runtmp.run_sourmash('tax', 'genome', '-g', perfect_g_csv, '--taxonomy-csv', tax)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == -1
    assert "ERROR: The tax summary of query 'test1' is 1.1, which is > 100% of the query!!" in runtmp.last_result.err


def test_genome_ani_threshold_input_errors(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather_ani.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    below_threshold = "-1"

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', tax, '--taxonomy-csv', tax,
                       '--ani-threshold', below_threshold)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Argument must be >0 and <1" in str(exc.value)

    above_threshold = "1.1"
    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', above_threshold)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Argument must be >0 and <1" in str(exc.value)

    not_a_float = "str"

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', not_a_float)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "ERROR: Must be a floating point number" in str(exc.value)


def test_genome_ani_threshold(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather_ani.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', "0.95")

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "WARNING: Please run gather with sourmash >= 4.4 to estimate query ANI at rank. Continuing without ANI..." not in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,0.9328896594471843' in c.last_result.out 

    # more lax threshold
    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', "0.9")

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0'  in c.last_result.out

    # too stringent of threshold (using rank)
    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', "1.0", '--rank', 'species')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "WARNING: classifying query test1 at desired rank species does not meet query ANI/AAI threshold 1.0" in c.last_result.err
    assert "test1,below_threshold,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0,0.9247805047263588" in c.last_result.out


def test_genome_ani_oldgather(runtmp):
    # Ignore ANI if we don't have the information we need to estimate it
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,' in c.last_result.out

    c.run_sourmash('tax', 'genome', '-g', g_csv, '--taxonomy-csv', tax,
                       '--ani-threshold', "0.95")

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "WARNING: Please run gather with sourmash >= 4.4 to estimate query ANI at rank. Continuing without ANI..." in c.last_result.err
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,family,0.116,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae,md5,test1.sig,0.073,582000.0,' in c.last_result.out


def test_genome_ani_lemonade_classify(runtmp):
    # test a complete MAG classification with lemonade MAG from STAMPS 2022
    # (real data!)
    c = runtmp

    ## first run gather
    genome = utils.get_test_data('tax/lemonade-MAG3.sig.gz')
    matches = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.zip')

    c.run_sourmash('gather', genome, matches,
                   '--threshold-bp=5000', '-o', 'gather.csv')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    this_gather_file = c.output('gather.csv')
    this_gather = open(this_gather_file).readlines()

    assert len(this_gather) == 4

    ## now run 'tax genome' with human output
    taxonomy_file = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.tax.csv')
    c.run_sourmash('tax', 'genome', '-g', this_gather_file, '-t', taxonomy_file,
                   '--ani', '0.8', '-F', 'human')

    output = c.last_result.out
    assert 'MAG3_1             5.3%     91.0%  d__Bacteria;p__Bacteroidota;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Prosthecochloris;s__Prosthecochloris vibrioformis' in output

    # aaand classify to lineage_csv
    c.run_sourmash('tax', 'genome', '-g', this_gather_file, '-t', taxonomy_file,
                   '--ani', '0.8', '-F', 'lineage_csv')

    output = c.last_result.out
    assert 'ident,superkingdom,phylum,class,order,family,genus,species' in output
    assert 'MAG3_1,d__Bacteria,p__Bacteroidota,c__Chlorobia,o__Chlorobiales,f__Chlorobiaceae,g__Prosthecochloris,s__Prosthecochloris vibrioformis' in output


def test_metagenome_no_gather_csv(runtmp):
    # test tax metagenome with no -g
    taxonomy_file = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.tax.csv')
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'metagenome', '-t', taxonomy_file)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)


def test_genome_no_gather_csv(runtmp):
    # test tax genome with no -g
    taxonomy_file = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.tax.csv')
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'genome', '-t', taxonomy_file)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)


def test_annotate_no_gather_csv(runtmp):
    # test tax annotate with no -g
    taxonomy_file = utils.get_test_data('tax/lemonade-MAG3.x.gtdb.matches.tax.csv')
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'annotate', '-t', taxonomy_file)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)


def test_annotate_0(runtmp):
    # test annotate basics
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    c.run_sourmash('tax', 'annotate', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', out_dir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))
    assert f"saving 'annotate' output to '{csvout}'" in runtmp.last_result.err

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[2]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in lin_gather_results[3]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[4]


def test_annotate_gzipped_gather(runtmp):
    # test annotate basics
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    # rewrite gather_csv as gzipped csv
    gz_gather = runtmp.output('test1.gather.csv.gz')
    with open(g_csv, 'rb') as f_in, gzip.open(gz_gather, 'wb') as f_out:
        f_out.writelines(f_in)

    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    c.run_sourmash('tax', 'annotate', '--gather-csv', gz_gather, '--taxonomy-csv', tax, '-o', out_dir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))
    assert f"saving 'annotate' output to '{csvout}'" in runtmp.last_result.err

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[2]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in lin_gather_results[3]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[4]


def test_annotate_gather_argparse(runtmp):
    # test annotate with two gather CSVs, second one empty, and --force.
    # this tests argparse handling w/extend.
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    g_empty_csv = runtmp.output('g_empty.csv')
    with open(g_empty_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_empty_csv)

    c.run_sourmash('tax', 'annotate', '--gather-csv', g_csv,
                   '-g', g_empty_csv, '--taxonomy-csv', tax, '-o', out_dir,
                   '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))
    assert f"saving 'annotate' output to '{csvout}'" in runtmp.last_result.err

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]


def test_annotate_0_db(runtmp):
    # test annotate with sqlite db
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.db')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    c.run_sourmash('tax', 'annotate', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', out_dir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))
    assert f"saving 'annotate' output to '{csvout}'" in runtmp.last_result.err

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[2]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in lin_gather_results[3]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[4]


def test_annotate_empty_gather_results(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    #creates empty gather result
    g_csv = runtmp.output('g.csv')
    with open(g_csv, "w") as fp:
        fp.write("")
    print("g_csv: ", g_csv)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'annotate', '-g', g_csv, '--taxonomy-csv', tax)

    assert f"Cannot read gather results from '{g_csv}'. Is file empty?" in str(exc.value)
    assert runtmp.last_result.status == -1


def test_annotate_bad_gather_header(runtmp):
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    bad_g_csv = runtmp.output('g.csv')

    #creates bad gather result
    bad_g = [x.replace("query_name", "nope") for x in open(g_csv, 'r')]
    with open(bad_g_csv, 'w') as fp:
        for line in bad_g:
            fp.write(line)
    print("bad_gather_results: \n", bad_g)

    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'annotate', '-g', bad_g_csv, '--taxonomy-csv', tax)

    assert f"Not all required gather columns are present in '{bad_g_csv}'." in str(exc.value)
    assert runtmp.last_result.status == -1


def test_annotate_empty_tax_lineage_input(runtmp):
    tax_empty = runtmp.output('t.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)


    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'annotate', '-g', g_csv, '--taxonomy-csv', tax_empty)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status != 0
    assert "cannot read taxonomy assignments from" in str(exc.value)


def test_annotate_empty_tax_lineage_input_recover_with_second_taxfile(runtmp):
    tax_empty = runtmp.output('t.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    runtmp.run_sourmash('tax', 'annotate', '-g', g_csv, '-t', tax_empty, '--taxonomy-csv', tax, '--force')

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0


def test_annotate_empty_tax_lineage_input_recover_with_second_taxfile_2(runtmp):
    # test with empty tax second, to check on argparse handling
    tax_empty = runtmp.output('t.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    runtmp.run_sourmash('tax', 'annotate', '-g', g_csv,
                        '--taxonomy-csv', tax, '-t', tax_empty, '--force')

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0


def test_tax_prepare_1_csv_to_csv(runtmp, keep_identifiers, keep_versions):
    # CSV -> CSV; same assignments
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    taxout = runtmp.output('out.csv')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                                taxout, '-F', 'csv', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                        taxout, '-F', 'csv', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)

    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_1_combine_csv(runtmp):
    # multiple CSVs to a single combined CSV
    tax1 = utils.get_test_data('tax/test.taxonomy.csv')
    tax2 = utils.get_test_data('tax/protozoa_genbank_lineage.csv')

    taxout = runtmp.output('out.csv')

    runtmp.sourmash('tax', 'prepare', '-t', tax1, tax2, '-F', 'csv',
                    '-o', taxout)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert not out
    assert "...loaded 8 entries" in err

    out = open(taxout).readlines()
    assert len(out) == 9


def test_tax_prepare_1_csv_to_csv_empty_ranks(runtmp, keep_identifiers, keep_versions):
    # CSV -> CSV; same assignments, even when trailing ranks are empty
    tax = utils.get_test_data('tax/test-empty-ranks.taxonomy.csv')
    taxout = runtmp.output('out.csv')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                                taxout, '-F', 'csv', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                        taxout, '-F', 'csv', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)

    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_1_csv_to_csv_empty_file(runtmp, keep_identifiers, keep_versions):
    # CSV -> CSV with an empty input file and --force
    # tests argparse extend
    tax = utils.get_test_data('tax/test-empty-ranks.taxonomy.csv')
    tax_empty = runtmp.output('t.csv')
    taxout = runtmp.output('out.csv')

    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                                taxout, '-F', 'csv', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-t', tax_empty, '-o',
                        taxout, '-F', 'csv', *args, '--force')
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)

    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_1_csv_to_csv_empty_ranks_2(runtmp, keep_identifiers, keep_versions):
    # CSV -> CSV; same assignments for situations with empty internal ranks
    tax = utils.get_test_data('tax/test-empty-ranks-2.taxonomy.csv')
    taxout = runtmp.output('out.csv')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                                taxout, '-F', 'csv', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                        taxout, '-F', 'csv', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)

    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_1_csv_to_csv_empty_ranks_3(runtmp, keep_identifiers, keep_versions):
    # CSV -> CSV; same assignments for situations with empty internal ranks
    tax = utils.get_test_data('tax/test-empty-ranks-3.taxonomy.csv')
    taxout = runtmp.output('out.csv')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                                taxout, '-F', 'csv', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o',
                        taxout, '-F', 'csv', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)

    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_2_csv_to_sql(runtmp, keep_identifiers, keep_versions):
    # CSV -> SQL; same assignments?
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    taxout = runtmp.output('out.db')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                                '-F', 'sql', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                        '-F', 'sql', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)
    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)

    # cannot overwrite -
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                            '-F', 'sql', *args)
    assert 'taxonomy table already exists' in str(exc.value)


def test_tax_prepare_2_csv_to_sql_empty_ranks(runtmp, keep_identifiers, keep_versions):
    # CSV -> SQL with some empty ranks in the taxonomy file
    tax = utils.get_test_data('tax/test-empty-ranks.taxonomy.csv')
    taxout = runtmp.output('out.db')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                                '-F', 'sql', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                        '-F', 'sql', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)
    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_3_db_to_csv(runtmp):
    # SQL -> CSV; same assignments
    taxcsv = utils.get_test_data('tax/test.taxonomy.csv')
    taxdb = utils.get_test_data('tax/test.taxonomy.db')
    taxout = runtmp.output('out.csv')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxdb,
                        '-o', taxout, '-F', 'csv')
    assert os.path.exists(taxout)
    with open(taxout) as fp:
        print(fp.read())

    db1 = tax_utils.MultiLineageDB.load([taxcsv],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)

    db2 = tax_utils.MultiLineageDB.load([taxout])
    db3 = tax_utils.MultiLineageDB.load([taxdb],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)
    assert set(db1) == set(db2)
    assert set(db1) == set(db3)


def test_tax_prepare_3_db_to_csv_gz(runtmp):
    # SQL -> CSV; same assignments
    taxcsv = utils.get_test_data('tax/test.taxonomy.csv')
    taxdb = utils.get_test_data('tax/test.taxonomy.db')
    taxout = runtmp.output('out.csv.gz')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxdb,
                        '-o', taxout, '-F', 'csv')
    assert os.path.exists(taxout)
    with gzip.open(taxout, 'rt') as fp:
        print(fp.read())

    db1 = tax_utils.MultiLineageDB.load([taxcsv],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)

    db2 = tax_utils.MultiLineageDB.load([taxout])
    db3 = tax_utils.MultiLineageDB.load([taxdb],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)
    assert set(db1) == set(db2)
    assert set(db1) == set(db3)


def test_tax_prepare_2_csv_to_sql_empty_ranks_2(runtmp, keep_identifiers, keep_versions):
    # CSV -> SQL with some empty internal ranks in the taxonomy file
    tax = utils.get_test_data('tax/test-empty-ranks-2.taxonomy.csv')
    taxout = runtmp.output('out.db')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                                '-F', 'sql', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                        '-F', 'sql', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)
    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_2_csv_to_sql_empty_ranks_3(runtmp, keep_identifiers, keep_versions):
    # CSV -> SQL with some empty internal ranks in the taxonomy file
    tax = utils.get_test_data('tax/test-empty-ranks-3.taxonomy.csv')
    taxout = runtmp.output('out.db')

    args = []
    if keep_identifiers:
        args.append('--keep-full-identifiers')
    if keep_versions:
        args.append('--keep-identifier-versions')

    # this is an error - can't strip versions if not splitting identifiers
    if keep_identifiers and not keep_versions:
        with pytest.raises(SourmashCommandFailed):
            runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                                '-F', 'sql', *args)
        return

    runtmp.run_sourmash('tax', 'prepare', '-t', tax, '-o', taxout,
                        '-F', 'sql', *args)
    assert os.path.exists(taxout)

    db1 = tax_utils.MultiLineageDB.load([tax],
                                        keep_full_identifiers=keep_identifiers,
                                        keep_identifier_versions=keep_versions)
    db2 = tax_utils.MultiLineageDB.load([taxout])

    assert set(db1) == set(db2)


def test_tax_prepare_3_db_to_csv_empty_ranks(runtmp):
    # SQL -> CSV; same assignments, with empty ranks
    taxcsv = utils.get_test_data('tax/test-empty-ranks.taxonomy.csv')
    taxdb = utils.get_test_data('tax/test-empty-ranks.taxonomy.db')
    taxout = runtmp.output('out.csv')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxdb,
                        '-o', taxout, '-F', 'csv')
    assert os.path.exists(taxout)
    with open(taxout) as fp:
        print(fp.read())

    db1 = tax_utils.MultiLineageDB.load([taxcsv],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)

    db2 = tax_utils.MultiLineageDB.load([taxout])
    db3 = tax_utils.MultiLineageDB.load([taxdb],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)
    assert set(db1) == set(db2)
    assert set(db1) == set(db3)


def test_tax_prepare_3_db_to_csv_empty_ranks_2(runtmp):
    # SQL -> CSV; same assignments, with empty ranks
    taxcsv = utils.get_test_data('tax/test-empty-ranks-2.taxonomy.csv')
    taxdb = utils.get_test_data('tax/test-empty-ranks-2.taxonomy.db')
    taxout = runtmp.output('out.csv')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxdb,
                        '-o', taxout, '-F', 'csv')
    assert os.path.exists(taxout)
    with open(taxout) as fp:
        print(fp.read())

    db1 = tax_utils.MultiLineageDB.load([taxcsv],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)

    db2 = tax_utils.MultiLineageDB.load([taxout])
    db3 = tax_utils.MultiLineageDB.load([taxdb],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)
    assert set(db1) == set(db2)
    assert set(db1) == set(db3)


def test_tax_prepare_3_db_to_csv_empty_ranks_3(runtmp):
    # SQL -> CSV; same assignments, with empty ranks
    taxcsv = utils.get_test_data('tax/test-empty-ranks-3.taxonomy.csv')
    taxdb = utils.get_test_data('tax/test-empty-ranks-3.taxonomy.db')
    taxout = runtmp.output('out.csv')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxdb,
                        '-o', taxout, '-F', 'csv')
    assert os.path.exists(taxout)
    with open(taxout) as fp:
        print(fp.read())

    db1 = tax_utils.MultiLineageDB.load([taxcsv],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)

    db2 = tax_utils.MultiLineageDB.load([taxout])
    db3 = tax_utils.MultiLineageDB.load([taxdb],
                                        keep_full_identifiers=False,
                                        keep_identifier_versions=False)
    assert set(db1) == set(db2)
    assert set(db1) == set(db3)


def test_tax_prepare_sqlite_lineage_version(runtmp):
    # test bad sourmash_internals version for SqliteLineage
    taxcsv = utils.get_test_data('tax/test.taxonomy.csv')
    taxout = runtmp.output('out.db')

    runtmp.run_sourmash('tax', 'prepare', '-t', taxcsv,
                        '-o', taxout, '-F', 'sql')
    assert os.path.exists(taxout)

    # set bad version
    conn = sqlite_utils.open_sqlite_db(taxout)
    c = conn.cursor()
    c.execute("UPDATE sourmash_internal SET value='0.9' WHERE key='SqliteLineage'")

    conn.commit()
    conn.close()

    with pytest.raises(IndexNotSupported):
        db = tax_utils.MultiLineageDB.load([taxout])


def test_tax_prepare_sqlite_no_lineage():
    # no lineage table at all
    sqldb = utils.get_test_data('sqlite/index.sqldb')

    with pytest.raises(ValueError):
        db = tax_utils.MultiLineageDB.load([sqldb])


def test_tax_grep_exists(runtmp):
    # test that 'tax grep' exists

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('tax', 'grep')

    err = runtmp.last_result.err
    assert 'usage:' in err


def test_tax_grep_search_shew(runtmp):
    # test 'tax grep Shew'
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', 'Shew', '-t', taxfile)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert lines[1][0] == 'GCF_000017325.1'
    assert lines[2][0] == 'GCF_000021665.1'
    assert len(lines) == 3

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert 'found 2 matches; saved identifiers to picklist' in err


def test_tax_grep_search_shew_out(runtmp):
    # test 'tax grep Shew', save result to a file
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', 'Shew', '-t', taxfile, '-o', 'pick.csv')

    err = runtmp.last_result.err

    out = open(runtmp.output('pick.csv')).read()
    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert lines[1][0] == 'GCF_000017325.1'
    assert lines[2][0] == 'GCF_000021665.1'
    assert len(lines) == 3

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert 'found 2 matches; saved identifiers to picklist' in err


def test_tax_grep_search_shew_sqldb_out(runtmp):
    # test 'tax grep Shew' on a sqldb, save result to a file
    taxfile = utils.get_test_data('tax/test.taxonomy.db')

    runtmp.sourmash('tax', 'grep', 'Shew', '-t', taxfile, '-o', 'pick.csv')

    err = runtmp.last_result.err

    out = open(runtmp.output('pick.csv')).read()
    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert lines[1][0] == 'GCF_000017325'
    assert lines[2][0] == 'GCF_000021665'
    assert len(lines) == 3

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert 'found 2 matches; saved identifiers to picklist' in err


def test_tax_grep_search_shew_lowercase(runtmp):
    # test 'tax grep shew' (lowercase), save result to a file
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', 'shew', '-t', taxfile, '-o', 'pick.csv')

    err = runtmp.last_result.err
    assert "searching 1 taxonomy files for 'shew'" in err
    assert 'found 0 matches; saved identifiers to picklist' in err

    runtmp.sourmash('tax', 'grep', '-i', 'shew',
                    '-t', taxfile, '-o', 'pick.csv')

    err = runtmp.last_result.err
    assert "searching 1 taxonomy files for 'shew'" in err
    assert 'found 2 matches; saved identifiers to picklist' in err

    out = open(runtmp.output('pick.csv')).read()
    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert lines[1][0] == 'GCF_000017325.1'
    assert lines[2][0] == 'GCF_000021665.1'
    assert len(lines) == 3


def test_tax_grep_search_shew_out_use_picklist(runtmp):
    # test 'tax grep Shew', output to a picklist, use picklist
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')
    dbfile = utils.get_test_data('tax/gtdb-tax-grep.sigs.zip')

    runtmp.sourmash('tax', 'grep', 'Shew', '-t', taxfile, '-o', 'pick.csv')

    runtmp.sourmash('sig', 'cat', dbfile, '--picklist',
                    'pick.csv:ident:ident', '-o', 'pick-out.zip')

    all_sigs = sourmash.load_file_as_index(dbfile)
    assert len(all_sigs) == 3

    pick_sigs = sourmash.load_file_as_index(runtmp.output('pick-out.zip'))
    assert len(pick_sigs) == 2

    names = [ ss.name.split()[0] for ss in pick_sigs.signatures() ]
    assert len(names) == 2
    assert 'GCF_000017325.1' in names
    assert 'GCF_000021665.1' in names


def test_tax_grep_search_shew_invert(runtmp):
    # test 'tax grep -v Shew'
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', '-v', 'Shew', '-t', taxfile)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "-v/--invert-match specified; returning only lineages that do not match." in err

    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert lines[1][0] == 'GCF_001881345.1'
    assert lines[2][0] == 'GCF_003471795.1'
    assert len(lines) == 5

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert 'found 4 matches; saved identifiers to picklist' in err

    all_names = set([ x[0] for x in lines ])
    assert 'GCF_000017325.1' not in all_names
    assert 'GCF_000021665.1' not in all_names


def test_tax_grep_search_shew_invert_select_phylum(runtmp):
    # test 'tax grep -v Shew -r phylum'
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', '-v', 'Shew', '-t', taxfile, '-r', 'phylum')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "-v/--invert-match specified; returning only lineages that do not match." in err
    assert "limiting matches to phylum"

    lines = [ x.strip() for x in out.splitlines() ]
    lines = [ x.split(',') for x in lines ]
    assert lines[0][0] == 'ident'
    assert len(lines) == 7

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert 'found 6 matches; saved identifiers to picklist' in err

    all_names = set([ x[0] for x in lines ])
    assert 'GCF_000017325.1' in all_names
    assert 'GCF_000021665.1' in all_names


def test_tax_grep_search_shew_invert_select_bad_rank(runtmp):
    # test 'tax grep -v Shew -r badrank' - should fail
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('tax', 'grep', '-v', 'Shew', '-t', taxfile,
                        '-r', 'badrank')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(err)
    assert 'error: argument -r/--rank: invalid choice:' in err


def test_tax_grep_search_shew_count(runtmp):
    # test 'tax grep Shew --count'
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'grep', 'Shew', '-t', taxfile, '-c')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert not out.strip()

    assert "searching 1 taxonomy files for 'Shew'" in err
    assert not 'found 2 matches; saved identifiers to picklist' in err


def test_tax_grep_multiple_csv(runtmp):
    # grep on multiple CSVs
    tax1 = utils.get_test_data('tax/test.taxonomy.csv')
    tax2 = utils.get_test_data('tax/protozoa_genbank_lineage.csv')

    taxout = runtmp.output('out.csv')

    runtmp.sourmash('tax', 'grep', "Toxo|Gamma",
                    '-t', tax1, tax2,
                    '-o', taxout)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert not out
    assert "found 4 matches" in err

    lines = open(taxout).readlines()
    assert len(lines) == 5

    names = set([ x.split(',')[0] for x in lines ])
    assert 'GCA_000256725' in names
    assert 'GCF_000017325.1' in names
    assert 'GCF_000021665.1' in names
    assert 'GCF_001881345.1' in names


def test_tax_grep_multiple_csv_empty_force(runtmp):
    # grep on multiple CSVs, one empty, with --force
    tax1 = utils.get_test_data('tax/test.taxonomy.csv')
    tax2 = utils.get_test_data('tax/protozoa_genbank_lineage.csv')
    tax_empty = runtmp.output('t.csv')

    taxout = runtmp.output('out.csv')
    with open(tax_empty, "w") as fp:
        fp.write("")
    print("t_csv: ", tax_empty)

    runtmp.sourmash('tax', 'grep', "Toxo|Gamma",
                    '-t', tax1, tax2, '-t', tax_empty,
                    '-o', taxout, '--force')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert not out
    assert "found 4 matches" in err

    lines = open(taxout).readlines()
    assert len(lines) == 5

    names = set([ x.split(',')[0] for x in lines ])
    assert 'GCA_000256725' in names
    assert 'GCF_000017325.1' in names
    assert 'GCF_000021665.1' in names
    assert 'GCF_001881345.1' in names


def test_tax_grep_duplicate_csv(runtmp):
    # grep on duplicates => should collapse to uniques on identifiers
    tax1 = utils.get_test_data('tax/test.taxonomy.csv')

    taxout = runtmp.output('out.csv')

    runtmp.sourmash('tax', 'grep', "Gamma",
                    '-t', tax1, tax1,
                    '-o', taxout)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert not out
    assert "found 3 matches" in err

    lines = open(taxout).readlines()
    assert len(lines) == 4

    names = set([ x.split(',')[0] for x in lines ])
    assert 'GCF_000017325.1' in names
    assert 'GCF_000021665.1' in names
    assert 'GCF_001881345.1' in names


def test_tax_summarize(runtmp):
    # test basic operation with summarize
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'summarize', taxfile)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "rank superkingdom:        1 distinct taxonomic lineages" in out
    assert "rank phylum:              2 distinct taxonomic lineages" in out
    assert "rank class:               2 distinct taxonomic lineages" in out
    assert "rank order:               2 distinct taxonomic lineages" in out
    assert "rank family:              3 distinct taxonomic lineages" in out
    assert "rank genus:               4 distinct taxonomic lineages" in out
    assert "rank species:             4 distinct taxonomic lineages" in out


def test_tax_summarize_multiple(runtmp):
    # test basic operation with summarize on multiple files
    tax1 = utils.get_test_data('tax/bacteria_refseq_lineage.csv')
    tax2 = utils.get_test_data('tax/protozoa_genbank_lineage.csv')

    runtmp.sourmash('tax', 'summarize', tax1, tax2)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "rank superkingdom:        2 distinct taxonomic lineages" in out
    assert "rank phylum:              3 distinct taxonomic lineages" in out
    assert "rank class:               4 distinct taxonomic lineages" in out
    assert "rank order:               4 distinct taxonomic lineages" in out
    assert "rank family:              5 distinct taxonomic lineages" in out
    assert "rank genus:               5 distinct taxonomic lineages" in out
    assert "rank species:             5 distinct taxonomic lineages" in out


def test_tax_summarize_empty_line(runtmp):
    # test basic operation with summarize on a file w/empty line
    taxfile = utils.get_test_data('tax/test-empty-line.taxonomy.csv')

    runtmp.sourmash('tax', 'summarize', taxfile)

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "rank superkingdom:        1 distinct taxonomic lineages" in out
    assert "rank phylum:              2 distinct taxonomic lineages" in out
    assert "rank class:               2 distinct taxonomic lineages" in out
    assert "rank order:               2 distinct taxonomic lineages" in out
    assert "rank family:              3 distinct taxonomic lineages" in out
    assert "rank genus:               4 distinct taxonomic lineages" in out
    assert "rank species:             4 distinct taxonomic lineages" in out


def test_tax_summarize_empty(runtmp):
    # test failure on empty file
    taxfile = runtmp.output('no-exist')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('tax', 'summarize', taxfile)

    out = runtmp.last_result.out
    err = runtmp.last_result.err
    assert "ERROR while loading taxonomies" in err


def test_tax_summarize_csv(runtmp):
    # test basic operation w/csv output
    taxfile = utils.get_test_data('tax/test.taxonomy.csv')

    runtmp.sourmash('tax', 'summarize', taxfile, '-o', 'ranks.csv')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "saved 18 lineage counts to 'ranks.csv'" in err

    csv_out = runtmp.output('ranks.csv')

    with sourmash_args.FileInputCSV(csv_out) as r:
        # count number across ranks as a cheap consistency check
        c = Counter()
        for row in r:
            val = row['lineage_count']
            c[val] += 1

        assert c['3'] == 7
        assert c['2'] == 5
        assert c['1'] == 5


def test_tax_summarize_on_annotate(runtmp):
    # test summarize on output of annotate basics
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    runtmp.run_sourmash('tax', 'annotate', '--gather-csv', g_csv, '--taxonomy-csv', tax, '-o', out_dir)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    # so far so good - now see if we can run summarize!

    runtmp.run_sourmash('tax', 'summarize', csvout)
    out = runtmp.last_result.out
    err = runtmp.last_result.err

    print(out)
    print(err)

    assert "number of distinct taxonomic lineages: 4" in out
    assert "rank superkingdom:        1 distinct taxonomic lineages" in out
    assert "rank phylum:              2 distinct taxonomic lineages" in out
    assert "rank class:               2 distinct taxonomic lineages" in out
    assert "rank order:               2 distinct taxonomic lineages" in out
    assert "rank family:              2 distinct taxonomic lineages" in out
    assert "rank genus:               3 distinct taxonomic lineages" in out
    assert "rank species:             3 distinct taxonomic lineages" in out


def test_tax_summarize_strain_csv(runtmp):
    # test basic operation w/csv output on taxonomy with strains
    taxfile = utils.get_test_data('tax/test-strain.taxonomy.csv')

    runtmp.sourmash('tax', 'summarize', taxfile, '-o', 'ranks.csv')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "saved 24 lineage counts to 'ranks.csv'" in err

    csv_out = runtmp.output('ranks.csv')

    with sourmash_args.FileInputCSV(csv_out) as r:
        # count number across ranks as a cheap consistency check
        c = Counter()
        for row in r:
            print(row)
            val = row['lineage_count']
            c[val] += 1

        print(list(c.most_common()))

        assert c['3'] == 7
        assert c['2'] == 5
        assert c['6'] == 1
        assert c['1'] == 11


def test_tax_summarize_strain_csv_with_lineages(runtmp):
    # test basic operation w/csv output on lineages-style file w/strain csv
    taxfile = utils.get_test_data('tax/test-strain.taxonomy.csv')
    lineage_csv = runtmp.output('lin-with-strains.csv')

    taxdb = tax_utils.LineageDB.load(taxfile)
    with open(lineage_csv, 'w', newline="") as fp:
        w = csv.writer(fp)
        w.writerow(['name', 'lineage'])
        for k, v in taxdb.items():
            linstr = lca_utils.display_lineage(v)
            w.writerow([k, linstr])

    runtmp.sourmash('tax', 'summarize', lineage_csv, '-o', 'ranks.csv')

    out = runtmp.last_result.out
    err = runtmp.last_result.err

    assert "number of distinct taxonomic lineages: 6" in out
    assert "saved 24 lineage counts to" in err

    csv_out = runtmp.output('ranks.csv')

    with sourmash_args.FileInputCSV(csv_out) as r:
        # count number across ranks as a cheap consistency check
        c = Counter()
        for row in r:
            print(row)
            val = row['lineage_count']
            c[val] += 1

        print(list(c.most_common()))

        assert c['3'] == 7
        assert c['2'] == 5
        assert c['6'] == 1
        assert c['1'] == 11
