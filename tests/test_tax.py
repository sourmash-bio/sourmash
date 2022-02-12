"""
Tests for the 'sourmash tax' command line and high level API.
"""
import os
import csv
import pytest

import sourmash_tst_utils as utils
from sourmash.tax import tax_utils
from sourmash_tst_utils import SourmashCommandFailed

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
    assert f"saving `csv_summary` output to {csvout}" in runtmp.last_result.err
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
    assert f"saving `krona` output to {csvout}" in runtmp.last_result.err

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
    assert f"saving `lineage_summary` output to {csvout}" in runtmp.last_result.err

    gn_lineage_summary = [x.rstrip().split('\t') for x in open(csvout)]
    print("species lineage summary results: \n", gn_lineage_summary)
    assert ['lineage', 'test1'] == gn_lineage_summary[0]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola', '0.027522935779816515'] == gn_lineage_summary[1]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella', '0.0885520542481053'] == gn_lineage_summary[2]
    assert ['d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia', '0.08815317112086159'] == gn_lineage_summary[3]
    assert ['unclassified', '0.7957718388512166']  == gn_lineage_summary[4]


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

    assert "of 6, missed 2 lineage assignments." in c.last_result.err
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

    assert "of 6, missed 0 lineage assignments." in c.last_result.err
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

    assert f'Cannot read gather results from {g_csv}. Is file empty?' in str(exc.value)
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

    assert f'Not all required gather columns are present in {bad_g_csv}.' in str(exc.value)
    assert runtmp.last_result.status == -1


def test_metagenome_empty_tax_lineage_input(runtmp):
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
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # different filename, contents identical to test1
    g_res2 = runtmp.output("test2.gather.csv")
    with open(g_res2, 'w') as fp:
        for line in open(g_res, 'r'):
            fp.write(line)

    c.run_sourmash('tax', 'metagenome',  '--gather-csv', g_res, g_res2,
                   '--taxonomy-csv', taxonomy_csv, '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert '--force is set, ignoring duplicate query.' in c.last_result.err
    assert 'No gather results loaded from ' in c.last_result.err
    assert 'loaded results from 1 gather CSVs' in c.last_result.err
    assert 'query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,superkingdom,0.204,d__Bacteria,md5,test1.sig,0.131,1024000' in c.last_result.out
    assert 'test1,superkingdom,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out
    assert 'test1,phylum,0.116,d__Bacteria;p__Bacteroidota,md5,test1.sig,0.073,582000' in c.last_result.out
    assert 'test1,phylum,0.088,d__Bacteria;p__Proteobacteria,md5,test1.sig,0.058,442000' in c.last_result.out
    assert 'test1,phylum,0.796,unclassified,md5,test1.sig,0.869,3990000' in c.last_result.out


def test_metagenome_gather_duplicate_filename(runtmp):
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

    assert f'Cannot read gather results from {g_csv}. Is file empty?' in str(exc.value)
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

    assert f'Not all required gather columns are present in {bad_g_csv}.' in str(exc.value)
    assert runtmp.last_result.status == -1


def test_genome_empty_tax_lineage_input(runtmp):
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

    assert f"saving `classification` output to {csvout}" in runtmp.last_result.err
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

    assert f"saving `krona` output to {csvout}" in runtmp.last_result.err
    assert c.last_result.status == 0
    kr_results = [x.rstrip().split('\t') for x in open(csvout)]
    print(kr_results)
    assert ['fraction', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']  == kr_results[0]
    assert ['0.0885520542481053', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Prevotella', 's__Prevotella copri'] == kr_results[1]


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


def test_genome_gather_duplicate_filename(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'genome', '--gather-csv', g_res, g_res, '--taxonomy-csv', taxonomy_csv,
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

    c.run_sourmash('tax', 'genome', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species', '--containment-threshold', '0', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank' in c.last_result.out
    assert 'test1,match,species,0.089,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri,md5,test1.sig,0.057,444000.0' in c.last_result.out
    assert '--force is set, ignoring duplicate query.' in c.last_result.err
    assert 'No gather results loaded from ' in c.last_result.err
    assert 'loaded results from 1 gather CSVs' in c.last_result.err


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


def test_genome_rank_duplicated_taxonomy_force(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
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
    assert f'Cannot read gather results from {empty_tax}. Is file empty?' in str(exc.value)
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
    assert 'loaded results from 1 gather CSVs' in c.last_result.err
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


def test_annotate_0(runtmp):
    # test annotate
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

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))
    assert f"saving `annotate` output to {csvout}" in runtmp.last_result.err

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[2]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in lin_gather_results[3]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[4]


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
    assert f"saving `annotate` output to {csvout}" in runtmp.last_result.err

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

    assert f'Cannot read gather results from {g_csv}. Is file empty?' in str(exc.value)
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

    assert f'Not all required gather columns are present in {bad_g_csv}.' in str(exc.value)
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


def test_tax_prepare_3_db_to_csv(runtmp):
    # CSV -> CSV; same assignments
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
