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

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert 'test1,superkingdom,0.131,d__Bacteria' in c.last_result.out
    assert "test1,phylum,0.073,d__Bacteria;p__Bacteroidota" in c.last_result.out
    assert "test1,phylum,0.058,d__Bacteria;p__Proteobacteria" in c.last_result.out
    assert "test1,class,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in c.last_result.out
    assert "test1,class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in c.last_result.out
    assert "test1,order,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in c.last_result.out
    assert "test1,order,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales" in c.last_result.out
    assert "test1,family,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae" in c.last_result.out
    assert "test1,family,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae" in c.last_result.out
    assert "test1,genus,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia" in c.last_result.out
    assert "test1,genus,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella" in c.last_result.out
    assert "test1,genus,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out
    assert "test1,species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in c.last_result.out
    assert "test1,species,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in c.last_result.out


def test_summarize_summary_csv_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    sum_csv = csv_base + ".summarized.csv"
    csvout = runtmp.output(sum_csv)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '-o', csv_base)

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    sum_gather_results = [x.rstrip() for x in open(csvout)]
    assert "query_name,rank,fraction,lineage" in sum_gather_results[0]
    assert 'test1,superkingdom,0.131,d__Bacteria' in sum_gather_results[1]
    assert "test1,phylum,0.073,d__Bacteria;p__Bacteroidota" in sum_gather_results[2]
    assert "test1,phylum,0.058,d__Bacteria;p__Proteobacteria" in sum_gather_results[3]
    assert "test1,class,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in sum_gather_results[4]
    assert "test1,class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in sum_gather_results[5]
    assert "test1,order,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in sum_gather_results[6]
    assert "test1,order,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales" in sum_gather_results[7]
    assert "test1,family,0.073,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae" in sum_gather_results[8]
    assert "test1,family,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae" in sum_gather_results[9]
    assert "test1,genus,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia" in sum_gather_results[10]
    assert "test1,genus,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella" in sum_gather_results[11]
    assert "test1,genus,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola" in sum_gather_results[12]
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in sum_gather_results[13]
    assert "test1,species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in sum_gather_results[14]
    assert "test1,species,0.016,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in sum_gather_results[15]


def test_summarize_krona_tsv_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    kr_csv = csv_base + ".krona.tsv"
    csvout = runtmp.output(kr_csv)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'krona', '--rank', 'genus')

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    gn_krona_results = [x.rstrip().split('\t') for x in open(csvout)]
    print("species krona results: \n", gn_krona_results)
    assert ['fraction', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus'] == gn_krona_results[0]
    assert ['0.05815279361459521', 'd__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacterales', 'f__Enterobacteriaceae', 'g__Escherichia']  == gn_krona_results[1]
    assert ['0.05701254275940707', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Prevotella'] == gn_krona_results[2]
    assert ['0.015637726014008795', 'd__Bacteria', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Phocaeicola'] == gn_krona_results[3]


def test_summarize_lineage_summary_out(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    lin_csv = csv_base + ".lineage_summary.tsv"
    csvout = runtmp.output(lin_csv)
    print("csvout: ", csvout)

    runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'lineage_summary', '--rank', 'genus')

    print(runtmp.last_result.status)
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0
    assert os.path.exists(csvout)

    gn_lineage_summary = [x.rstrip().split('\t') for x in open(csvout)]
    print("species lineage summary results: \n", gn_lineage_summary)
    assert ['lineage', 'test1'] == gn_lineage_summary[0]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola', '0.015637726014008795'] == gn_lineage_summary[1]
    assert ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella', '0.05701254275940707'] == gn_lineage_summary[2]
    assert ['d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia', '0.05815279361459521'] == gn_lineage_summary[3]


def test_summarize_no_taxonomy_fail(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(ValueError) as exc:
        c.run_sourmash('tax', 'summarize', g_csv)
    assert "error: the following arguments are required: -t/--taxonomy-csv" in str(exc.value)


def test_summarize_no_rank_lineage_summary(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(ValueError) as exc:
        runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'lineage_summary')
    assert "Rank (--rank) is required for krona and lineage_summary output formats." in str(exc.value)


def test_summarize_no_rank_krona(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(ValueError) as exc:
        runtmp.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'krona')
    assert "Rank (--rank) is required for krona and lineage_summary output formats." in str(exc.value)


def test_classify_no_rank_krona(runtmp):
    c = runtmp
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"

    with pytest.raises(ValueError) as exc:
        runtmp.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', tax, '-o', csv_base, '--output-format', 'krona')
    assert "Rank (--rank) is required for krona output format." in str(exc.value)


def test_summarize_duplicated_taxonomy_fail(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(Exception) as exc:
        c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', duplicated_csv)
        assert str(exc.value == "multiple lineages for identifier GCF_001881345")


def test_summarize_duplicated_taxonomy_force(runtmp):
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', duplicated_csv, '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    # same as stdout test - just check the first few lines
    assert c.last_result.status == 0
    assert "rank,fraction,lineage" in c.last_result.out
    assert 'superkingdom,0.131,d__Bacteria' in c.last_result.out
    assert "phylum,0.073,d__Bacteria;p__Bacteroidota" in c.last_result.out
    assert "phylum,0.058,d__Bacteria;p__Proteobacteria" in c.last_result.out


def test_summarize_missing_taxonomy(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        subset.write("\n".join(tax[:4]))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', subset_csv)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_003471795" in c.last_result.err
    assert "rank,fraction,lineage" in c.last_result.out

    assert "superkingdom,0.124,d__Bacteria" in c.last_result.out
    assert "phylum,0.066,d__Bacteria;p__Bacteroidota" in c.last_result.out
    assert "phylum,0.058,d__Bacteria;p__Proteobacteria" in c.last_result.out
    assert "class,0.066,d__Bacteria;p__Bacteroidota;c__Bacteroidia" in c.last_result.out
    assert "class,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria" in c.last_result.out
    assert "order,0.066,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales" in c.last_result.out


def test_summarize_missing_taxonomy_fail(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        subset.write("\n".join(tax[:4]))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'summarize', g_csv, '--taxonomy-csv', subset_csv, '--fail-on-missing-taxonomy', fail_ok=True)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)
    assert "The following are missing from the taxonomy information: GCF_003471795" in c.last_result.err
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in c.last_result.err
    assert c.last_result.status == -1


def test_classify_rank_stdout_0(runtmp):
    # test basic classify
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')

    c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_rank_csv_0(runtmp):
    # test basic classify - output csv
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csv_base = "out"
    cl_csv = csv_base + ".classifications.csv"
    csvout = runtmp.output(cl_csv)
    print("csvout: ", csvout)

    c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', tax,
                   '--rank', 'species', '-o', csv_base)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    cl_results = [x.rstrip() for x in open(csvout)]
    assert "query_name,rank,fraction,lineage" in cl_results[0]
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in cl_results[1]


def test_classify_gather_from_file_rank(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'classify', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_gather_from_file_two_files(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')

    # make test2 results (identical to test1 except query_name)
    g_res2 = runtmp.output("test2.gather.csv")
    test2_results = [x.replace("test1", "test2") for x in open(g_res, 'r')]
    with open(g_res2, 'w') as fp:
        for line in test2_results:
            fp.write(line)

    # write test1 and test2 files to a text file for input
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res2}\n")

    c.run_sourmash('tax', 'classify', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out
    assert "test2,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_gather_from_file_duplicate(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'classify', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_gather_cli_and_from_file(runtmp):
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

    c.run_sourmash('tax', 'classify', g_res, '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out
    assert "test2,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_gather_cli_and_from_file_duplicate(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")

    # also write test1 csv to a text file for input
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'classify', g_res, '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--rank', 'species')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_gather_from_file_threshold_0(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    g_res = utils.get_test_data('tax/test1.gather.csv')
    g_from_file = runtmp.output("tmp-from-file.txt")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{g_res}\n")

    c.run_sourmash('tax', 'classify', '--from-file', g_from_file, '--taxonomy-csv', taxonomy_csv,
                   '--containment-threshold', '0')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_rank_duplicated_taxonomy_fail(runtmp):
    # test basic summarize
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(Exception) as exc:
        c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', duplicated_csv,
                       '--rank', 'species')
        assert str(exc.value == "multiple lineages for identifier GCF_001881345")


def test_classify_rank_duplicated_taxonomy_force(runtmp):
    # test basic summarize
    c = runtmp
    # write temp taxonomy with duplicates
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', duplicated_csv,
                   '--rank', 'species', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_classify_missing_taxonomy_ignore_threshold(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', subset_csv, '--containment-threshold', '0')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in c.last_result.out


def test_classify_missing_taxonomy_ignore_rank(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', subset_csv, '--rank', 'species')
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert "query_name,rank,fraction,lineage" in c.last_result.out
    assert "test1,species,0.057,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in c.last_result.out


def test_classify_missing_taxonomy_fail_threshold(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', subset_csv,
                       '--fail-on-missing-taxonomy', '--containment-threshold', '0', fail_ok=True)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in c.last_result.err
    assert c.last_result.status == -1


def test_classify_missing_taxonomy_fail_rank(runtmp):
    c = runtmp
    # write temp taxonomy with missing entry
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    subset_csv = runtmp.output("subset_taxonomy.csv")
    with open(subset_csv, 'w') as subset:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax = [tax[0]] + tax[2:] # remove the best match (1st tax entry)
        subset.write("\n".join(tax))

    g_csv = utils.get_test_data('tax/test1.gather.csv')

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', g_csv, '--taxonomy-csv', subset_csv,
                       '--fail-on-missing-taxonomy', '--rank', 'species', fail_ok=True)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert "The following are missing from the taxonomy information: GCF_001881345" in c.last_result.err
    assert "Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy." in c.last_result.err
    assert c.last_result.status == -1


def test_classify_empty_gather_results_with_header_single(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results = [x for x in open(g_csv, 'r')]
    empty_tax_with_header = runtmp.output('tax_header.csv')
    # write temp empty gather results (header only)
    with open(empty_tax_with_header, "w") as fp:
        fp.write(gather_results[0])

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', empty_tax_with_header, '--taxonomy-csv', taxonomy_csv, fail_ok=True)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f'No gather results loaded from {empty_tax_with_header}.' in c.last_result.err
    assert 'Exiting.' in c.last_result.err


def test_classify_empty_gather_results_single(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results
    empty_tax = runtmp.output('tax_header.csv')
    with open(empty_tax, "w") as fp:
        fp.write("")

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', empty_tax, '--taxonomy-csv', taxonomy_csv, fail_ok=True)


    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f'No gather results loaded from {empty_tax}.' in c.last_result.err
    assert 'Exiting.' in c.last_result.err


def test_classify_empty_gather_results_single_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results (header only)
    empty_tax = runtmp.output('tax_header.csv')
    with open(empty_tax, "w") as fp:
        fp.write("")

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', empty_tax, '--taxonomy-csv', taxonomy_csv,
                       '--force', fail_ok=True)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f'--force is set. Attempting to continue to next set of gather results.' in c.last_result.err
    assert f'No results for classification. Exiting.' in c.last_result.err


def test_classify_empty_gather_results_with_empty_csv_force(runtmp):
    c = runtmp
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    # write temp empty gather results
    empty_tax = runtmp.output('tax_empty.txt')
    with open(empty_tax, "w") as fp:
        fp.write("")

    g_from_file = runtmp.output("tmp-from-csv.csv")
    with open(g_from_file, 'w') as f_csv:
        f_csv.write(f"{empty_tax}\n")

    with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
        c.run_sourmash('tax', 'classify', empty_tax, '--from-file', g_from_file,
                       '--taxonomy-csv', taxonomy_csv, '--rank', 'species', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == -1
    assert f'--force is set. Attempting to continue to next set of gather results.' in c.last_result.err
    assert 'No results for classification. Exiting.' in c.last_result.err


def test_classify_empty_gather_results_with_csv_force(runtmp):
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

    #with pytest.raises(ValueError) as exc: # should fail_ok handle this instead? Why ValueError?
    c.run_sourmash('tax', 'classify', empty_tax, '--from-file', g_from_file,
                   '--taxonomy-csv', taxonomy_csv, '--rank', 'species', '--force')

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert f'--force is set. Attempting to continue to next set of gather results.' in c.last_result.err
    assert f'loaded 1 gather files for classification' in c.last_result.err
    assert "test1,species,0.058,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in c.last_result.out


def test_label_0(runtmp):
    # test label
    c = runtmp

    g_csv = utils.get_test_data('tax/test1.gather.csv')
    tax = utils.get_test_data('tax/test.taxonomy.csv')
    csvout = runtmp.output("test1.gather.with-lineages.csv")
    out_dir = os.path.dirname(csvout)

    c.run_sourmash('tax', 'label', g_csv, '--taxonomy-csv', tax, '-o', out_dir)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    lin_gather_results = [x.rstrip() for x in open(csvout)]
    print("\n".join(lin_gather_results))

    assert "lineage" in lin_gather_results[0]
    assert "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli" in lin_gather_results[1]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[2]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus" in lin_gather_results[3]
    assert "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri" in lin_gather_results[4]

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
#def test_classify_bad_gather_results():
#    pass
#def test_classify_bad_lineage_input():
#    pass

