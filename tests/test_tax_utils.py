"""
Tests for functions in taxonomy submodule.
"""
import pytest

import sourmash
import sourmash_tst_utils as utils

from sourmash.tax import tax_utils
from sourmash.tax.tax_utils import (ascending_taxlist, get_ident, load_gather_results,
                                    summarize_gather_at, find_missing_identities,
                                    write_summary, load_gather_files_from_csv,
                                    write_classifications,
                                    make_krona_header, format_for_krona, write_krona,
                                    agg_sumgather_csvs_by_lineage, write_lineage_sample_frac)

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import LineagePair

from sourmash.lca.command_index import load_taxonomy_assignments

# utility functions for testing
def make_mini_gather_results(g_infolist):
    # make mini gather_results
    min_header = ["name","match_ident","f_unique_weighted"]
    gather_results = []
    for g_info in g_infolist:
        inf = dict(zip(min_header, g_info))
        gather_results.append(inf)
    return gather_results


def make_mini_taxonomy(tax_info):
    #pass in list of tuples: (name, lineage)
    taxD = {}
    for (name,lin) in tax_info:
        taxD[name] = lca_utils.make_lineage(lin)
    return taxD


## tests
def test_ascending_taxlist_1():
    assert list(ascending_taxlist()) ==  ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def test_ascending_taxlist_2():
    assert list(ascending_taxlist(include_strain=False)) ==  ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def test_get_ident_default():
    ident = "GCF_001881345.1"
    n_id = tax_utils.get_ident(ident)
    assert n_id == "GCF_001881345"


def test_get_ident_split_but_keep_version():
    ident = "GCF_001881345.1"
    n_id = tax_utils.get_ident(ident, keep_identifier_versions=True)
    assert n_id == "GCF_001881345.1"


def test_get_ident_no_split():
    ident = "GCF_001881345.1 secondname"
    n_id = tax_utils.get_ident(ident, split_identifiers=False)
    assert n_id == "GCF_001881345.1 secondname"


def test_load_gatherfiles_from_csv():
    from_csv = utils.get_test_data('tax/from-csv.csv')
    gather_files, seen_idents = load_gather_files_from_csv(from_csv)
    print("gather_files: ", gather_files)
    assert len(gather_files) == 1
    assert gather_files == [('test1', 'test1.gather.csv')]
    assert "test1" in seen_idents


def test_load_gather_results():
    gather_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results = tax_utils.load_gather_results(gather_csv)
    assert len(gather_results) == 4


# this function is in lca.command_index for now, but not tested there
def test_load_taxonomy_assignments():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign, num_rows = load_taxonomy_assignments(taxonomy_csv)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1']
    assert num_rows == 4 # should have read 4 rows


def test_load_taxonomy_assignments_split_id():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign, num_rows = load_taxonomy_assignments(taxonomy_csv, split_identifiers=True)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345', 'GCF_009494285', 'GCF_013368705', 'GCF_003471795']
    assert num_rows == 4 # should have read 4 rows


def test_load_taxonomy_assignments_with_ncbi_id(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    upd_csv = runtmp.output("updated_taxonomy.csv")
    with open(upd_csv, 'w') as new_tax:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        ncbi_id = "ncbi_id after_space"
        fake_lin = [ncbi_id] + ["sk", "phy", "cls", "ord", "fam", "gen", "sp"]
        ncbi_tax = ",".join(fake_lin)
        tax.append(ncbi_tax)
        new_tax.write("\n".join(tax))

    tax_assign, num_rows = load_taxonomy_assignments(upd_csv)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', "ncbi_id after_space"]
    assert num_rows == 5  # should have read 5 rows


def test_load_taxonomy_assignments_split_id_ncbi(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    upd_csv = runtmp.output("updated_taxonomy.csv")
    with open(upd_csv, 'w') as new_tax:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        ncbi_id = "ncbi_id after_space"
        fake_lin = [ncbi_id] + ["sk", "phy", "cls", "ord", "fam", "gen", "sp"]
        ncbi_tax = ",".join(fake_lin)
        tax.append(ncbi_tax)
        new_tax.write("\n".join(tax))

    tax_assign, num_rows = load_taxonomy_assignments(upd_csv, split_identifiers=True)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345', 'GCF_009494285', 'GCF_013368705', 'GCF_003471795', "ncbi_id"]
    assert num_rows == 5 # should have read 5 rows


def test_load_taxonomy_assignments_duplicate(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    with pytest.raises(Exception) as exc:
        tax_assign, num_rows = load_taxonomy_assignments(duplicated_csv)
        assert str(exc.value == "multiple lineages for identifier GCF_001881345.1")


def test_load_taxonomy_assignments_duplicate_force(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    # now force
    tax_assign, num_rows = load_taxonomy_assignments(duplicated_csv, force=True)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1']
    assert num_rows == 5 # should have read 5 rows


def test_find_missing_identities():
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    n, ids = find_missing_identities(g_res, taxD)
    print("n_missing: ", n)
    print("ids_missing: ", ids)
    assert n == 1
    assert ids == ["gB"]


def test_summarize_gather_at_0():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 1.0)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.0)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.5)]


def test_summarize_gather_at_1():
    """test two matches, diff f_unique_weighted"""
    # make mini gather_results
    gA = ["gA","0.5","0.6"]
    gB = ["gB","0.3","0.1"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 0.7)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),0.7)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.6),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.1)]


def test_summarize_gather_at_over100percent_f_unique_weighted():
    """gather matches that add up to >100% f_unique_weighted"""
    ## @NTP:  currently passes, we should probably make this fail
    # make mini gather_results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.6"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 1.1)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.1)]
    cl_sum = summarize_gather_at("class", taxD, g_res)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='d')),0.6),
                      ((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5)]


def test_summarize_gather_at_missing_ignore():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res, skip_idents=['gB'])
    print("sk_sum: ", sk_sum)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 0.5)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res, skip_idents=['gB'])
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),0.5)]
    cl_sum = summarize_gather_at("class", taxD, g_res, skip_idents=['gB'])
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5)]


def test_summarize_gather_at_missing_fail():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    # run summarize_gather_at and check results!
    with pytest.raises(KeyError) as exc:
        sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
        assert exc.value == "ident gB is not in the taxonomy database."


def test_summarize_gather_at_best_only_0():
    """test two matches, diff f_unique_weighted"""
    # make mini gather_results
    gA = ["gA","0.5","0.6"]
    gB = ["gB","0.3","0.1"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res, best_only=True)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 0.7)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res, best_only=True)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),0.7)]
    cl_sum = summarize_gather_at("class", taxD, g_res, best_only=True)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.6)]


def test_summarize_gather_at_best_only_equal_choose_first():
    """test two matches, equal f_unique_weighted. best_only chooses first"""
    # make mini gather_results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res, best_only=True)
    assert sk_sum == [((LineagePair(rank='superkingdom', name='a'),), 1.0)]
    phy_sum = summarize_gather_at("phylum", taxD, g_res, best_only=True)
    assert phy_sum == [((LineagePair(rank='superkingdom', name='a'),
                         LineagePair(rank='phylum', name='b')),1.0)]
    cl_sum = summarize_gather_at("class", taxD, g_res, best_only=True)
    assert cl_sum == [((LineagePair(rank='superkingdom', name='a'),
                        LineagePair(rank='phylum', name='b'),
                        LineagePair(rank='class', name='c')),0.5)]


def test_write_summary_csv(runtmp):
    """test summary csv write function"""

    sum_gather = {'superkingdom': [((LineagePair(rank='superkingdom', name='a'),), 1.0)],
                  'phylum': [((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='b')), 1.0)]}

    outs= runtmp.output("outsum.csv")
    with open(outs, 'w') as out_fp:
        write_summary(sum_gather, out_fp)

    sr = [x.rstrip().split(',') for x in open(outs, 'r')]
    print("gather_summary_results_from_file: \n", sr)
    assert sr[0] ==  ['rank', 'fraction', 'lineage']
    assert sr[1] ==  ['superkingdom', '1.000', 'a']
    assert sr[2] ==  ['phylum', '1.000', 'a;b']


def test_write_classification_csv(runtmp):
    """test classification csv write function"""

    classif = {'superkingdom': [("x",((LineagePair(rank='superkingdom', name='a'),), 1.0))],
                  'phylum': [("y", ((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='b')), 1.0))]}

    outc= runtmp.output("outclass.csv")
    with open(outc, 'w') as out_fp:
        write_classifications(classif, out_fp)

    cr = [x.rstrip().split(',') for x in open(outc, 'r')]
    print("classification_summary_results_from_file: \n", cr)
    assert cr[0] == ['query_name', 'classification_rank', 'fraction_matched_at_rank', 'lineage']
    assert cr[1] == ['superkingdom', 'x', '1.000', 'a']
    assert cr[2] == ['phylum', 'y', '1.000', 'a;b']


def test_make_krona_header_0():
    hd = make_krona_header("species")
    print("header: ", hd)
    assert hd == ("fraction", "superkingdom", "phylum", "class", "order", "family", "genus", "species")


def test_make_krona_header_1():
    hd = make_krona_header("order")
    print("header: ", hd)
    assert hd == ("fraction", "superkingdom", "phylum", "class", "order")


def test_make_krona_header_strain():
    hd = make_krona_header("strain", include_strain=True)
    print("header: ", hd)
    assert hd == ("fraction", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain")


def test_make_krona_header_fail():
    with pytest.raises(ValueError) as exc:
        hd = make_krona_header("strain")
        assert str(exc.value) == "Rank strain not present in available ranks"


def test_format_for_krona_0():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # check krona format and check results!
    sk_sum = summarize_gather_at("superkingdom", taxD, g_res)
    krona_res = format_for_krona("superkingdom", {"superkingdom": sk_sum})
    print("krona_res: ", krona_res)
    assert krona_res == [(1.0, 'a')]

    phy_sum = summarize_gather_at("phylum", taxD, g_res)
    krona_res = format_for_krona("phylum", {"phylum": phy_sum})
    print("krona_res: ", krona_res)
    assert krona_res == [(1.0, 'a', 'b')]


def test_format_for_krona_1():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # summarize with all ranks
    sum_res = {}
    #for rank in lca_utils.taxlist(include_strain=False):
    for rank in ['superkingdom', 'phylum', 'class']:
        sum_res[rank] = summarize_gather_at(rank, taxD, g_res)
    print('summarized gather: ', sum_res)
    # check krona format
    sk_krona = format_for_krona("superkingdom", sum_res)
    print("sk_krona: ", sk_krona)
    assert sk_krona == [(1.0, 'a')]
    phy_krona = format_for_krona("phylum", sum_res)
    print("phy_krona: ", phy_krona)
    assert phy_krona ==  [(1.0, 'a', 'b')]
    cl_krona = format_for_krona("class", sum_res)
    print("cl_krona: ", cl_krona)
    assert cl_krona ==  [(0.5, 'a', 'b', 'c'), (0.5, 'a', 'b', 'd')]


def test_format_for_krona_best_only():
    """test two matches, equal f_unique_weighted"""
    # make gather results
    gA = ["gA","0.5","0.5"]
    gB = ["gB","0.3","0.5"]
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # summarize with all ranks
    sum_res = {}
    #for rank in lca_utils.taxlist(include_strain=False):
    for rank in ['superkingdom', 'phylum', 'class']:
        sum_res[rank] = summarize_gather_at(rank, taxD, g_res, best_only=True)
    print('summarized gather: ', sum_res)
    # check krona format
    sk_krona = format_for_krona("superkingdom", sum_res)
    print("sk_krona: ", sk_krona)
    assert sk_krona == [(1.0, 'a')]
    phy_krona = format_for_krona("phylum", sum_res)
    print("phy_krona: ", phy_krona)
    assert phy_krona ==  [(1.0, 'a', 'b')]
    cl_krona = format_for_krona("class", sum_res)
    print("cl_krona: ", cl_krona)
    assert cl_krona ==  [(0.5, 'a', 'b', 'c')]


def test_write_krona(runtmp):
    """test two matches, equal f_unique_weighted"""
    class_krona_results =  [(0.5, 'a', 'b', 'c'), (0.5, 'a', 'b', 'd')]
    outk= runtmp.output("outkrona.tsv")
    with open(outk, 'w') as out_fp:
        write_krona("class", class_krona_results, out_fp)

    kr = [x.strip().split('\t') for x in open(outk, 'r')]
    print("krona_results_from_file: \n", kr)
    assert kr[0] == ["fraction", "superkingdom", "phylum", "class"]
    assert kr[1] == ["0.5", "a", "b", "c"]
    assert kr[2] == ["0.5", "a", "b", "d"]


def test_agg_sumgather_csvs_by_lineage(runtmp):
    # some summarized gather dicts
    sum_gather1 = {'superkingdom': [((LineagePair(rank='superkingdom', name='a'),), 0.5)],
                  'phylum': [((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='b')), 0.5)]}

    sum_gather2 = {'superkingdom': [((LineagePair(rank='superkingdom', name='a'),), 0.7)],
                  'phylum': [((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='c')), 0.7)]}

    # write summarized gather results csvs
    sg1= runtmp.output("sample1.csv")
    with open(sg1, 'w') as out_fp:
        write_summary(sum_gather1, out_fp)

    sg2= runtmp.output("sample2.csv")
    with open(sg2, 'w') as out_fp:
        write_summary(sum_gather2, out_fp)

    # test agg_summarized_gather_csvs_by_lineage_at_rank
    linD, sample_names = agg_sumgather_csvs_by_lineage([sg1,sg2], rank="phylum")
    print("lineage dict: \n", linD)
    assert linD == {'a;b': {'sample1': '0.500', 'sample2': 0.0}, 'a;c': {'sample1': 0.0, 'sample2': '0.700'}}
    assert sample_names == ['sample1', 'sample2']
    linD, sample_names = agg_sumgather_csvs_by_lineage([sg1,sg2], rank="superkingdom")
    print("lineage dict: \n", linD)
    assert linD == {'a': {'sample1': '0.500' ,'sample2': '0.700'}}
    assert sample_names == ['sample1', 'sample2']


def test_write_lineage_sample_frac(runtmp):
    outfrac = runtmp.output('outfrac.csv')
    sample_names = ['sample1', 'sample2']
    sk_linD = {'a': {'sample1': '0.500' ,'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, sk_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a', '0.500', '0.700']]

    phy_linD = {'a;b': {'sample1': '0.500', 'sample2': '0'}, 'a;c': {'sample1': '0', 'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, phy_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a;b', '0.500', '0'],  ['a;c', '0', '0.700']]

def test_agg_sumgather_csvs_by_lineage_improper_rank(runtmp):
    # some summarized gather dicts
    sum_gather1 = {'superkingdom': [((LineagePair(rank='superkingdom', name='a'),), 0.5)],
                  'phylum': [((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='b')), 0.5)]}
    sum_gather2 = {'superkingdom': [((LineagePair(rank='superkingdom', name='a'),), 0.7)],
                  'phylum': [((LineagePair(rank='superkingdom', name='a'),
                               LineagePair(rank='phylum', name='c')), 0.7)]}

    # write summarized gather results csvs
    sg1= runtmp.output("sample1.csv")
    with open(sg1, 'w') as out_fp:
        write_summary(sum_gather1, out_fp)

    sg2= runtmp.output("sample2.csv")
    with open(sg2, 'w') as out_fp:
        write_summary(sum_gather2, out_fp)

    # test agg_summarized_gather_csvs_by_lineage_at_rank
    with pytest.raises(ValueError) as exc:
        linD, sample_names = agg_sumgather_csvs_by_lineage([sg1,sg2], rank="strain")
        print("ValueError: ", exc.value)
        assert exc.value == "Rank strain not available."
