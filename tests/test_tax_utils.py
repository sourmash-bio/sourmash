"""
Tests for functions in taxonomy submodule.
"""

import pytest
import os
from os.path import basename
import gzip

import sourmash_tst_utils as utils

from sourmash.tax.tax_utils import (ascending_taxlist, get_ident, load_gather_results,
                                    summarize_gather_at, find_missing_identities,
                                    write_summary, MultiLineageDB,
                                    collect_gather_csvs, check_and_load_gather_csvs,
                                    SummarizedGatherResult, ClassificationResult,
                                    QueryInfo, GatherRow, TaxResult, QueryTaxResult,
                                    BaseLineageInfo, RankLineageInfo, LineagePair,
                                    write_classifications,
                                    aggregate_by_lineage_at_rank,
                                    make_krona_header, format_for_krona, write_krona,
                                    combine_sumgather_csvs_by_lineage, write_lineage_sample_frac,
                                    LineageDB, LineageDB_Sqlite,
                                    SumGathInf, ClassInf, QInfo)

# import lca utils as needed for now
from sourmash.lca import lca_utils

# utility functions for testing
def make_mini_gather_results(g_infolist, include_ksize_and_scaled=False):
    # make mini gather_results
    min_header = ["query_name", "name", "match_ident", "f_unique_to_query", "query_md5", "query_filename", "f_unique_weighted", "unique_intersect_bp", "remaining_bp"]
    if include_ksize_and_scaled:
        min_header.extend(['ksize', 'scaled'])
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


def make_GatherRow(gather_dict=None, exclude_cols=[]):
    """Load artificial gather row (dict) into GatherRow class"""
    # default contains just the essential cols
    gatherD = {'query_name': 'q1',
               'query_md5': 'md5',
               'query_filename': 'query_fn',
               'name': 'gA',
               'f_unique_weighted': 0.2,
               'f_unique_to_query': 0.1,
               'query_bp':100,
               'unique_intersect_bp': 20,
               'remaining_bp': 1,
               'ksize': 31,
               'scaled': 1}
    if gather_dict is not None:
        gatherD.update(gather_dict)
    for col in exclude_cols:
        gatherD.pop(col)
    gatherRaw = GatherRow(**gatherD)
    return gatherRaw


def make_TaxResult(gather_dict=None, taxD=None, keep_full_ident=False, keep_ident_version=False, skip_idents=None):
    """Make TaxResult from artificial gather row (dict)"""
    gRow = make_GatherRow(gather_dict)
    taxres = TaxResult(raw=gRow, keep_full_identifiers=keep_full_ident, keep_identifier_versions=keep_ident_version)
    if taxD is not None:
        taxres.get_match_lineage(tax_assignments=taxD, skip_idents=skip_idents)
    return taxres


def make_QueryTaxResults(gather_info, taxD=None, single_query=False, keep_full_ident=False, keep_ident_version=False,
                        skip_idents=None, summarize=False, classify=False, classify_rank=None, c_thresh=0.1, ani_thresh=None):
    """Make QueryTaxResult(s) from artificial gather information, formatted as list of gather rows (dicts)"""
    gather_results = {}
    this_querytaxres = None
    for gather_infoD in gather_info:
        taxres = make_TaxResult(gather_infoD, taxD=taxD,  keep_full_ident=keep_full_ident,
                                keep_ident_version=keep_ident_version, skip_idents=skip_idents)
        query_name = taxres.query_name
        # add to matching QueryTaxResult or create new one
        if not this_querytaxres or not this_querytaxres.is_compatible(taxres):
            # get existing or initialize new
            this_querytaxres = gather_results.get(query_name, QueryTaxResult(taxres.query_info))
        this_querytaxres.add_taxresult(taxres)
#        print('missed_ident?', taxres.missed_ident)
        gather_results[query_name] = this_querytaxres
    if summarize:
        for query_name, qres in gather_results.items():
            qres.build_summarized_result()
    if classify:
        for query_name, qres in gather_results.items():
            qres.build_classification_result(rank=classify_rank, containment_threshold=c_thresh, ani_threshold=ani_thresh)
    # for convenience: If working with single query, just return that QueryTaxResult.
    if single_query:
        if len(gather_results.keys()) > 1:
            raise ValueError("You passed in results for more than one query")
        else:
            return next(iter(gather_results.values()))
    return gather_results


## tests
def test_ascending_taxlist_1():
    assert list(ascending_taxlist()) ==  ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def test_ascending_taxlist_2():
    assert list(ascending_taxlist(include_strain=False)) ==  ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def test_GatherRow_old_gather():
    # gather does not contain query_name column
    gA = {"name": "gA.1 name"}
    with pytest.raises(TypeError) as exc:
        make_GatherRow(gA, exclude_cols=['query_bp'])
    print(str(exc))
    assert "__init__() missing 1 required positional argument: 'query_bp'" in str(exc)


def test_get_ident_default():
    ident = "GCF_001881345.1"
    n_id = get_ident(ident)
    assert n_id == "GCF_001881345"


def test_TaxResult_get_ident_default():
    gA = {"name": "GCF_001881345.1"}  # gather result with match name as GCF_001881345.1
    taxres = make_TaxResult(gA)
    print(taxres.match_ident)
    assert taxres.match_ident == "GCF_001881345"


def test_get_ident_split_but_keep_version():
    ident = "GCF_001881345.1 secondname"
    n_id = get_ident(ident, keep_identifier_versions=True)
    assert n_id == "GCF_001881345.1"


def test_TaxResult_get_ident_split_but_keep_version():
    gA = {"name": "GCF_001881345.1 secondname"}
    taxres = make_TaxResult(gA, keep_ident_version=True)
    print("raw ident: ", taxres.raw.name)
    print("keep_full?: ", taxres.keep_full_identifiers)
    print("keep_version?: ",taxres.keep_identifier_versions)
    print("final ident: ", taxres.match_ident)
    assert taxres.match_ident == "GCF_001881345.1"


def test_get_ident_no_split():
    ident = "GCF_001881345.1 secondname"
    n_id = get_ident(ident, keep_full_identifiers=True)
    assert n_id == "GCF_001881345.1 secondname"


def test_TaxResult_get_ident_keep_full():
    gA = {"name": "GCF_001881345.1 secondname"}
    taxres = make_TaxResult(gA, keep_full_ident=True)
    print("raw ident: ", taxres.raw.name)
    print("keep_full?: ", taxres.keep_full_identifiers)
    print("keep_version?: ",taxres.keep_identifier_versions)
    print("final ident: ", taxres.match_ident)
    assert taxres.match_ident == "GCF_001881345.1 secondname"


def test_collect_gather_csvs(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    from_file = runtmp.output("tmp-from-file.txt")
    with open(from_file, 'w') as fp:
        fp.write(f"{g_csv}\n")

    gather_files = collect_gather_csvs([g_csv], from_file=from_file)
    print("gather_files: ", gather_files)
    assert len(gather_files) == 1
    assert basename(gather_files[0]) == 'test1.gather.csv'


def test_check_and_load_gather_csvs_empty(runtmp):
    g_res = runtmp.output('empty.gather.csv')
    with open(g_res, 'w') as fp:
        fp.write("")
    csvs = [g_res]
    # load taxonomy csv
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv], keep_full_identifiers=1)

    print(tax_assign)
    # check gather results and missing ids
    with pytest.raises(Exception) as exc:
        gather_results, ids_missing, n_missing, header = check_and_load_gather_csvs(csvs, tax_assign)
    assert "Cannot read gather results from" in str(exc.value)


def test_check_and_load_gather_csvs_with_empty_force(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    #  make gather results with taxonomy name not in tax_assign
    g_res2 = runtmp.output('gA.gather.csv')
    g_results = [x.replace("GCF_001881345.1", "gA") for x in open(g_csv, 'r')]
    with open(g_res2, 'w') as fp:
        for line in g_results:
            fp.write(line)
    # make empty gather results
    g_res3 = runtmp.output('empty.gather.csv')
    with open(g_res3, 'w') as fp:
        fp.write("")

    csvs = [g_res2, g_res3]

    # load taxonomy csv
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    print(tax_assign)
    # check gather results and missing ids
    gather_results, ids_missing, n_missing, header = check_and_load_gather_csvs(csvs, tax_assign, force=True)
    assert len(gather_results) == 4
    print("n_missing: ", n_missing)
    print("ids_missing: ", ids_missing)
    assert n_missing == 1
    assert ids_missing == {"gA"}


def test_check_and_load_gather_lineage_csvs_empty(runtmp):
    # try loading an empty annotated gather file
    g_res = runtmp.output('empty.gather-tax.csv')
    with open(g_res, 'w') as fp:
        fp.write("")

    with pytest.raises(ValueError) as exc:
        tax_assign = LineageDB.load_from_gather_with_lineages(g_res)
    assert "cannot read taxonomy assignments" in str(exc.value)


def test_check_and_load_gather_lineage_csvs_bad_header(runtmp):
    # test on file with wrong headers
    g_res = runtmp.output('empty.gather-tax.csv')
    with open(g_res, 'w', newline="") as fp:
        fp.write("x,y,z")

    with pytest.raises(ValueError) as exc:
        tax_assign = LineageDB.load_from_gather_with_lineages(g_res)
    assert "Expected headers 'name' and 'lineage' not found. Is this a with-lineages file?" in str(exc.value)


def test_check_and_load_gather_lineage_csvs_dne(runtmp):
    # test loading with-lineage file that does not exist
    g_res = runtmp.output('empty.gather-tax.csv')

    with pytest.raises(ValueError) as exc:
        tax_assign = LineageDB.load_from_gather_with_lineages(g_res)
    assert "does not exist" in str(exc.value)


def test_check_and_load_gather_lineage_csvs_isdir(runtmp):
    # test loading a with-lineage file that is actually a directory
    g_res = runtmp.output('empty.gather-tax.csv')
    os.mkdir(g_res)

    with pytest.raises(ValueError) as exc:
        tax_assign = LineageDB.load_from_gather_with_lineages(g_res)
    assert "is a directory" in str(exc.value)


def test_check_and_load_gather_csvs_fail_on_missing(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')
    # make gather results with taxonomy name not in tax_assign
    g_res2 = runtmp.output('gA.gather.csv')
    g_results = [x.replace("GCF_001881345.1", "gA") for x in open(g_csv, 'r')]
    with open(g_res2, 'w') as fp:
        for line in g_results:
            fp.write(line)

    csvs = [g_res2]

    # load taxonomy csv
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv], keep_full_identifiers=1)
    print(tax_assign)
    # check gather results and missing ids
    with pytest.raises(ValueError) as exc:
        gather_results, ids_missing, n_missing, header = check_and_load_gather_csvs(csvs, tax_assign, fail_on_missing_taxonomy=True, force=True)
    assert "Failing on missing taxonomy" in str(exc)


def test_load_gather_results():
    gather_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results, header, seen_queries = load_gather_results(gather_csv)
    assert len(gather_results) == 4


def test_load_gather_results_gzipped(runtmp):
    gather_csv = utils.get_test_data('tax/test1.gather.csv')

    # rewrite gather_csv as gzipped csv
    gz_gather = runtmp.output('g.csv.gz')
    with open(gather_csv, 'rb') as f_in, gzip.open(gz_gather, 'wb') as f_out:
        f_out.writelines(f_in)
    gather_results, header, seen_queries = load_gather_results(gz_gather)
    assert len(gather_results) == 4


def test_load_gather_results_bad_header(runtmp):
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    bad_g_csv = runtmp.output('g.csv')

    #creates bad gather result
    bad_g = [x.replace("f_unique_to_query", "nope") for x in open(g_csv, 'r')]
    with open(bad_g_csv, 'w') as fp:
        for line in bad_g:
            fp.write(line)
    print("bad_gather_results: \n", bad_g)

    with pytest.raises(ValueError) as exc:
        gather_results, header = load_gather_results(bad_g_csv)
    assert f"Not all required gather columns are present in '{bad_g_csv}'." in str(exc.value)


def test_load_gather_results_empty(runtmp):
    empty_csv = runtmp.output('g.csv')

    #creates empty gather result
    with open(empty_csv, 'w') as fp:
        fp.write('')

    with pytest.raises(ValueError) as exc:
        gather_results, header = load_gather_results(empty_csv)
    assert f"Cannot read gather results from '{empty_csv}'. Is file empty?" in str(exc.value)


def test_load_taxonomy_csv():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv])
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1']
    assert len(tax_assign) == 6 # should have read 6 rows


def test_load_taxonomy_csv_gzip(runtmp):
    # test loading a gzipped taxonomy csv file
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_gz = runtmp.output('tax.csv.gz')

    with gzip.open(tax_gz, 'wt') as outfp:
        with open(taxonomy_csv, 'rt') as infp:
            data = infp.read()
        outfp.write(data)

    tax_assign = MultiLineageDB.load([tax_gz])
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1']
    assert len(tax_assign) == 6 # should have read 6 rows


def test_load_taxonomy_csv_split_id():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv], keep_full_identifiers=0,
                                     keep_identifier_versions=False)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345', 'GCF_009494285', 'GCF_013368705', 'GCF_003471795', 'GCF_000017325', 'GCF_000021665']
    assert len(tax_assign) == 6 # should have read 6 rows


def test_load_taxonomy_csv_with_ncbi_id(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    upd_csv = runtmp.output("updated_taxonomy.csv")
    with open(upd_csv, 'w') as new_tax:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        ncbi_id = "ncbi_id after_space"
        fake_lin = [ncbi_id] + ["sk", "phy", "cls", "ord", "fam", "gen", "sp"]
        ncbi_tax = ",".join(fake_lin)
        tax.append(ncbi_tax)
        new_tax.write("\n".join(tax))

    tax_assign = MultiLineageDB.load([upd_csv], keep_full_identifiers=True)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1', "ncbi_id after_space"]
    assert len(tax_assign) == 7  # should have read 7 rows


def test_load_taxonomy_csv_split_id_ncbi(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    upd_csv = runtmp.output("updated_taxonomy.csv")
    with open(upd_csv, 'w') as new_tax:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        ncbi_id = "ncbi_id after_space"
        fake_lin = [ncbi_id] + ["sk", "phy", "cls", "ord", "fam", "gen", "sp"]
        ncbi_tax = ",".join(fake_lin)
        tax.append(ncbi_tax)
        new_tax.write("\n".join(tax))

    tax_assign = MultiLineageDB.load([upd_csv], keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345', 'GCF_009494285', 'GCF_013368705', 'GCF_003471795', 'GCF_000017325', 'GCF_000021665', "ncbi_id"]
    assert len(tax_assign) == 7 # should have read 7 rows

    # check for non-sensical args.
    with pytest.raises(ValueError):
        tax_assign = MultiLineageDB.load([upd_csv], keep_full_identifiers=1,
                                         keep_identifier_versions=False)


def test_load_taxonomy_csv_duplicate(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1] + 'FOO') # add first tax_assign again
        print(tax[-1])
        dup.write("\n".join(tax))

    with pytest.raises(Exception) as exc:
        MultiLineageDB.load([duplicated_csv])

    assert "cannot read taxonomy assignments" in str(exc.value)
    assert "multiple lineages for identifier GCF_001881345.1" in str(exc.value)


def test_load_taxonomy_csv_duplicate_force(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    duplicated_csv = runtmp.output("duplicated_taxonomy.csv")
    with open(duplicated_csv, 'w') as dup:
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        tax.append(tax[1]) # add first tax_assign again
        dup.write("\n".join(tax))

    # now force
    tax_assign = MultiLineageDB.load([duplicated_csv], force=True)

    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1']


def test_find_missing_identities():
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    ids = find_missing_identities(g_res, taxD)
    print("ids_missing: ", ids)
    assert ids == {"gB"}


def test_summarize_gather_at_0():
    """test two matches, equal f_unique_to_query"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)

    # superkingdom
    assert len(sk_sum) == 1
    print("superkingdom summarized gather: ", sk_sum[0])
    assert sk_sum[0].query_name == "queryA"
    assert sk_sum[0].query_md5 == "queryA_md5"
    assert sk_sum[0].query_filename == "queryA.sig"
    assert sk_sum[0].rank == 'superkingdom'
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 1.0
    assert sk_sum[0].f_weighted_at_rank == 1.0
    assert sk_sum[0].bp_match_at_rank == 100

    # phylum
    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res)
    print("phylum summarized gather: ", phy_sum[0])
    assert len(phy_sum) == 1
    assert phy_sum[0].query_name == "queryA"
    assert phy_sum[0].query_md5 == "queryA_md5"
    assert phy_sum[0].query_filename == "queryA.sig"
    assert phy_sum[0].rank == 'phylum'
    assert phy_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),lca_utils.LineagePair(rank='phylum', name='b'))
    assert phy_sum[0].fraction == 1.0
    assert phy_sum[0].f_weighted_at_rank == 1.0
    assert phy_sum[0].bp_match_at_rank == 100
    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res)
    assert len(cl_sum) == 2
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].query_name == "queryA"
    assert cl_sum[0].query_md5 == "queryA_md5"
    assert cl_sum[0].query_filename == "queryA.sig"
    assert cl_sum[0].rank == 'class'
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[0].fraction == 0.5
    assert cl_sum[0].f_weighted_at_rank == 0.5
    assert cl_sum[0].bp_match_at_rank == 50
    assert cl_sum[1].rank == 'class'
    assert cl_sum[1].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='d'))
    assert cl_sum[1].fraction == 0.5
    assert cl_sum[1].f_weighted_at_rank == 0.5
    assert cl_sum[1].bp_match_at_rank == 50


def test_summarize_gather_at_1():
    """test two matches, diff f_unique_to_query"""
    # make mini gather_results
    ksize=31
    scaled=10
    gA = ["queryA", "gA","0.5","0.6", "queryA_md5", "queryA.sig", '0.5', '60', '40', ksize, scaled]
    gB = ["queryA", "gB","0.3","0.1", "queryA_md5", "queryA.sig", '0.1', '10', '90', ksize, scaled]
    g_res = make_mini_gather_results([gA,gB], include_ksize_and_scaled=True)

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res, estimate_query_ani=True)

    # superkingdom
    assert len(sk_sum) == 2
    print("\nsuperkingdom summarized gather 0: ", sk_sum[0])
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 0.7
    assert sk_sum[0].bp_match_at_rank == 70
    print("superkingdom summarized gather 1: ", sk_sum[1])
    assert sk_sum[1].lineage == ()
    assert round(sk_sum[1].fraction, 1) == 0.3
    assert sk_sum[1].bp_match_at_rank == 30
    assert sk_sum[0].query_ani_at_rank == 0.9885602934376099
    assert sk_sum[1].query_ani_at_rank == None

    # phylum
    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res, estimate_query_ani=False)
    print("phylum summarized gather 0: ", phy_sum[0])
    print("phylum summarized gather 1: ", phy_sum[1])
    assert len(phy_sum) == 2
    assert phy_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),lca_utils.LineagePair(rank='phylum', name='b'))
    assert phy_sum[0].fraction == 0.7
    assert phy_sum[0].f_weighted_at_rank == 0.6
    assert phy_sum[0].bp_match_at_rank == 70
    assert phy_sum[1].lineage == ()
    assert round(phy_sum[1].fraction, 1) == 0.3
    assert phy_sum[1].bp_match_at_rank == 30
    assert phy_sum[0].query_ani_at_rank == None
    assert phy_sum[1].query_ani_at_rank == None
    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res, estimate_query_ani=True)
    assert len(cl_sum) == 3
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[0].fraction == 0.6
    assert cl_sum[0].f_weighted_at_rank == 0.5
    assert cl_sum[0].bp_match_at_rank == 60
    assert cl_sum[0].query_ani_at_rank == 0.9836567776983505

    assert cl_sum[1].rank == 'class'
    assert cl_sum[1].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='d'))
    assert cl_sum[1].fraction == 0.1
    assert cl_sum[1].f_weighted_at_rank == 0.1
    assert cl_sum[1].bp_match_at_rank == 10
    assert cl_sum[1].query_ani_at_rank == 0.9284145445194744
    assert cl_sum[2].lineage == ()
    assert round(cl_sum[2].fraction, 1) == 0.3
    assert cl_sum[2].query_ani_at_rank == None


def test_summarize_gather_at_perfect_match():
    """test 100% gather match (f_unique_to_query == 1)"""
    # make mini gather_results
    gA = ["queryA", "gA","0.5","1.0", "queryA_md5", "queryA.sig", '0.5', '100', '0']
    gB = ["queryA", "gB","0.3","0.0", "queryA_md5", "queryA.sig", '0.5', '0', '100']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)
    # superkingdom
    assert len(sk_sum) == 1
    print("superkingdom summarized gather: ", sk_sum[0])
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 1.0


def test_summarize_gather_at_over100percent_f_unique_to_query():
    """gather matches that add up to >100% f_unique_to_query"""
    # make mini gather_results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.6", "queryA_md5", "queryA.sig", '0.5', '60', '40']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    with pytest.raises(ValueError) as exc:
        sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)
    assert "The tax summary of query 'queryA' is 1.1, which is > 100% of the query!!" in str(exc)

    # phylum
    with pytest.raises(ValueError) as exc:
        phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res)
    assert "The tax summary of query 'queryA' is 1.1, which is > 100% of the query!!" in str(exc)

    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res)
    assert len(cl_sum) == 2
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='d'))
    assert cl_sum[0].fraction == 0.6
    assert cl_sum[0].bp_match_at_rank == 60
    assert cl_sum[1].rank == 'class'
    assert cl_sum[1].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[1].fraction == 0.5
    assert cl_sum[1].bp_match_at_rank == 50


def test_summarize_gather_at_missing_ignore():
    """test two matches, ignore missing taxonomy"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    # run summarize_gather_at and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res, skip_idents=['gB'])
    # superkingdom
    assert len(sk_sum) == 2
    print("superkingdom summarized gather: ", sk_sum[0])
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 0.5
    assert sk_sum[0].bp_match_at_rank == 50
    assert sk_sum[1].lineage == ()
    assert sk_sum[1].fraction == 0.5
    assert sk_sum[1].bp_match_at_rank == 50

    # phylum
    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res, skip_idents=['gB'])
    print("phylum summarized gather: ", phy_sum[0])
    assert len(phy_sum) == 2
    assert phy_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),lca_utils.LineagePair(rank='phylum', name='b'))
    assert phy_sum[0].fraction == 0.5
    assert phy_sum[0].bp_match_at_rank == 50
    assert phy_sum[1].lineage == ()
    assert phy_sum[1].fraction == 0.5
    assert phy_sum[1].bp_match_at_rank == 50
    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res, skip_idents=['gB'])
    assert len(cl_sum) == 2
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[0].fraction == 0.5
    assert cl_sum[0].bp_match_at_rank == 50
    assert cl_sum[1].lineage == ()
    assert cl_sum[1].fraction == 0.5
    assert cl_sum[1].bp_match_at_rank == 50


def test_summarize_gather_at_missing_fail():
    """test two matches, fail on missing taxonomy"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    # run summarize_gather_at and check results!
    with pytest.raises(ValueError) as exc:
        sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)
    assert "ident gB is not in the taxonomy database." in str(exc.value)


def test_summarize_gather_at_best_only_0():
    """test two matches, diff f_unique_to_query"""
    # make mini gather_results
    ksize =31
    scaled=10
    gA = ["queryA", "gA","0.5","0.6", "queryA_md5", "queryA.sig", '0.5', '60', '40', ksize, scaled]
    gB = ["queryA", "gB","0.3","0.1", "queryA_md5", "queryA.sig", '0.5', '10', '90', ksize, scaled]
    g_res = make_mini_gather_results([gA,gB],include_ksize_and_scaled=True)

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])
    # run summarize_gather_at and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res, best_only=True,estimate_query_ani=True)
    # superkingdom
    assert len(sk_sum) == 1
    print("superkingdom summarized gather: ", sk_sum[0])
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 0.7
    assert sk_sum[0].bp_match_at_rank == 70
    print("superk ANI:",sk_sum[0].query_ani_at_rank)
    assert sk_sum[0].query_ani_at_rank == 0.9885602934376099

    # phylum
    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res, best_only=True,estimate_query_ani=True)
    print("phylum summarized gather: ", phy_sum[0])
    assert len(phy_sum) == 1
    assert phy_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),lca_utils.LineagePair(rank='phylum', name='b'))
    assert phy_sum[0].fraction == 0.7
    assert phy_sum[0].bp_match_at_rank == 70
    print("phy ANI:",phy_sum[0].query_ani_at_rank)
    assert phy_sum[0].query_ani_at_rank == 0.9885602934376099
    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res, best_only=True, estimate_query_ani=True)
    assert len(cl_sum) == 1
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[0].fraction == 0.6
    assert cl_sum[0].bp_match_at_rank == 60
    print("cl ANI:",cl_sum[0].query_ani_at_rank)
    assert cl_sum[0].query_ani_at_rank == 0.9836567776983505


def test_summarize_gather_at_best_only_equal_choose_first():
    """test two matches, equal f_unique_to_query. best_only chooses first"""
    # make mini gather_results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # run summarize_gather_at and check results!
    # class
    cl_sum, _, _ = summarize_gather_at("class", taxD, g_res, best_only=True)
    assert len(cl_sum) == 1
    print("class summarized gather: ", cl_sum)
    assert cl_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                 lca_utils.LineagePair(rank='phylum', name='b'),
                                 lca_utils.LineagePair(rank='class', name='c'))
    assert cl_sum[0].fraction == 0.5
    assert cl_sum[0].bp_match_at_rank == 50


def test_write_summary_csv(runtmp):
    """test summary csv write function"""

    sum_gather = {'superkingdom': [SumGathInf(query_name='queryA', rank='superkingdom', fraction=1.0,
                                                          query_md5='queryA_md5', query_filename='queryA.sig',
                                                          f_weighted_at_rank=1.0, bp_match_at_rank=100,
                                                          lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),),
                                                          query_ani_at_rank=None,
                                                          total_weighted_hashes=0)],
                  'phylum':  [SumGathInf(query_name='queryA', rank='phylum', fraction=1.0,
                                                     query_md5='queryA_md5', query_filename='queryA.sig',
                                                     f_weighted_at_rank=1.0, bp_match_at_rank=100,
                                                     lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),
                                                              lca_utils.LineagePair(rank='phylum', name='b')),
                                                     query_ani_at_rank=None,
                                                     total_weighted_hashes=0)]}

    outs= runtmp.output("outsum.csv")
    with open(outs, 'w') as out_fp:
        write_summary(sum_gather, out_fp)

    sr = [x.rstrip().split(',') for x in open(outs, 'r')]
    print("gather_summary_results_from_file: \n", sr)
    assert ['query_name', 'rank', 'fraction', 'lineage', 'query_md5', 'query_filename', 'f_weighted_at_rank', 'bp_match_at_rank', 'query_ani_at_rank', 'total_weighted_hashes'] == sr[0]
    assert ['queryA', 'superkingdom', '1.0', 'a', 'queryA_md5', 'queryA.sig', '1.0', '100', '', '0'] == sr[1]
    assert ['queryA', 'phylum', '1.0', 'a;b', 'queryA_md5', 'queryA.sig', '1.0', '100','','0'] == sr[2]


def test_write_classification(runtmp):
    """test classification csv write function"""
    classif = ClassInf('queryA', 'match', 'phylum', 1.0,
                                    (lca_utils.LineagePair(rank='superkingdom', name='a'),
                                     lca_utils.LineagePair(rank='phylum', name='b')),
                                     'queryA_md5', 'queryA.sig', 1.0, 100,
                                     query_ani_at_rank=None)

    classification = {'phylum': [classif]}

    outs= runtmp.output("outsum.csv")
    with open(outs, 'w') as out_fp:
        write_classifications(classification, out_fp)

    sr = [x.rstrip().split(',') for x in open(outs, 'r')]
    print("gather_classification_results_from_file: \n", sr)
    assert ['query_name', 'status', 'rank', 'fraction', 'lineage', 'query_md5', 'query_filename', 'f_weighted_at_rank', 'bp_match_at_rank', 'query_ani_at_rank'] == sr[0]
    assert ['queryA', 'match', 'phylum', '1.0', 'a;b', 'queryA_md5', 'queryA.sig', '1.0', '100', ''] == sr[1]


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
        make_krona_header("strain")
    assert "Rank strain not present in available ranks" in str(exc.value)


def test_aggregate_by_lineage_at_rank_by_query():
    """test two queries, aggregate lineage at rank for each"""
    # make gather results
    gA = ["queryA","gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '100', '100']
    gB = ["queryA","gB","0.3","0.4", "queryA_md5", "queryA.sig", '0.5', '60', '140']
    gC = ["queryB","gB","0.3","0.3", "queryB_md5", "queryB.sig", '0.5', '60', '140']
    g_res = make_mini_gather_results([gA,gB,gC])

    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # aggregate by lineage at rank
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)
    print("superkingdom summarized gather results:", sk_sum)
    assert len(sk_sum) ==4
    assert sk_sum[0].query_name == "queryA"
    assert sk_sum[0].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[0].fraction == 0.9
    assert sk_sum[0].bp_match_at_rank == 160
    # check for unassigned for queryA
    assert sk_sum[1].query_name == "queryA"
    assert sk_sum[1].lineage == ()
    assert sk_sum[1].bp_match_at_rank == 40
    assert round(sk_sum[1].fraction,1) == 0.1
    # queryB
    assert sk_sum[2].query_name == "queryB"
    assert sk_sum[2].lineage == (lca_utils.LineagePair(rank='superkingdom', name='a'),)
    assert sk_sum[2].fraction == 0.3
    assert sk_sum[2].bp_match_at_rank == 60
    # check for unassigned for queryA
    assert sk_sum[3].query_name == "queryB"
    assert sk_sum[3].lineage == ()
    assert sk_sum[3].fraction == 0.7
    assert sk_sum[3].bp_match_at_rank == 140
    sk_lin_sum, query_names, num_queries = aggregate_by_lineage_at_rank(sk_sum, by_query=True)
    print("superkingdom lineage summary:", sk_lin_sum, '\n')
    assert sk_lin_sum == {(lca_utils.LineagePair(rank='superkingdom', name='a'),): {'queryA': 0.9, 'queryB': 0.3},
                          (): {'queryA': 0.09999999999999998, 'queryB': 0.7}}
    assert num_queries == 2
    assert query_names == ['queryA', 'queryB']

    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res)
    print("phylum summary:", phy_sum, ']\n')
    phy_lin_sum, query_names, num_queries = aggregate_by_lineage_at_rank(phy_sum, by_query=True)
    print("phylum lineage summary:", phy_lin_sum, '\n')
    assert phy_lin_sum ==  {(lca_utils.LineagePair(rank='superkingdom', name='a'), lca_utils.LineagePair(rank='phylum', name='b')): {'queryA': 0.5},
                            (lca_utils.LineagePair(rank='superkingdom', name='a'), lca_utils.LineagePair(rank='phylum', name='c')): {'queryA': 0.4, 'queryB': 0.3},
                            (): {'queryA': 0.09999999999999998, 'queryB': 0.7}}
    assert num_queries == 2
    assert query_names == ['queryA', 'queryB']


def test_format_for_krona_0():
    """test format for krona, equal matches"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # check krona format and check results!
    sk_sum, _, _ = summarize_gather_at("superkingdom", taxD, g_res)
    print("superkingdom summarized gather results:", sk_sum)
    krona_res = format_for_krona("superkingdom", {"superkingdom": sk_sum})
    print("krona_res: ", krona_res)
    assert krona_res == [(1.0, 'a')]

    phy_sum, _, _ = summarize_gather_at("phylum", taxD, g_res)
    krona_res = format_for_krona("phylum", {"phylum": phy_sum})
    print("krona_res: ", krona_res)
    assert krona_res == [(1.0, 'a', 'b')]


def test_format_for_krona_1():
    """test format for krona at each rank"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # summarize with all ranks
    sum_res = {}
    #for rank in lca_utils.taxlist(include_strain=False):
    for rank in ['superkingdom', 'phylum', 'class']:
        sum_res[rank], _, _ = summarize_gather_at(rank, taxD, g_res)
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
    """test two matches, equal f_unique_to_query"""
    # make gather results
    gA = ["queryA", "gA","0.5","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    gB = ["queryA", "gB","0.3","0.5", "queryA_md5", "queryA.sig", '0.5', '50', '50']
    g_res = make_mini_gather_results([gA,gB])

    # make mini taxonomy
    gA_tax = ("gA", "a;b;c")
    gB_tax = ("gB", "a;b;d")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    # summarize with all ranks
    sum_res = {}
    #for rank in lca_utils.taxlist(include_strain=False):
    for rank in ['superkingdom', 'phylum', 'class']:
        sum_res[rank], _, _ = summarize_gather_at(rank, taxD, g_res, best_only=True)
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
    """test two matches, equal f_unique_to_query"""
    class_krona_results =  [(0.5, 'a', 'b', 'c'), (0.5, 'a', 'b', 'd')]
    outk= runtmp.output("outkrona.tsv")
    with open(outk, 'w') as out_fp:
        write_krona("class", class_krona_results, out_fp)

    kr = [x.strip().split('\t') for x in open(outk, 'r')]
    print("krona_results_from_file: \n", kr)
    assert kr[0] == ["fraction", "superkingdom", "phylum", "class"]
    assert kr[1] == ["0.5", "a", "b", "c"]
    assert kr[2] == ["0.5", "a", "b", "d"]


def test_combine_sumgather_csvs_by_lineage(runtmp):
    # some summarized gather dicts
    sum_gather1 = {'superkingdom': [SumGathInf(query_name='queryA', rank='superkingdom', fraction=0.5,
                                                          query_md5='queryA_md5', query_filename='queryA.sig',
                                                          f_weighted_at_rank=1.0, bp_match_at_rank=100,
                                                          lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),),
                                                           query_ani_at_rank=None,
                                                           total_weighted_hashes=0)],
                  'phylum':  [SumGathInf(query_name='queryA', rank='phylum', fraction=0.5,
                                                     query_md5='queryA_md5', query_filename='queryA.sig',
                                                     f_weighted_at_rank=0.5, bp_match_at_rank=50,
                                                     lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),
                                                              lca_utils.LineagePair(rank='phylum', name='b')),
                                                     query_ani_at_rank=None,
                                                     total_weighted_hashes=0)]}
    sum_gather2 = {'superkingdom': [SumGathInf(query_name='queryB', rank='superkingdom', fraction=0.7,
                                                          query_md5='queryB_md5', query_filename='queryB.sig',
                                                          f_weighted_at_rank=0.7, bp_match_at_rank=70,
                                                          lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),),
                                                           query_ani_at_rank=None,
                                                           total_weighted_hashes=0)],
                  'phylum':  [SumGathInf(query_name='queryB', rank='phylum', fraction=0.7,
                                                     query_md5='queryB_md5', query_filename='queryB.sig',
                                                     f_weighted_at_rank=0.7, bp_match_at_rank=70,
                                                     lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),
                                                              lca_utils.LineagePair(rank='phylum', name='c')),
                                                     query_ani_at_rank=None,
                                                     total_weighted_hashes=0)]}

    # write summarized gather results csvs
    sg1= runtmp.output("sample1.csv")
    with open(sg1, 'w') as out_fp:
        write_summary(sum_gather1, out_fp)

    sg2= runtmp.output("sample2.csv")
    with open(sg2, 'w') as out_fp:
        write_summary(sum_gather2, out_fp)

    # test combine_summarized_gather_csvs_by_lineage_at_rank
    linD, query_names = combine_sumgather_csvs_by_lineage([sg1,sg2], rank="phylum")
    print("lineage_dict", linD)
    assert linD == {'a;b': {'queryA': '0.5'}, 'a;c': {'queryB': '0.7'}}
    assert query_names == ['queryA', 'queryB']
    linD, query_names = combine_sumgather_csvs_by_lineage([sg1,sg2], rank="superkingdom")
    print("lineage dict: \n", linD)
    assert linD, query_names == {'a': {'queryA': '0.5', 'queryB': '0.7'}}
    assert query_names == ['queryA', 'queryB']


def test_write_lineage_sample_frac(runtmp):
    outfrac = runtmp.output('outfrac.csv')
    sample_names = ['sample1', 'sample2']
    sk_linD = {'a': {'sample1': '0.500' ,'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, sk_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a', '0.500', '0.700']]

    phy_linD = {'a;b': {'sample1': '0.500'}, 'a;c': {'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, phy_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a;b', '0.500', '0'],  ['a;c', '0', '0.700']]


def test_write_lineage_sample_frac_format_lineage(runtmp):
    outfrac = runtmp.output('outfrac.csv')
    sample_names = ['sample1', 'sample2']
    sk_lineage = lca_utils.make_lineage('a')
    print(sk_lineage)
    sk_linD = {sk_lineage: {'sample1': '0.500' ,'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, sk_linD, out_fp, format_lineage=True)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a', '0.500', '0.700']]

    phy_lineage = lca_utils.make_lineage('a;b')
    print(phy_lineage)
    phy2_lineage = lca_utils.make_lineage('a;c')
    print(phy2_lineage)
    phy_linD = {phy_lineage: {'sample1': '0.500'}, phy2_lineage: {'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, phy_linD, out_fp, format_lineage=True)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a;b', '0.500', '0'],  ['a;c', '0', '0.700']]


def test_combine_sumgather_csvs_by_lineage_improper_rank(runtmp):
    # some summarized gather dicts
    sum_gather1 = {'superkingdom': [SumGathInf(query_name='queryA', rank='superkingdom', fraction=0.5,
                                                          query_md5='queryA_md5', query_filename='queryA.sig',
                                                          f_weighted_at_rank=0.5, bp_match_at_rank=50,
                                                          lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),),
                                                           query_ani_at_rank=None,
                                                           total_weighted_hashes=0)],
                  'phylum':  [SumGathInf(query_name='queryA', rank='phylum', fraction=0.5,
                                                     query_md5='queryA_md5', query_filename='queryA.sig',
                                                     f_weighted_at_rank=0.5, bp_match_at_rank=50,
                                                     lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),
                                                              lca_utils.LineagePair(rank='phylum', name='b')),
                                                     query_ani_at_rank=None,
                                                     total_weighted_hashes=0)]}
    sum_gather2 = {'superkingdom': [SumGathInf(query_name='queryB', rank='superkingdom', fraction=0.7,
                                                          query_md5='queryB_md5', query_filename='queryB.sig',
                                                          f_weighted_at_rank=0.7, bp_match_at_rank=70,
                                                          lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),),
                                                           query_ani_at_rank=None,
                                                           total_weighted_hashes=0)],
                  'phylum':  [SumGathInf(query_name='queryB', rank='phylum', fraction=0.7,
                                                     query_md5='queryB_md5', query_filename='queryB.sig',
                                                     f_weighted_at_rank=0.7, bp_match_at_rank=70,
                                                     lineage=(lca_utils.LineagePair(rank='superkingdom', name='a'),
                                                              lca_utils.LineagePair(rank='phylum', name='c')),
                                                     query_ani_at_rank=None,
                                                     total_weighted_hashes=0)]}

    # write summarized gather results csvs
    sg1= runtmp.output("sample1.csv")
    with open(sg1, 'w') as out_fp:
        write_summary(sum_gather1, out_fp)

    sg2= runtmp.output("sample2.csv")
    with open(sg2, 'w') as out_fp:
        write_summary(sum_gather2, out_fp)

    # test combine_summarized_gather_csvs_by_lineage_at_rank
    with pytest.raises(ValueError) as exc:
        linD, sample_names = combine_sumgather_csvs_by_lineage([sg1,sg2], rank="strain")
        print("ValueError: ", exc.value)
    assert "Rank strain not available." in str(exc.value)


def test_tax_multi_load_files(runtmp):
    # test loading various good and bad files
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    taxonomy_csv2 = utils.get_test_data('tax/test-strain.taxonomy.csv')
    badcsv = utils.get_test_data('tax/47+63_x_gtdb-rs202.gather.csv')

    db = MultiLineageDB.load([taxonomy_csv])
    assert len(db) == 6
    assert 'strain' not in db.available_ranks

    db = MultiLineageDB.load([taxonomy_csv2])
    assert len(db) == 6
    assert 'strain' in db.available_ranks
    assert db['GCF_001881345.1'][0].rank == 'superkingdom'

    # load a string rather than a list
    with pytest.raises(TypeError):
        MultiLineageDB.load(badcsv)

    # load a bad CSV
    with pytest.raises(ValueError):
        MultiLineageDB.load([badcsv])

    # load a directory
    with pytest.raises(ValueError):
        MultiLineageDB.load([runtmp.output('')])

    # file does not exist
    with pytest.raises(ValueError):
        MultiLineageDB.load([runtmp.output('no-such-file')])


def test_tax_sql_load_new_file(runtmp):
    # test loading a newer-format sql file with sourmash_internals table
    taxonomy_db = utils.get_test_data('sqlite/test.taxonomy.db')

    db = MultiLineageDB.load([taxonomy_db])
    print(list(db.keys()))
    assert len(db) == 6
    assert 'strain' not in db.available_ranks
    assert db['GCF_001881345'][0].rank == 'superkingdom'


def test_tax_multi_load_files_shadowed(runtmp):
    # test loading various good and bad files
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    taxonomy_csv2 = utils.get_test_data('tax/test-strain.taxonomy.csv')
    taxonomy_db = utils.get_test_data('tax/test.taxonomy.db')

    db = MultiLineageDB.load([taxonomy_csv, taxonomy_csv2, taxonomy_db],
                             keep_full_identifiers=False,
                             keep_identifier_versions=False)
    assert len(db.shadowed_identifiers()) == 6

    # we should have everything including strain
    assert set(lca_utils.taxlist()) == set(db.available_ranks)

    db = MultiLineageDB.load([taxonomy_csv, taxonomy_db],
                             keep_full_identifiers=False,
                             keep_identifier_versions=False)
    assert len(db.shadowed_identifiers()) == 6
    assert set(lca_utils.taxlist(include_strain=False)) == set(db.available_ranks)


def test_tax_multi_save_files(runtmp, keep_identifiers, keep_versions):
    # test save
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    if keep_identifiers and not keep_versions:
        with pytest.raises(ValueError):
            db = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=keep_identifiers,
                                     keep_identifier_versions=keep_versions)
        return

    db = MultiLineageDB.load([taxonomy_csv],
                             keep_full_identifiers=keep_identifiers,
                             keep_identifier_versions=keep_versions)

    out_db = runtmp.output('out.db')
    out_csv = runtmp.output('out.csv')
    out2_csv = runtmp.output('out2.csv')

    # can't save to fp with sql
    with open(out_csv, 'wt') as fp:
        with pytest.raises(ValueError):
            db.save(fp, 'sql')

    # these should all work...
    with open(out_csv, 'wt') as fp:
        db.save(fp, 'csv')

    db.save(out2_csv, 'csv')
    db.save(out_db, 'sql')

    # ...and be equal
    db1 = db.load([out_db])
    db2 = db.load([out_csv])
    db3 = db.load([out2_csv])

    def strip_strain(it):
        for k, v in it:
            if v[-1].rank == 'strain':
                v = v[:-1]
            yield k, v

    import pprint
    db_items = list(strip_strain(db.items()))
    db1_items = list(strip_strain(db1.items()))
    db2_items = list(strip_strain(db2.items()))
    db3_items = list(strip_strain(db3.items()))
    pprint.pprint(db_items)
    print('XXX')
    pprint.pprint(list(db1_items))
    print('XXX')
    pprint.pprint(list(db2_items))

    assert set(db_items) == set(db1_items)
    assert set(db_items) == set(db2_items)
    assert set(db_items) == set(db3_items)


def test_lineage_db_csv_load(runtmp):
    # test LineageDB.load
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    taxonomy_csv2 = utils.get_test_data('tax/test-strain.taxonomy.csv')
    badcsv = utils.get_test_data('tax/47+63_x_gtdb-rs202.gather.csv')
    badcsv2 = utils.get_test_data('tax/test-missing-ranks.taxonomy.csv')

    db = LineageDB.load(taxonomy_csv)
    assert len(db) == 6
    assert 'strain' not in db.available_ranks

    db = LineageDB.load(taxonomy_csv2)
    assert len(db) == 6
    assert 'strain' in db.available_ranks

    # load the wrong kind of csv
    with pytest.raises(ValueError):
        LineageDB.load(badcsv)

    # load a bad CSV
    with pytest.raises(ValueError):
        LineageDB.load(badcsv2)

    # load a directory
    with pytest.raises(ValueError):
        LineageDB.load(runtmp.output(''))

    # file does not exist
    with pytest.raises(ValueError):
        LineageDB.load(runtmp.output('no-such-file'))

    # construct a CSV with bad headers
    with open(runtmp.output('xxx.csv'), 'w', newline="") as fp:
        fp.write('x,y,z\n')
    with pytest.raises(ValueError):
        LineageDB.load(runtmp.output('xxx.csv'))


def test_lineage_db_sql_load(runtmp):
    # test LineageDB_sqlite.load
    taxonomy_db = utils.get_test_data('tax/test.taxonomy.db')
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')

    db = LineageDB_Sqlite.load(taxonomy_db)
    assert bool(db)
    assert len(db) == 6
    db.available_ranks
    assert 'strain' not in db.available_ranks
    assert db['GCF_001881345'][0].rank == 'superkingdom'
    with pytest.raises(KeyError):
        db['foo']

    # load any kind of CSV
    with pytest.raises(ValueError):
        LineageDB_Sqlite.load(taxonomy_csv)

    # load a directory
    with pytest.raises(ValueError):
        LineageDB_Sqlite.load(runtmp.output(''))

    # file does not exist
    with pytest.raises(ValueError):
        LineageDB_Sqlite.load(runtmp.output('no-such-file'))


def test_LineagePair():
    lin = LineagePair(rank="rank1", name='name1')
    print(lin)


def test_LineagePair_1():
    lin1 = LineagePair(rank="rank1", name='name1', taxid=1)
    print(lin1)


def test_BaseLineageInfo_init_empty():
    ranks=["A", "B", "C"]
    taxinf = BaseLineageInfo(ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['', '', ''] # this is a bit odd, but it's what preserves empty ranks...
    print(taxinf.filled_lineage)
    assert taxinf.filled_lineage == ()
    assert taxinf.lowest_lineage_name == ""
    assert taxinf.lowest_lineage_taxid == ""
    assert taxinf.filled_ranks == ()
    assert taxinf.lowest_rank == None
    assert taxinf.display_lineage() == ""
    assert taxinf.display_lineage(null_as_unclassified=True) == "unclassified"


def test_BaseLineageInfo_init_lineage_str():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    taxinf = BaseLineageInfo(lineage_str=x, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c']
    print(taxinf.filled_lineage)
    assert taxinf.filled_lineage == (LineagePair(rank='A', name='a', taxid=None),
                                     LineagePair(rank='B', name='b', taxid=None),
                                     LineagePair(rank='C', name='c', taxid=None))
    assert taxinf.lowest_lineage_name == "c"
    assert taxinf.lowest_rank == "C"

def test_BaseLineageInfo_init_lineage_str_comma_sep():
    x = "a,b,c"
    ranks=["A", "B", "C"]
    taxinf = BaseLineageInfo(lineage_str=x, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c']
    print(taxinf.filled_lineage)
    assert taxinf.lowest_lineage_name == "c"


def test_BaseLineageInfo_init_lineage_tups():
    ranks=["A", "B", "C"]
    lin_tups = (LineagePair(rank="A", name='a'), LineagePair(rank="C", name='b'))
    taxinf = BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', '', 'b']


def test_BaseLineageInfo_init_lca_lineage_tups():
    ranks=["A", "B", "C"]
    lin_tups = (lca_utils.LineagePair(rank="A", name='a'), lca_utils.LineagePair(rank="C", name='b'))
    taxinf = BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', '', 'b']


def test_BaseLineageInfo_init_lineage_dict_fail():
    ranks=["A", "B", "C"]
    lin_tups = (lca_utils.LineagePair(rank="A", name='a'), lca_utils.LineagePair(rank="C", name='b'))
    with pytest.raises(ValueError) as exc:
        taxinf = BaseLineageInfo(ranks=ranks, lineage_dict=lin_tups)
    print(str(exc))

    assert "is not dictionary" in str(exc)


def test_BaseLineageInfo_init_lineage_dict  ():
    x = {'rank1': 'name1', 'rank2': 'name2'}
    taxinf = BaseLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', 'name2']


def test_BaseLineageInfo_init_lineage_dict_withtaxid():
    x = {'rank1': {'name': 'name1', 'taxid': 1}, 'rank2': {'name':'name2', 'taxid': 2}}
    taxinf = BaseLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', 'name2']
    assert taxinf.zip_taxid()== ['1', '2']
    assert taxinf.lowest_lineage_taxid == 2
    assert taxinf.lowest_lineage_name == "name2"


def test_BaseLineageInfo_init_lineage_str_lineage_dict_test_eq():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    rankD = {"A": "a", "B": "b", "C": "c"}
    lin1 = BaseLineageInfo(lineage_str=x, ranks=ranks)
    lin2 = BaseLineageInfo(lineage_dict=rankD, ranks=ranks)
    assert lin1 == lin2


def test_BaseLineageInfo_init_no_ranks():
    x = "a;b;c"
    rankD = {"superkingdom": "a", "phylum": "b", "class": "c"}
    lin_tups = (LineagePair(rank="rank2", name='name1'), LineagePair(rank="rank1", name='name1'))
    with pytest.raises(TypeError) as exc:
        BaseLineageInfo(lineage_str=x)
    print(exc)
    assert "__init__() missing 1 required positional argument: 'ranks'" in str(exc)
    with pytest.raises(TypeError) as exc:
        BaseLineageInfo(lineage_dict=rankD)
    print(exc)
    assert "__init__() missing 1 required positional argument: 'ranks'" in str(exc)
    with pytest.raises(TypeError) as exc:
        BaseLineageInfo(lineage=lin_tups)
    print(exc)
    assert "__init__() missing 1 required positional argument: 'ranks'" in str(exc)


def test_BaseLineageInfo_init_with_wrong_ranks():
    ranks=["A", "B", "C"]
    lin_tups = [LineagePair(rank="rank1", name='name1')]
    linD = {"rank1": "a"}
    with pytest.raises(ValueError) as exc:
        BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(str(exc))
    assert "Rank 'rank1' not present in A, B, C" in str(exc)
    with pytest.raises(ValueError) as exc:
        BaseLineageInfo(lineage_dict=linD, ranks=ranks)
    print(str(exc))
    assert "Rank 'rank1' not present in A, B, C" in str(exc)


def test_BaseLineageInfo_init_not_lineagepair():
    ranks=["A", "B", "C"]
    lin_tups = (("rank1", "name1"),)
    with pytest.raises(ValueError) as exc:
        BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(str(exc))
    assert "is not LineagePair" in str(exc)


def test_RankLineageInfo_taxlist():
    taxinf = RankLineageInfo()
    taxranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
    assert taxinf.taxlist == taxranks
    assert taxinf.ascending_taxlist == taxranks[::-1]


def test_RankLineageInfo_init_lineage_str():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c', '', '', '', '', '']

def test_RankLineageInfo_init_lineage_str_with_ranks_as_list():
    x = "a;b;c"
    taxranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxinf = RankLineageInfo(lineage_str=x, ranks=taxranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', 'c', '', '', '', '']


def test_RankLineageInfo_init_lineage_tups():
    x = (LineagePair(rank="superkingdom", name='a'), LineagePair(rank="phylum", name='b'))
    taxinf = RankLineageInfo(lineage=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict():
    x = {"superkingdom":'a',"phylum":'b'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_missing_rank():
    x = {'superkingdom': 'name1', 'class': 'name2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']
    assert taxinf.zip_lineage(truncate_empty=True)== ['name1', '', 'name2']


def test_RankLineageInfo_init_lineage_dict_missing_rank_withtaxid():
    x = {'superkingdom': {'name': 'name1', 'taxid': 1}, 'class': {'name':'name2', 'taxid': 2}}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']
    assert taxinf.zip_taxid()== ['1', '', '2', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_str_lineage_dict_test_eq():
    x = "a;b;c"
    rankD = {"superkingdom": "a", "phylum": "b", "class": "c"}
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_dict=rankD)
    print("lin1: ", lin1)
    print("lin2: ", lin2)
    assert lin1 == lin2

def test_RankLineageInfo_init_lineage_str_1_truncate():
    x = "a;b;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', 'c']


def test_RankLineageInfo_init_lineage_str_2():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', 'c' '', '', '', '', '']


def test_RankLineageInfo_init_lineage_str_2_truncate():
    x = "a;b;;c"
    taxinf = RankLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage(truncate_empty=True)== ['a', 'b', '', 'c']


def test_RankLineageInfo_init_lineage_with_incorrect_rank():
    x = [ LineagePair('superkingdom', 'a'), LineagePair("NotARank", ''), LineagePair('class', 'c') ]
    with pytest.raises(ValueError) as exc:
        RankLineageInfo(lineage=x)
    print(str(exc))
    assert f"Rank 'NotARank' not present in " in str(exc)


def test_zip_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage() == ['a', 'b', '', '', '', '', '', '']


def test_zip_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    print("ranks: ", taxinf.ranks)
    print("zipped lineage: ", taxinf.zip_lineage(truncate_empty=True))
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', 'b']


def test_zip_lineage_3():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage() == ['a', '', 'c', '', '', '', '', '']


def test_zip_lineage_3_truncate():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_zip_lineage_4():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.zip_lineage(truncate_empty=True) == ['a', '', 'c']


def test_display_lineage_1():
    x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.display_lineage() == "a;b"


def test_display_lineage_2():
    x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    taxinf = RankLineageInfo(lineage=x)
    assert taxinf.display_lineage() == "a;;c"


def test_display_taxid_1():
    x = [ LineagePair('superkingdom', 'a', 1), LineagePair('phylum', 'b', 2) ]
    taxinf = RankLineageInfo(lineage=x)
    print(taxinf)
    assert taxinf.display_taxid() == "1;2"

def test_display_taxid_2():
    x = [ LineagePair('superkingdom', 'name1', 1), LineagePair(None, ''), LineagePair    ('class', 'name2',2) ]
    taxinf = RankLineageInfo(lineage=x)
    print(taxinf)
    assert taxinf.display_taxid() == "1;;2"


def test_is_lineage_match_1():
    # basic behavior: match at order and above, but not at family or below.
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    assert lin1.is_lineage_match(lin2, 'superkingdom')
    assert lin2.is_lineage_match(lin1, 'superkingdom')
    assert lin1.is_lineage_match(lin2, 'phylum')
    assert lin2.is_lineage_match(lin1, 'phylum')
    assert lin1.is_lineage_match(lin2, 'class')
    assert lin2.is_lineage_match(lin1, 'class')
    assert lin1.is_lineage_match(lin2, 'order')
    assert lin2.is_lineage_match(lin1, 'order')

    assert not lin1.is_lineage_match(lin2, 'family')
    assert not lin2.is_lineage_match(lin1, 'family')
    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_2():
    # match at family, and above, levels; no genus or species to match
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    assert lin1.is_lineage_match(lin2, 'superkingdom')
    assert lin2.is_lineage_match(lin1, 'superkingdom')
    assert lin1.is_lineage_match(lin2, 'phylum')
    assert lin2.is_lineage_match(lin1, 'phylum')
    assert lin1.is_lineage_match(lin2, 'class')
    assert lin2.is_lineage_match(lin1, 'class')
    assert lin1.is_lineage_match(lin2, 'order')
    assert lin2.is_lineage_match(lin1, 'order')
    assert lin1.is_lineage_match(lin2, 'family')
    assert lin2.is_lineage_match(lin1, 'family')

    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_3():
    # one lineage is empty
    lin1 = RankLineageInfo()
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')

    assert not lin1.is_lineage_match(lin2, 'superkingdom')
    assert not lin2.is_lineage_match(lin1, 'superkingdom')
    assert not lin1.is_lineage_match(lin2, 'phylum')
    assert not lin2.is_lineage_match(lin1, 'phylum')
    assert not lin1.is_lineage_match(lin2, 'class')
    assert not lin2.is_lineage_match(lin1, 'class')
    assert not lin1.is_lineage_match(lin2, 'order')
    assert not lin2.is_lineage_match(lin1, 'order')
    assert not lin1.is_lineage_match(lin2, 'family')
    assert not lin2.is_lineage_match(lin1, 'family')
    assert not lin1.is_lineage_match(lin2, 'genus')
    assert not lin2.is_lineage_match(lin1, 'genus')
    assert not lin1.is_lineage_match(lin2, 'species')
    assert not lin2.is_lineage_match(lin1, 'species')


def test_is_lineage_match_incorrect_ranks():
    #test comparison with incompatible ranks
    taxranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e', ranks=taxranks[::-1])
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    with pytest.raises(ValueError) as exc:
        lin1.is_lineage_match(lin2, 'superkingdom')
    print(str(exc))
    assert 'Cannot compare lineages from taxonomies with different ranks.' in str(exc)


def test_is_lineage_match_improper_rank():
    #test comparison with incompatible ranks
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    with pytest.raises(ValueError) as exc:
        lin1.is_lineage_match(lin2, 'NotARank')
    print(str(exc))
    assert "Desired Rank 'NotARank' not available for this lineage" in str(exc)


def test_pop_to_rank_1():
    # basic behavior - pop to order?
    lin1 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d')
    lin2 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d;f__f')

    print(lin1)
    popped = lin2.pop_to_rank('order')
    print(popped)
    assert popped == lin1


def test_pop_to_rank_2():
    # what if we're already above rank?
    lin2 = RankLineageInfo(lineage_str='d__a;p__b;c__c;o__d;f__f')
    print(lin2.pop_to_rank('species'))
    assert lin2.pop_to_rank('species') == lin2


def test_pop_to_rank_rank_not_avail():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(ValueError) as exc:
        lin1.pop_to_rank("NotARank")
    print(str(exc))
    assert "Desired Rank 'NotARank' not available for this lineage" in str(exc)


def test_lineage_at_rank_norank():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(TypeError) as exc:
        lin1.lineage_at_rank()
    print(str(exc))
    assert "lineage_at_rank() missing 1 required positional argument: 'rank'" in str(exc)


def test_lineage_at_rank_rank_not_avail():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    with pytest.raises(ValueError) as exc:
        lin1.lineage_at_rank("NotARank")
    print(str(exc))
    assert "Desired Rank 'NotARank' not available for this lineage" in str(exc)


def test_lineage_at_rank_1():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage_at_rank('superkingdom'))
    
    assert lin1.lineage_at_rank('superkingdom') == (LineagePair(rank='superkingdom', name='d__a', taxid=None),)
    print(lin1.lineage_at_rank('class'))
    assert lin1.lineage_at_rank('class') == (LineagePair(rank='superkingdom', name='d__a', taxid=None),
                                             LineagePair(rank='phylum', name='p__b', taxid=None),
                                             LineagePair(rank='class', name='c__c', taxid=None))


def test_lineage_at_rank_below_rank():
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage_at_rank('superkingdom'))
    # if rank is not provided, we only return the filled lineage, to follow original pop_to_rank behavior.

    print(lin1.lineage_at_rank('genus'))
    assert lin1.lineage_at_rank('genus') == (LineagePair(rank='superkingdom', name='d__a', taxid=None),
                                             LineagePair(rank='phylum', name='p__b', taxid=None),
                                             LineagePair(rank='class', name='c__c', taxid=None),
                                             LineagePair(rank='order', name='o__d', taxid=None),
                                             LineagePair(rank='family', name='f__f', taxid=None))
