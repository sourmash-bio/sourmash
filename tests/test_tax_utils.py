"""
Tests for functions in taxonomy submodule.
"""

import pytest
from pytest import approx
import os
from os.path import basename
import gzip

import sourmash_tst_utils as utils

from sourmash.tax.tax_utils import (ascending_taxlist, get_ident, load_gather_results,
                                    collect_gather_csvs, check_and_load_gather_csvs,
                                    LineagePair, QueryInfo, GatherRow, TaxResult, QueryTaxResult,
                                    SummarizedGatherResult, ClassificationResult, AnnotateTaxResult,
                                    BaseLineageInfo, RankLineageInfo, LINLineageInfo,
                                    aggregate_by_lineage_at_rank, format_for_krona,
                                    write_krona, write_lineage_sample_frac, read_lingroups,
                                    LineageTree, LineageDB, LineageDB_Sqlite, MultiLineageDB)

# utility functions for testing
def make_mini_taxonomy(tax_info, LIN=False):
    #pass in list of tuples: (name, lineage)
    taxD = {}
    for (name, lin) in tax_info:
        if LIN:
            lineage = LINLineageInfo(lineage_str=lin)
        else:
            lineage = RankLineageInfo(lineage_str=lin)
        taxD[name] = lineage.filled_lineage
    return taxD

def make_mini_taxonomy_with_taxids(tax_info, LIN=False):
    taxD = {}
    for (name, lin, taxids) in tax_info:
        if LIN:
            lineage = LINLineageInfo(lineage_str=lin)
        else:
            ranks = RankLineageInfo.ranks
            txs = taxids.split(';')
            lns = lin.split(';')
            lineage_tups = []
            for n, taxname in enumerate(lns):
                rk = ranks[n]
                tx = txs[n]
                this_lineage = LineagePair(rk, name=taxname, taxid=tx)
                lineage_tups.append(this_lineage)
            lineage = RankLineageInfo(lineage=lineage_tups)
        taxD[name] = lineage.filled_lineage
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


def make_TaxResult(gather_dict=None, taxD=None, keep_full_ident=False, keep_ident_version=False, skip_idents=None, LIN=False):
    """Make TaxResult from artificial gather row (dict)"""
    gRow = make_GatherRow(gather_dict)
    taxres = TaxResult(raw=gRow, keep_full_identifiers=keep_full_ident,
                       keep_identifier_versions=keep_ident_version, lins=LIN)
    if taxD is not None:
        taxres.get_match_lineage(tax_assignments=taxD, skip_idents=skip_idents)
    return taxres


def make_QueryTaxResults(gather_info, taxD=None, single_query=False, keep_full_ident=False, keep_ident_version=False,
                        skip_idents=None, summarize=False, classify=False, classify_rank=None, c_thresh=0.1, ani_thresh=None,
                        LIN=False):
    """Make QueryTaxResult(s) from artificial gather information, formatted as list of gather rows (dicts)"""
    gather_results = {}
    this_querytaxres = None
    for gather_infoD in gather_info:
        taxres = make_TaxResult(gather_infoD, taxD=taxD,  keep_full_ident=keep_full_ident,
                                keep_ident_version=keep_ident_version, skip_idents=skip_idents, LIN=LIN)
        query_name = taxres.query_name
        # add to matching QueryTaxResult or create new one
        if not this_querytaxres or not this_querytaxres.is_compatible(taxres):
            # get existing or initialize new
            this_querytaxres = gather_results.get(query_name, QueryTaxResult(taxres.query_info, lins=LIN))
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


def test_QueryInfo_basic():
    "basic functionality of QueryInfo dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    assert qInf.query_name == 'q1'
    assert isinstance(qInf.query_n_hashes, int)
    assert isinstance(qInf.ksize, int)
    assert isinstance(qInf.scaled, int)
    assert qInf.total_weighted_hashes == 200
    assert qInf.total_weighted_bp == 2000


def test_QueryInfo_no_hash_info():
    "QueryInfo dataclass for older gather results without query_n_hashes or total_weighted_hashes"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',ksize=31,scaled=10)
    assert qInf.query_name == 'q1'
    assert qInf.query_n_hashes == 0
    assert qInf.total_weighted_hashes == 0
    assert qInf.total_weighted_bp == 0


def test_QueryInfo_missing():
    "check that required args"
    with pytest.raises(TypeError) as exc:
        QueryInfo(query_name='q1', query_filename='f1',query_bp='100',query_n_hashes='10',ksize=31,scaled=10, total_weighted_hashes=200)
    print(str(exc))
    assert "missing 1 required positional argument: 'query_md5'" in str(exc)


def test_SummarizedGatherResult():
    "basic functionality of SummarizedGatherResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    print(sgr)
    assert sgr.rank=='phylum'
    sumD = sgr.as_summary_dict(query_info=qInf)
    print(sumD)
    assert sumD == {'rank': 'phylum', 'fraction': "0.2", 'lineage': 'a;b', 'f_weighted_at_rank': "0.3",
                    'bp_match_at_rank': "30", 'query_ani_at_rank': None, 'query_name': 'q1',
                    'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}
    hD = sgr.as_human_friendly_dict(query_info=qInf)
    print(hD)
    assert hD == {'rank': 'phylum', 'fraction': '0.200', 'lineage': 'a;b', 'f_weighted_at_rank': '30.0%',
                  'bp_match_at_rank': "30", 'query_ani_at_rank': '-    ', 'query_name': 'q1',
                  'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}
    krD = sgr.as_kreport_dict(query_info=qInf)
    print(krD)
    assert krD == {'ncbi_taxid': None, 'sci_name': 'b', 'rank_code': 'P', 'num_bp_assigned': "0",
                   'percent_containment': '30.00', 'num_bp_contained': "600"}
    lD = sgr.as_lineage_dict(ranks = RankLineageInfo().ranks, query_info=qInf)
    print(lD)
    assert lD == {'ident': 'q1', 'superkingdom': 'a', 'phylum': 'b', 'class': '', 'order': '',
                  'family': '', 'genus': '', 'species': '', 'strain': ''}
    cami = sgr.as_cami_bioboxes()
    print(cami)
    assert cami == [None, 'phylum', None, 'a|b', '30.00']


def test_SummarizedGatherResult_withtaxids():
    "basic functionality of SummarizedGatherResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    lin = [LineagePair(rank='superkingdom', name='a', taxid='1'), LineagePair(rank='phylum', name='b', taxid=2)]
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage=lin),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    print(sgr)
    assert sgr.rank=='phylum'
    sumD = sgr.as_summary_dict(query_info=qInf)
    print(sumD)
    assert sumD == {'rank': 'phylum', 'fraction': "0.2", 'lineage': 'a;b', 'f_weighted_at_rank': "0.3",
                    'bp_match_at_rank': "30", 'query_ani_at_rank': None, 'query_name': 'q1',
                    'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}
    hD = sgr.as_human_friendly_dict(query_info=qInf)
    print(hD)
    assert hD == {'rank': 'phylum', 'fraction': '0.200', 'lineage': 'a;b', 'f_weighted_at_rank': '30.0%',
                  'bp_match_at_rank': "30", 'query_ani_at_rank': '-    ', 'query_name': 'q1',
                  'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}
    krD = sgr.as_kreport_dict(query_info=qInf)
    print(krD)
    assert krD == {'ncbi_taxid': '2', 'sci_name': 'b', 'rank_code': 'P', 'num_bp_assigned': "0",
                   'percent_containment': '30.00', 'num_bp_contained': "600"}
    lD = sgr.as_lineage_dict(ranks = RankLineageInfo().ranks, query_info=qInf)
    print(lD)
    assert lD == {'ident': 'q1', 'superkingdom': 'a', 'phylum': 'b', 'class': '', 'order': '',
                  'family': '', 'genus': '', 'species': '', 'strain': ''}
    cami = sgr.as_cami_bioboxes()
    print(cami)
    assert cami == ['2', 'phylum', '1|2', 'a|b', '30.00']


def test_SummarizedGatherResult_LINs():
    "SummarizedGatherResult with LINs"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=LINLineageInfo(lineage_str="0;0;1"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)

    lgD = sgr.as_lingroup_dict(query_info=qInf, lg_name="lg_name")
    print(lgD)
    assert lgD == {'name': "lg_name", "lin": "0;0;1",
                   'percent_containment': '30.00', 'num_bp_contained': "600"}
    lgD = sgr.as_lingroup_dict(query_info=qInf, lg_name="lg_name")
    print(lgD)
    assert lgD == {'name': "lg_name", "lin": "0;0;1",
                   'percent_containment': '30.00', 'num_bp_contained': "600"}
    with pytest.raises(ValueError) as exc:
        sgr.as_kreport_dict(query_info=qInf)
    print(str(exc))
    assert "Cannot produce 'kreport' with LIN taxonomy." in str(exc)
    with pytest.raises(ValueError) as exc:
        sgr.as_cami_bioboxes()
    print(str(exc))
    assert "Cannot produce 'bioboxes' with LIN taxonomy." in str(exc)


def test_SummarizedGatherResult_set_query_ani():
    "Check ANI estimation within SummarizedGatherResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    sgr.set_query_ani(query_info=qInf)
    print(sgr.query_ani_at_rank)
    assert sgr.query_ani_at_rank == approx(0.949,  rel=1e-3)
    # ANI can be calculated with query_bp OR query_n_hashes. Remove each and check the results are identical
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes=0,ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    sgr.set_query_ani(query_info=qInf)
    print(sgr.query_ani_at_rank)
    assert sgr.query_ani_at_rank == approx(0.949,  rel=1e-3)
    # try without query_bp
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp=0,
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    sgr.set_query_ani(query_info=qInf)
    print(sgr.query_ani_at_rank)
    assert sgr.query_ani_at_rank == approx(0.949,  rel=1e-3)


def test_SummarizedGatherResult_greater_than_1():
    "basic functionality of SummarizedGatherResult dataclass"
    # fraction > 1
    with pytest.raises(ValueError) as exc:
        SummarizedGatherResult(rank="phylum", fraction=0.3, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=1.2, bp_match_at_rank=30)
    print(str(exc))
    assert "> 100% of the query!" in str(exc)
    # f_weighted > 1
    with pytest.raises(ValueError) as exc:
        SummarizedGatherResult(rank="phylum", fraction=1.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    print(str(exc))
    assert "> 100% of the query!" in str(exc)


def test_SummarizedGatherResult_0_fraction():
    with pytest.raises(ValueError) as exc:
        SummarizedGatherResult(rank="phylum", fraction=-.1, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    err_msg = "Summarized fraction is <=0% of the query! This should not occur."
    assert err_msg in str(exc)
    #assert cr.status == 'nomatch'
    
    with pytest.raises(ValueError) as exc:
        SummarizedGatherResult(rank="phylum", fraction=.1, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0, bp_match_at_rank=30)
    print(str(exc))
    assert err_msg in str(exc)


def test_SummarizedGatherResult_species_kreport():
    "basic functionality of SummarizedGatherResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="species", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b;c;d;e;f;g"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    print(sgr)
    assert sgr.rank=='species'
    krD = sgr.as_kreport_dict(query_info=qInf)
    print(krD)
    assert krD == {'ncbi_taxid': None, 'sci_name': 'g', 'rank_code': 'S', 'num_bp_assigned': "600",
                   'percent_containment': '30.00', 'num_bp_contained': "600"}


def test_SummarizedGatherResult_summary_dict_limit_float():
    "basic functionality of SummarizedGatherResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    sgr = SummarizedGatherResult(rank="phylum", fraction=0.123456, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.345678, bp_match_at_rank=30)
    print(sgr)
    assert sgr.rank=='phylum'
    sumD = sgr.as_summary_dict(query_info=qInf)
    print(sumD)
    assert sumD == {'rank': 'phylum', 'fraction': "0.123456", 'lineage': 'a;b', 'f_weighted_at_rank': "0.345678",
                    'bp_match_at_rank': "30", 'query_ani_at_rank': None, 'query_name': 'q1',
                    'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}
    
    sumD = sgr.as_summary_dict(query_info=qInf, limit_float=True)
    print(sumD)
    assert sumD == {'rank': 'phylum', 'fraction': "0.123", 'lineage': 'a;b', 'f_weighted_at_rank': "0.346",
                    'bp_match_at_rank': "30", 'query_ani_at_rank': None, 'query_name': 'q1',
                    'query_md5': 'md5', 'query_filename': 'f1', 'total_weighted_hashes': "200"}


def test_ClassificationResult():
    "basic functionality of ClassificationResult dataclass"
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    cr = ClassificationResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                              f_weighted_at_rank=0.3, bp_match_at_rank=30, query_ani_at_rank=0.97)
    cr.set_status(query_info=qInf, containment_threshold=0.1)
    assert cr.status == 'match'
    print(cr.query_ani_at_rank)
    assert cr.query_ani_at_rank == approx(0.949,  rel=1e-3)
    cr.set_status(query_info=qInf, containment_threshold=0.35)
    assert cr.status == 'below_threshold'
    lD = cr.as_lineage_dict(ranks = RankLineageInfo().ranks, query_info=qInf)
    print(lD)
    assert lD == {'ident': 'q1', 'superkingdom': 'a', 'phylum': 'b', 'class': '', 'order': '',
                  'family': '', 'genus': '', 'species': '', 'strain': ''}


def test_ClassificationResult_greater_than_1():
    "basic functionality of SummarizedGatherResult dataclass"
    # fraction > 1
    with pytest.raises(ValueError) as exc:
        ClassificationResult(rank="phylum", fraction=0.3, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=1.2, bp_match_at_rank=30)
    print(str(exc))
    assert "> 100% of the query!" in str(exc)
    # f_weighted > 1
    with pytest.raises(ValueError) as exc:
        ClassificationResult(rank="phylum", fraction=1.2, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    print(str(exc))
    assert "> 100% of the query!" in str(exc)


def test_ClassificationResult_0_fraction():
    with pytest.raises(ValueError) as exc:
        ClassificationResult(rank="phylum", fraction=-.1, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0.3, bp_match_at_rank=30)
    err_msg = "Summarized fraction is <=0% of the query! This should not occur."
    assert err_msg in str(exc)
    #assert cr.status == 'nomatch'
    
    with pytest.raises(ValueError) as exc:
        ClassificationResult(rank="phylum", fraction=.1, lineage=RankLineageInfo(lineage_str="a;b"),
                                 f_weighted_at_rank=0, bp_match_at_rank=30)
    print(str(exc))
    assert err_msg in str(exc)


def test_ClassificationResult_build_krona_result():
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    cr = ClassificationResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                              f_weighted_at_rank=0.3, bp_match_at_rank=30, query_ani_at_rank=0.97)
    #cr.set_status(query_info=qInf, rank='phylum')
    kr, ukr = cr.build_krona_result(rank='phylum')
    print(kr)
    assert kr == (0.2, 'a', 'b')
    print(ukr)
    assert ukr == (0.8, 'unclassified', 'unclassified')  


def test_ClassificationResult_build_krona_result_no_rank():
    qInf = QueryInfo(query_name='q1', query_md5='md5', query_filename='f1',query_bp='100',
                     query_n_hashes='10',ksize='31',scaled='10', total_weighted_hashes='200')
    cr = ClassificationResult(rank="phylum", fraction=0.2, lineage=RankLineageInfo(lineage_str="a;b"),
                              f_weighted_at_rank=0.3, bp_match_at_rank=30, query_ani_at_rank=0.97)
    cr.set_status(query_info=qInf, containment_threshold=0.1)


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


def test_AnnotateTaxResult_get_ident_default():
    gA = {"name": "GCF_001881345.1"}  # gather result with match name as GCF_001881345.1
    taxres = AnnotateTaxResult(raw=gA)
    print(taxres.match_ident)
    assert taxres.match_ident == "GCF_001881345"


def test_AnnotateTaxResult_get_ident_idcol():
    gA = {"name": "n1", "match_name": "n2", "ident": "n3", "accession": "n4"}  # gather result with match name as GCF_001881345.1
    taxres = AnnotateTaxResult(raw=gA)
    print(taxres.match_ident)
    assert taxres.match_ident == "n1"
    taxres = AnnotateTaxResult(raw=gA, id_col="match_name")
    print(taxres.match_ident)
    assert taxres.match_ident == "n2"
    taxres = AnnotateTaxResult(raw=gA, id_col="ident")
    print(taxres.match_ident)
    assert taxres.match_ident == "n3"
    taxres = AnnotateTaxResult(raw=gA, id_col="accession")
    print(taxres.match_ident)
    assert taxres.match_ident == "n4"


def test_AnnotateTaxResult_get_ident_idcol_fail():
    gA = {"name": "n1", "match_name": "n2", "ident": "n3", "accession": "n4"}  # gather result with match name as GCF_001881345.1
    with pytest.raises(ValueError) as exc:
        AnnotateTaxResult(raw=gA, id_col="NotACol")
    print(str(exc))
    assert "ID column 'NotACol' not found." in str(exc)


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


def test_AnnotateTaxResult_get_ident_split_but_keep_version():
    gA = {"name": "GCF_001881345.1 secondname"}
    taxres = AnnotateTaxResult(gA, keep_identifier_versions=True)
    print("raw ident: ", taxres.raw['name'])
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


def test_AnnotateTaxResult_get_ident_keep_full():
    gA = {"name": "GCF_001881345.1 secondname"}
    taxres = AnnotateTaxResult(gA, keep_full_identifiers=True)
    print("raw ident: ", taxres.raw['name'])
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
        check_and_load_gather_csvs(csvs, tax_assign)
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
    gather_results = check_and_load_gather_csvs(csvs, tax_assign, force=True)
    assert len(gather_results) == 1
    q_res = gather_results[0]
    assert len(q_res.raw_taxresults) == 4
    assert q_res.n_missed == 1
    assert 'gA' in q_res.missed_idents
    assert q_res.n_skipped == 0


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
        check_and_load_gather_csvs(csvs, tax_assign, fail_on_missing_taxonomy=True, force=True)
    assert "Failing, as requested via --fail-on-missing-taxonomy" in str(exc)


def test_load_gather_results():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    gather_csv = utils.get_test_data('tax/test1.gather.csv')
    gather_results, header = load_gather_results(gather_csv, tax_assignments=tax_assign)
    assert len(gather_results) == 1
    for query_name, res in gather_results.items():
        assert query_name == 'test1'
        assert len(res.raw_taxresults) == 4


def test_load_gather_results_gzipped(runtmp):
    gather_csv = utils.get_test_data('tax/test1.gather.csv')
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    gather_csv = utils.get_test_data('tax/test1.gather.csv')

    # rewrite gather_csv as gzipped csv
    gz_gather = runtmp.output('g.csv.gz')
    with open(gather_csv, 'rb') as f_in, gzip.open(gz_gather, 'wb') as f_out:
        f_out.writelines(f_in)
    #gather_results, header, seen_queries = load_gather_results(gz_gather)
    gather_results, header = load_gather_results(gz_gather, tax_assignments=tax_assign)
    assert len(gather_results) == 1
    for query_name, res in gather_results.items():
        assert query_name == 'test1'
        assert len(res.raw_taxresults) == 4


def test_load_gather_results_bad_header(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    g_csv = utils.get_test_data('tax/test1.gather.csv')

    bad_g_csv = runtmp.output('g.csv')

    #creates bad gather result
    bad_g = [x.replace("f_unique_to_query", "nope") for x in open(g_csv, 'r')]
    with open(bad_g_csv, 'w') as fp:
        for line in bad_g:
            fp.write(line)
    print("bad_gather_results: \n", bad_g)

    with pytest.raises(ValueError) as exc:
        gather_results, header = load_gather_results(bad_g_csv, tax_assignments=tax_assign)
    assert f"'{bad_g_csv}' is missing columns needed for taxonomic summarization" in str(exc.value)


def test_load_gather_results_empty(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv],
                                     keep_full_identifiers=False,
                                     keep_identifier_versions=False)
    empty_csv = runtmp.output('g.csv')

    #creates empty gather result
    with open(empty_csv, 'w') as fp:
        fp.write('')

    with pytest.raises(ValueError) as exc:
        gather_results, header = load_gather_results(empty_csv, tax_assignments=tax_assign)
    assert f"Cannot read gather results from '{empty_csv}'. Is file empty?" in str(exc.value)


def test_load_taxonomy_csv():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv])
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1']
    assert len(tax_assign) == 6 # should have read 6 rows


def test_load_taxonomy_csv_LIN():
    taxonomy_csv = utils.get_test_data('tax/test.LIN-taxonomy.csv')
    tax_assign = MultiLineageDB.load([taxonomy_csv], lins=True)
    print("taxonomy assignments: \n", tax_assign)
    assert list(tax_assign.keys()) == ['GCF_001881345.1', 'GCF_009494285.1', 'GCF_013368705.1', 'GCF_003471795.1', 'GCF_000017325.1', 'GCF_000021665.1']
    #assert list(tax_assign.keys()) == ["GCF_000010525.1", "GCF_000007365.1", "GCF_000007725.1", "GCF_000009605.1", "GCF_000021065.1", "GCF_000021085.1"]
    assert len(tax_assign) == 6 # should have read 6 rows
    print(tax_assign.available_ranks)
    assert tax_assign.available_ranks == {str(x) for x in range(0,20)}


def test_load_taxonomy_csv_LIN_fail():
    taxonomy_csv = utils.get_test_data('tax/test.taxonomy.csv')
    with pytest.raises(ValueError) as exc:
        MultiLineageDB.load([taxonomy_csv], lins=True)
    assert f"'lin' column not found: cannot read LIN taxonomy assignments from {taxonomy_csv}." in str(exc.value)


def test_load_taxonomy_csv_LIN_mismatch_in_taxfile(runtmp):
    taxonomy_csv = utils.get_test_data('tax/test.LIN-taxonomy.csv')
    mimatchLIN_csv = runtmp.output('mmLIN-taxonomy.csv')
    with open(mimatchLIN_csv, 'w') as mm:
        tax21=[]
        tax = [x.rstrip() for x in open(taxonomy_csv, 'r')]
        for n, taxline in enumerate(tax):
            if n == 2: # add ;0 to a LIN
                taxlist = taxline.split(',')
                taxlist[1] += ';0' # add 21st position to LIN
                tax21.append(",".join(taxlist))
            else:
                tax21.append(taxline)
        mm.write("\n".join(tax21))
    with pytest.raises(ValueError) as exc:
        MultiLineageDB.load([mimatchLIN_csv], lins=True)
    assert "For taxonomic summarization, all LIN assignments must use the same number of LIN positions." in str(exc.value)


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


def test_format_for_krona_summarization():
    """test format for krona"""
    # make gather results
     # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.2,'f_unique_to_query': 0.2,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, summarize=True, single_query=True)
    kres, header = format_for_krona([q_res], 'superkingdom')
    assert header == ['fraction', 'superkingdom']
    print("krona_res: ", kres)
    assert kres == [(0.5, 'a'), (0.5, 'unclassified')]
    kres, header = format_for_krona([q_res], 'phylum')
    assert header == ['fraction', 'superkingdom', 'phylum']
    assert kres == [(0.3, 'a', 'c'), (0.2, 'a', 'b'), (0.5, 'unclassified', 'unclassified')]


def test_format_for_krona_classification():
    """test format for krona"""
    # make gather results
     # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.2,'f_unique_to_query': 0.2,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, classify=True, single_query=True)
    kres, header = format_for_krona([q_res], 'superkingdom', classification=True)
    assert header == ['fraction', 'superkingdom']
    print("krona_res: ", kres)
    assert kres == [(0.5, 'a')]#, (0.5, 'unclassified')]
    kres, header = format_for_krona([q_res], 'phylum', classification=True)
    assert header == ['fraction', 'superkingdom', 'phylum']
    assert kres == [(0.3, 'a', 'c')]#, (0.7, 'unclassified', 'unclassified')]


def test_format_for_krona_improper_rank():
    """test format for krona"""
    # make gather results
     # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.2,'f_unique_to_query': 0.2,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, summarize=True, single_query=True)
    with pytest.raises(ValueError) as exc:
        format_for_krona([q_res], 'NotARank')
    print(str(exc))
    assert "Rank 'NotARank' not present in summarized ranks." in str(exc)


def test_format_for_krona_summarization_two_queries():
    """test format for krona with multiple queries (normalize by n_queries)"""
    # make gather results
     # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax,gB_tax])

    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.2,'f_unique_to_query': 0.2,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30},
                      {'query_name': 'queryB', "name": 'gB', 'f_unique_weighted': 0.5,'f_unique_to_query': 0.5,'unique_intersect_bp': 50}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, summarize=True)
    kres, header = format_for_krona(list(gres.values()), 'superkingdom')
    assert header == ['fraction', 'superkingdom']
    print("krona_res: ", kres)
    assert kres == [(0.5, 'a'), (0.5, 'unclassified')]
    kres, header = format_for_krona(list(gres.values()), 'phylum')
    assert header == ['fraction', 'superkingdom', 'phylum']
    assert kres == [(0.4, 'a', 'c'), (0.1, 'a', 'b'), (0.5, 'unclassified', 'unclassified')]


def test_write_krona(runtmp):
    """test two matches, equal f_unique_to_query"""
    krona_results =  [(0.5, 'a', 'b', 'c'), (0.5, 'a', 'b', 'd')]
    header = ['fraction', 'superkingdom', 'phylum', 'class']
    outk= runtmp.output("outkrona.tsv")
    with open(outk, 'w') as out_fp:
        write_krona(header, krona_results, out_fp)

    kr = [x.strip().split('\t') for x in open(outk, 'r')]
    print("krona_results_from_file: \n", kr)
    assert kr[0] == ["fraction", "superkingdom", "phylum", "class"]
    assert kr[1] == ["0.5", "a", "b", "c"]
    assert kr[2] == ["0.5", "a", "b", "d"]


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
    sk_lineage='a'
    print(sk_lineage)
    sk_linD = {sk_lineage: {'sample1': '0.500' ,'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, sk_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a', '0.500', '0.700']]

    phy_lineage='a;b'
    print(phy_lineage)
    phy2_lineage = 'a;c'
    print(phy2_lineage)
    phy_linD = {phy_lineage: {'sample1': '0.500'}, phy2_lineage: {'sample2': '0.700'}}
    with open(outfrac, 'w') as out_fp:
        write_lineage_sample_frac(sample_names, phy_linD, out_fp)

    frac_lines = [x.strip().split('\t') for x in open(outfrac, 'r')]
    print("csv_lines: ", frac_lines)
    assert frac_lines == [['lineage', 'sample1', 'sample2'], ['a;b', '0.500', '0'],  ['a;c', '0', '0.700']]


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
    assert set(RankLineageInfo().taxlist) == set(db.available_ranks)

    db = MultiLineageDB.load([taxonomy_csv, taxonomy_db],
                             keep_full_identifiers=False,
                             keep_identifier_versions=False)
    assert len(db.shadowed_identifiers()) == 6
    assert set(RankLineageInfo().taxlist[:-1]) == set(db.available_ranks)


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
    assert lin.rank=="rank1"
    assert lin.name =="name1"
    assert lin.taxid==None


def test_LineagePair_1():
    lin = LineagePair(rank="rank1", name='name1', taxid=1)
    assert lin.rank=="rank1"
    assert lin.name =="name1"
    assert lin.taxid==1
    print(lin)


def test_BaseLineageInfo_init_empty():
    ranks=["A", "B", "C"]
    taxinf = BaseLineageInfo(ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['', '', ''] # this is a bit odd, but it's what preserves empty ranks...
    print(taxinf.filled_lineage)
    assert taxinf.filled_lineage == ()
    assert taxinf.lowest_lineage_name == None
    assert taxinf.lowest_lineage_taxid == None
    assert taxinf.filled_ranks == ()
    assert taxinf.name_at_rank("A") == None
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
    assert taxinf.name_at_rank("A") == "a"


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
    lin_tups = (LineagePair(rank="A", name='a'), LineagePair(rank="C", name='b'))
    taxinf = BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', '', 'b']


def test_BaseLineageInfo_init_no_ranks():
    x = "a;b;c"
    rankD = {"superkingdom": "a", "phylum": "b", "class": "c"}
    lin_tups = (LineagePair(rank="rank2", name='name1'), LineagePair(rank="rank1", name='name1'))
    with pytest.raises(TypeError) as exc:
        BaseLineageInfo(lineage_str=x)
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


def test_BaseLineageInfo_init_not_lineagepair():
    ranks=["A", "B", "C"]
    lin_tups = (("rank1", "name1"),)
    with pytest.raises(ValueError) as exc:
        BaseLineageInfo(lineage=lin_tups, ranks=ranks)
    print(str(exc))
    assert "is not tax_utils LineagePair" in str(exc)


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


def test_LINLineageInfo_init_empty():
    taxinf = LINLineageInfo()
    assert taxinf.n_lin_positions == 0
    assert taxinf.zip_lineage()== []
    assert taxinf.display_lineage()== ""
    assert taxinf.filled_ranks == ()
    assert taxinf.n_filled_pos == 0


def test_LINLineageInfo_init_n_pos():
    n_pos = 5
    taxinf = LINLineageInfo(n_lin_positions=n_pos)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.n_lin_positions == 5
    assert taxinf.zip_lineage()== ['', '', '', '', '']
    assert taxinf.filled_ranks == ()
    assert taxinf.n_filled_pos == 0


def test_LINLineageInfo_init_n_pos_and_lineage_str():
    x = "0;0;1"
    n_pos = 5
    taxinf = LINLineageInfo(lineage_str=x, n_lin_positions=n_pos)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.n_lin_positions == 5
    assert taxinf.zip_lineage()== ['0', '0', '1', '', '']
    assert taxinf.filled_ranks == ("0","1","2")
    assert taxinf.n_filled_pos == 3


def test_LINLineageInfo_init_n_pos_and_lineage_str_fail():
    x = "0;0;1"
    n_pos = 2
    with pytest.raises(ValueError) as exc:
        LINLineageInfo(lineage_str=x, n_lin_positions=n_pos)
    print(str(exc))
    assert "Provided 'n_lin_positions' has fewer positions than provided 'lineage_str'." in str(exc)


def test_LINLineageInfo_init_lineage_str_only():
    x = "0,0,1"
    taxinf = LINLineageInfo(lineage_str=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.n_lin_positions == 3
    assert taxinf.zip_lineage()== ['0', '0', '1']
    assert taxinf.filled_ranks == ("0","1","2")
    assert taxinf.n_filled_pos == 3


def test_LINLineageInfo_init_not_lineagepair():
    lin_tups = (("rank1", "name1"),)
    with pytest.raises(ValueError) as exc:
        LINLineageInfo(lineage=lin_tups)
    print(str(exc))
    assert "is not tax_utils LineagePair" in str(exc)


def test_LINLineageInfo_init_lineagepair():
    lin_tups = (LineagePair("rank1", "name1"), LineagePair("rank2", None),)
    taxinf = LINLineageInfo(lineage=lin_tups)
    print(taxinf.lineage)
    assert taxinf.n_lin_positions == 2
    assert taxinf.zip_lineage()== ["name1", ""]
    assert taxinf.zip_lineage(truncate_empty=True)== ["name1"]
    assert taxinf.filled_ranks == ("rank1",)
    assert taxinf.ranks == ("rank1", "rank2")
    assert taxinf.n_filled_pos == 1


def test_lca_LINLineageInfo_diff_n_pos():
    x = "0;0;1"
    y = '0'
    lin1 = LINLineageInfo(lineage_str=x)
    lin2 = LINLineageInfo(lineage_str=y)
    assert lin1.is_compatible(lin2)
    assert lin2.is_compatible(lin1)
    lca_from_lin1 = lin1.find_lca(lin2)
    lca_from_lin2 = lin2.find_lca(lin1)
    assert lca_from_lin1 == lca_from_lin2
    assert lca_from_lin1.display_lineage(truncate_empty=True) == "0"


def test_lca_LINLineageInfo_no_lca():
    x = "0;0;1"
    y = '12;0;1'
    lin1 = LINLineageInfo(lineage_str=x)
    lin2 = LINLineageInfo(lineage_str=y)
    assert lin1.is_compatible(lin2)
    assert lin2.is_compatible(lin1)
    lca_from_lin1 = lin1.find_lca(lin2)
    lca_from_lin2 = lin2.find_lca(lin1)
    assert lca_from_lin1 == lca_from_lin2 == None


def test_lca_RankLineageInfo_no_lca():
    x = "a;b;c"
    y = 'd;e;f;g'
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    assert lin1.is_compatible(lin2)
    assert lin2.is_compatible(lin1)
    lca_from_lin1 = lin1.find_lca(lin2)
    lca_from_lin2 = lin2.find_lca(lin1)
    assert lca_from_lin1 == lca_from_lin2 == None


def test_incompatibility_LINLineageInfo_RankLineageInfo():
    x="a;b;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = LINLineageInfo(lineage_str=x)
    assert not lin1.is_compatible(lin2)
    assert not lin2.is_compatible(lin1)


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


def test_RankLineageInfo_init_lineage_dict_fail():
    ranks=["A", "B", "C"]
    lin_tups = (LineagePair(rank="A", name='a'), LineagePair(rank="C", name='b'))
    with pytest.raises(ValueError) as exc:
        taxinf = RankLineageInfo(ranks=ranks, lineage_dict=lin_tups)
    print(str(exc))

    assert "is not dictionary" in str(exc)


def test_RankLineageInfo_init_lineage_dict():
    x = {'rank1': 'name1', 'rank2': 'name2'}
    taxinf = RankLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', 'name2']


def test_RankLineageInfo_init_lineage_dict_default_ranks():
    x = {"superkingdom":'a',"phylum":'b'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print(taxinf.lineage)
    print(taxinf.lineage_str)
    assert taxinf.zip_lineage()== ['a', 'b', '', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_withtaxpath():
    x = {'rank1': 'name1', 'rank2': 'name2', 'taxpath': "1|2"}
    taxinf = RankLineageInfo(lineage_dict=x, ranks=["rank1", "rank2"])
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    print("zipped taxids: ", taxinf.zip_taxid())
    assert taxinf.zip_lineage()== ['name1', 'name2']
    assert taxinf.zip_taxid()== ['1', '2']
    assert taxinf.lowest_lineage_taxid == "2"
    assert taxinf.lowest_lineage_name == "name2"


def test_RankLineageInfo_init_lineage_str_lineage_dict_test_eq():
    x = "a;b;c"
    ranks=["A", "B", "C"]
    rankD = {"A": "a", "B": "b", "C": "c"}
    lin1 = RankLineageInfo(lineage_str=x, ranks=ranks)
    lin2 = RankLineageInfo(lineage_dict=rankD, ranks=ranks)
    assert lin1 == lin2


def test_RankLineageInfo_init_lineage_dict_missing_rank():
    x = {'superkingdom': 'name1', 'class': 'name2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']
    assert taxinf.zip_lineage(truncate_empty=True)== ['name1', '', 'name2']


def test_RankLineageInfo_init_lineage_dict_missing_rank_with_taxpath():
    x = {'superkingdom': 'name1', 'class': 'name2', 'taxpath': '1||2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', 'name2', '', '', '', '', '']
    assert taxinf.zip_taxid()== ['1', '', '2', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_name_taxpath_mismatch():
    # If there's no name, we don't report the taxpath, because lineage is not "filled".
    # Is this desired behavior?
    x = {'superkingdom': 'name1', 'taxpath': '1||2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    assert taxinf.zip_lineage()== ['name1', '', '', '', '', '', '', '']
    assert taxinf.zip_taxid()== ['1', '', '', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_name_taxpath_missing_taxids():
    # If there's no name, we don't report the taxpath, because lineage is not "filled".
    # Is this desired behavior?
    x = {'superkingdom': 'name1', 'phylum': "name2", "class": "name3", 'taxpath': '|2'}
    taxinf = RankLineageInfo(lineage_dict=x)
    print("ranks: ", taxinf.ranks)
    print("lineage: ", taxinf.lineage)
    print("zipped lineage: ", taxinf.zip_lineage())
    print("zipped taxids: ", taxinf.zip_taxid())
    assert taxinf.zip_lineage()== ['name1', 'name2', 'name3', '', '', '', '', '']
    assert taxinf.zip_taxid()== ['', '2', '', '', '', '', '', '']


def test_RankLineageInfo_init_lineage_dict_taxpath_too_long():
    x = {'superkingdom': 'name1', 'class': 'name2', 'taxpath': '1||2||||||||||'}
    with pytest.raises(ValueError) as exc:
        RankLineageInfo(lineage_dict=x)
    print(str(exc))
    assert f"Number of NCBI taxids (13) exceeds number of ranks (8)" in str(exc)


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
    assert lin1.is_compatible(lin2)
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

    lca_from_lin1 = lin1.find_lca(lin2)
    print(lca_from_lin1.display_lineage())
    lca_from_lin2 = lin2.find_lca(lin1)
    assert lca_from_lin1 == lca_from_lin2
    assert lca_from_lin1.display_lineage() == "d__a;p__b;c__c;o__d"
    


def test_is_lineage_match_2():
    # match at family, and above, levels; no genus or species to match
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    assert lin1.is_compatible(lin2)
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

    lca_from_lin1 = lin1.find_lca(lin2)
    print(lca_from_lin1.display_lineage())
    lca_from_lin2 = lin2.find_lca(lin1)
    assert lca_from_lin1 == lca_from_lin2
    assert lca_from_lin1.display_lineage() == "d__a;p__b;c__c;o__d;f__f"


def test_is_lineage_match_3():
    # one lineage is empty
    lin1 = RankLineageInfo()
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')

    assert lin1.is_compatible(lin2)
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
    assert not lin1.is_compatible(lin2)
    with pytest.raises(ValueError) as exc:
        lin1.is_lineage_match(lin2, 'superkingdom')
    print(str(exc))
    assert 'Cannot compare lineages from taxonomies with different ranks.' in str(exc)


def test_is_lineage_match_improper_rank():
    #test comparison with incompatible ranks
    lin1 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__e')
    lin2 = RankLineageInfo(lineage_str = 'd__a;p__b;c__c;o__d;f__f')
    print(lin1.lineage)
    assert lin1.is_compatible(lin2)
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


def test_TaxResult_get_match_lineage_1():
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD)
    assert taxres.lineageInfo.display_lineage() == "a;b;c"


def test_AnnotateTaxResult_get_match_lineage_1():
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = AnnotateTaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD)
    assert taxres.lineageInfo.display_lineage() == "a;b;c"
    assert taxres.row_with_lineages() == {"name": "gA.1 name", "lineage": "a;b;c"}


def test_TaxResult_get_match_lineage_skip_ident():
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gA'])
    print("skipped_ident?: ", taxres.skipped_ident)
    print("missed_ident?: ", taxres.missed_ident)
    assert taxres.skipped_ident == True
    assert taxres.lineageInfo == RankLineageInfo()
    assert taxres.lineageInfo.display_lineage() == ""
    assert taxres.lineageInfo.display_lineage(null_as_unclassified=True) == "unclassified"


def test_TaxResult_get_match_lineage_missed_ident_fail_on_missing():
    gA_tax = ("gA.1", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gB'])
    print("skipped_ident?: ", taxres.skipped_ident)
    print("missed_ident?: ", taxres.missed_ident)
    assert taxres.skipped_ident == False
    assert taxres.missed_ident == True
    assert taxres.lineageInfo == RankLineageInfo()
    assert taxres.lineageInfo.display_lineage() == ""
    assert taxres.lineageInfo.display_lineage(null_as_unclassified=True) == "unclassified"


def test_TaxResult_get_match_lineage_missed_ident_fail_on_missing():
    gA_tax = ("gA.1", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])

    gA = {"name": "gA.1 name"}
    taxres = make_TaxResult(gA)
    with pytest.raises(ValueError) as exc:
        taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gB'], fail_on_missing_taxonomy=True)
    print(str(exc))
    assert "Error: ident 'gA' is not in the taxonomy database." in str(exc)


def test_QueryTaxResult():
    "basic functionality: initialize and add a taxresult"
    tax_info = [("gA", "a;b;c")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    taxres = make_TaxResult(taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    assert q_res.ranks == []
    assert q_res.ascending_ranks == []
    q_res.add_taxresult(taxres)
    # check that new querytaxres is compatible with taxres
    assert q_res.is_compatible(taxres)
    # check that a few thngs were set properly and/or are not yet set.
    assert q_res.query_name == "q1"
    assert q_res.query_info.query_bp == 100
    assert len(q_res.raw_taxresults) == 1
    assert q_res.skipped_idents == set()
    assert q_res.missed_idents == set()
    assert q_res.summarized_lineage_results == {}
    taxranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
    assert q_res.ranks == taxranks
    assert q_res.ascending_ranks == taxranks[::-1]


def test_QueryTaxResult_add_incompatible():
    "initialize and try to add incompatible taxresult"
    tax_info = [("gA", "a;b;c")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    taxres = make_TaxResult(taxD=taxD)
    taxres2 = make_TaxResult({'query_name': 'q2'}, taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # check that new querytaxres is compatible with taxres and not taxres2
    assert q_res.is_compatible(taxres)
    assert not q_res.is_compatible(taxres2)
    q_res.add_taxresult(taxres)
    with pytest.raises(ValueError) as exc:
        q_res.add_taxresult(taxres2)
    print(str(exc))
    assert "Error: Cannot add TaxResult: query information does not match." in str(exc)


def test_QueryTaxResult_add_without_tax_info():
    "initialize and add a taxresult with missed ident"
    taxres = make_TaxResult() # do not add taxonomic info
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    print("attempted to add lineage info?: ", taxres.match_lineage_attempted)
    with pytest.raises(ValueError) as exc:
        q_res.add_taxresult(taxres)
    print(str(exc))
    assert "Error: Cannot add TaxResult. Please use get_match_lineage() to add taxonomic lineage information first." in str(exc)
    
    
def test_QueryTaxResult_add_skipped_ident():
    "initialize and add a taxresult with skipped ident"
    gA_tax = ("gA", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])
    taxres = make_TaxResult(taxD=taxD, skip_idents = ['gA'])
#    taxres.get_match_lineage(tax_assignments=taxD, skip_idents=['gA'])
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    q_res.add_taxresult(taxres)
    assert len(q_res.skipped_idents) == 1
    assert len(q_res.raw_taxresults) == 1
    assert q_res.missed_idents == set()
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_add_missed_ident():
    "initialize and add a taxresult with missed ident"
    gA_tax = ("gB", "a;b;c")
    taxD = make_mini_taxonomy([gA_tax])
    taxres = make_TaxResult(taxD=taxD)
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # add taxonomic info to taxres
    q_res.add_taxresult(taxres)
    assert len(q_res.missed_idents) == 1
    assert len(q_res.raw_taxresults) == 1
    assert q_res.skipped_idents == set()
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_track_missed_and_skipped():
    "make sure missed and skipped idents are being tracked"
    # make taxonomy
    tax_info = [("gA", "a;b;c"), ("gB", "a;b;d")]
    taxD = make_mini_taxonomy(tax_info=tax_info)
    # make results
    taxres = make_TaxResult()
    taxres2 = make_TaxResult({"name": 'gB'}) # skipped
    taxres3 = make_TaxResult({"name": 'gB'}) # skipped
    taxres4 = make_TaxResult({"name": 'gC'}) # skipped
    taxres5 = make_TaxResult({"name": 'gD'}) # missed
    taxres6 = make_TaxResult({"name": 'gE'}) # missed
    # initialize
    q_res = QueryTaxResult(taxres.query_info)
    # add taxonomic info to taxres, add to q_res
    for n, tr in enumerate([taxres, taxres2, taxres3, taxres4, taxres5, taxres6]):
        tr.get_match_lineage(tax_assignments=taxD, skip_idents=['gB', 'gC'])
        print("num: ", n)
        print("skipped?: ", tr.skipped_ident)
        print("missed?: ", tr.missed_ident)
        q_res.add_taxresult(tr)
    assert len(q_res.raw_taxresults) == 6
    print(q_res.n_skipped)
    print(q_res.n_missed)
    assert q_res.n_missed == 2
    assert q_res.n_skipped == 3
    assert 'gB' in q_res.skipped_idents
    assert len(q_res.skipped_idents) == 2
    assert 'gD' in q_res.missed_idents
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_track_missed_and_skipped_using_fn():
    "make sure missed and skipped idents are being tracked. Same as above but use helper fn."
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}, {"name": 'gB'}, {"name": 'gC'}, {"name": 'gD'}, {"name": 'gE'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, skip_idents=['gB', 'gC'])
    # should have 6 results for default query 'q1' 
    print(gres.keys())
    q_res = next(iter(gres.values()))
    assert len(q_res.raw_taxresults) == 6
    print(q_res.n_skipped)
    print(q_res.n_missed)
    assert q_res.n_missed == 2
    assert q_res.n_skipped == 3
    assert 'gB' in q_res.skipped_idents
    assert len(q_res.skipped_idents) == 2
    assert 'gD' in q_res.missed_idents
    assert q_res.summarized_lineage_results == {}


def test_QueryTaxResult_summarize_up_ranks_1():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    #print(q_res.sum_uniq_weighted.values())
    #print(q_res.sum_uniq_weighted['superkingdom'])
    assert list(q_res.sum_uniq_weighted.keys()) == ['class', 'phylum', 'superkingdom']
    assert q_res.sum_uniq_weighted['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.4)}
    assert q_res.sum_uniq_to_query['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.2)}
    assert q_res.sum_uniq_bp['superkingdom'] == {RankLineageInfo(lineage_str="a"): 40}
    assert q_res.sum_uniq_weighted['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.4)}
    assert q_res.sum_uniq_to_query['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.2)}
    assert q_res.sum_uniq_bp['phylum'] == {RankLineageInfo(lineage_str="a;b"): 40}
    assert q_res.sum_uniq_weighted['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.2),
                                                RankLineageInfo(lineage_str="a;b;d"): approx(0.2)}
    assert q_res.sum_uniq_to_query['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.1),
                                                RankLineageInfo(lineage_str="a;b;d"): approx(0.1)}
    assert q_res.sum_uniq_bp['class'] == {RankLineageInfo(lineage_str="a;b;c"): 20,
                                          RankLineageInfo(lineage_str="a;b;d"): 20}


def test_QueryTaxResult_summarize_up_ranks_2():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.3)}
    assert q_res.sum_uniq_to_query['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.15)}
    assert q_res.sum_uniq_bp['superkingdom'] == {RankLineageInfo(lineage_str="a"): 30}
    assert q_res.sum_uniq_weighted['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.3)}
    assert q_res.sum_uniq_to_query['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.15)}
    assert q_res.sum_uniq_bp['phylum'] == {RankLineageInfo(lineage_str="a;b"): 30}
    assert q_res.sum_uniq_weighted['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.2),
                                                RankLineageInfo(lineage_str="a;b;d"): approx(0.1)}
    assert q_res.sum_uniq_to_query['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.1),
                                                RankLineageInfo(lineage_str="a;b;d"): approx(0.05)}
    assert q_res.sum_uniq_bp['class'] == {RankLineageInfo(lineage_str="a;b;c"): 20,
                                          RankLineageInfo(lineage_str="a;b;d"): 10}


def test_QueryTaxResult_summarize_up_ranks_missing_lineage():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    #print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.2)}
    assert q_res.sum_uniq_to_query['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.1)}
    assert q_res.sum_uniq_bp['superkingdom'] == {RankLineageInfo(lineage_str="a"): 20}
    assert q_res.sum_uniq_weighted['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.2)}
    assert q_res.sum_uniq_to_query['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.1)}
    assert q_res.sum_uniq_bp['phylum'] == {RankLineageInfo(lineage_str="a;b"): 20}
    assert q_res.sum_uniq_weighted['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.2)}
    assert q_res.sum_uniq_to_query['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.1)}
    assert q_res.sum_uniq_bp['class'] == {RankLineageInfo(lineage_str="a;b;c"): 20}


def test_QueryTaxResult_summarize_up_ranks_skipped_lineage():
    "basic functionality: summarize up ranks"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, skip_idents=['gB'])
    assert len(gres.keys()) == 1
    q_res = next(iter(gres.values()))
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == ['class', 'phylum', 'superkingdom']
    #print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['superkingdom'])
    assert q_res.sum_uniq_weighted['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.2)}
    assert q_res.sum_uniq_to_query['superkingdom'] == {RankLineageInfo(lineage_str="a"): approx(0.1)}
    assert q_res.sum_uniq_bp['superkingdom'] == {RankLineageInfo(lineage_str="a"): 20}
    assert q_res.sum_uniq_weighted['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.2)}
    assert q_res.sum_uniq_to_query['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.1)}
    assert q_res.sum_uniq_bp['phylum'] == {RankLineageInfo(lineage_str="a;b"): 20}
    assert q_res.sum_uniq_weighted['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.2)}
    assert q_res.sum_uniq_to_query['class'] == {RankLineageInfo(lineage_str="a;b;c"): approx(0.1)}
    assert q_res.sum_uniq_bp['class'] == {RankLineageInfo(lineage_str="a;b;c"): 20}


def test_QueryTaxResult_summarize_up_ranks_perfect_match():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{'f_unique_to_query': 1.0}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    assert len(q_res.raw_taxresults) == 1
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_to_query['superkingdom'])
    assert list(q_res.sum_uniq_to_query['superkingdom'].values()) == [1.0]
    assert 'gA' in q_res.perfect_match


def test_QueryTaxResult_summarize_up_ranks_already_summarized():
    "summarize up ranks: error, already summarized"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{'f_unique_to_query': 1.0}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    with pytest.raises(ValueError) as exc:
        q_res.summarize_up_ranks()
    print(str(exc))
    assert "Error: already summarized" in str(exc)
    

def test_QueryTaxResult_summarize_up_ranks_already_summarized_force():
    "summarize up ranks: already summarized but force"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks()
    q_res.summarize_up_ranks(force_resummarize=True)
    assert list(q_res.sum_uniq_weighted.keys()) == ['class', 'phylum', 'superkingdom']

    #check that all results are still good 
    assert len(q_res.raw_taxresults) == 2
    assert q_res.sum_uniq_weighted['superkingdom'] ==  {RankLineageInfo(lineage_str="a"): approx(0.3)}
    assert q_res.sum_uniq_weighted['phylum'] ==  {RankLineageInfo(lineage_str="a;b"): approx(0.3)}
    assert q_res.sum_uniq_to_query['phylum'] ==  {RankLineageInfo(lineage_str="a;b"): approx(0.15)}
    assert q_res.sum_uniq_bp['phylum'] ==  {RankLineageInfo(lineage_str="a;b"): 30}
    assert q_res.sum_uniq_to_query['class'] ==  {RankLineageInfo(lineage_str="a;b;c"): approx(0.1),
                                                 RankLineageInfo(lineage_str="a;b;d"): approx(0.05)}
    assert q_res.sum_uniq_weighted['class'] ==  {RankLineageInfo(lineage_str="a;b;c"): approx(0.2),
                                                 RankLineageInfo(lineage_str="a;b;d"): approx(0.1)}
    assert q_res.sum_uniq_bp['class'] ==  {RankLineageInfo(lineage_str="a;b;c"): 20,
                                           RankLineageInfo(lineage_str="a;b;d"): 10}


def test_QueryTaxResult_summarize_up_ranks_single_rank():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    q_res.summarize_up_ranks(single_rank='phylum')
    assert len(q_res.raw_taxresults) == 2
    assert list(q_res.sum_uniq_weighted.keys()) == ['phylum']
    print(q_res.sum_uniq_weighted.keys())
    print(q_res.sum_uniq_weighted.values())
    print(q_res.sum_uniq_weighted['phylum'])
    assert q_res.sum_uniq_weighted['phylum'] == {RankLineageInfo(lineage_str="a;b"): approx(0.3)}
    assert list(q_res.sum_uniq_to_query['phylum'].values()) == [approx(0.15)]                                                    
    assert list(q_res.sum_uniq_bp['phylum'].values()) == [30]                                                    
    assert q_res.summarized_ranks == ['phylum']

def test_QueryTaxResult_summarize_up_ranks_single_rank_not_available():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    with pytest.raises(ValueError) as exc:
        q_res.summarize_up_ranks(single_rank='NotARank')
    print(str(exc))
    assert "Error: rank 'NotARank' not in available ranks (strain, species, genus, family, order, class, phylum, superkingdom)" in str(exc)


def test_QueryTaxResult_summarize_up_ranks_single_rank_not_filled():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_weighted': 0.1,'f_unique_to_query': 0.05,'unique_intersect_bp': 10,}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    with pytest.raises(ValueError) as exc:
        q_res.summarize_up_ranks(single_rank='species')
    print(str(exc))
    assert "Error: rank 'species' was not available for any matching lineages." in str(exc)


def test_QueryTaxResult_build_summarized_result_1():
    "basic functionality: build summarized_result"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    sk = [SummarizedGatherResult(rank='superkingdom', fraction=0.2, f_weighted_at_rank=0.4, 
                             lineage=RankLineageInfo(lineage_str='a'),
                             bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2)),
          SummarizedGatherResult(rank='superkingdom', fraction=0.8, f_weighted_at_rank=0.6,
                             lineage=RankLineageInfo(), bp_match_at_rank=60, query_ani_at_rank=None)]
    print(q_res.summarized_lineage_results['superkingdom'])
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(rank='phylum', fraction=0.2, f_weighted_at_rank=0.4, 
                                   lineage=RankLineageInfo(lineage_str='a;b'),
                                   bp_match_at_rank=40, query_ani_at_rank=approx(0.95, rel=1e-2)),
            SummarizedGatherResult(rank='phylum', fraction=0.8, f_weighted_at_rank=0.6,
                                   lineage=RankLineageInfo(), bp_match_at_rank=60, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(rank='class', fraction=0.1, f_weighted_at_rank=0.2, 
                                 lineage=RankLineageInfo(lineage_str='a;b;c'), 
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.93, rel=1e-2)),
          SummarizedGatherResult(rank='class', fraction=0.1, f_weighted_at_rank=0.2,
                                 lineage=RankLineageInfo(lineage_str='a;b;d'),
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.93, rel=1e-2)),
          SummarizedGatherResult(rank='class', fraction=0.8, f_weighted_at_rank=0.6,
                                 lineage=RankLineageInfo(), bp_match_at_rank=60, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.4)
    assert q_res.total_f_classified['class'] == approx(0.2)
    assert q_res.total_bp_classified['superkingdom'] == 40


def test_QueryTaxResult_build_summarized_result_2():
    """test two queries, build summarized result for each"""
    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax, gB_tax])
    # make gather results
    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.5,'f_unique_to_query': 0.5,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.4,'f_unique_to_query': 0.3,'unique_intersect_bp': 30},
                      {'query_name': 'queryB', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD)
    
    for query_name, q_res in gres.items():
        q_res.build_summarized_result() # summarize and build result
        sk = q_res.summarized_lineage_results['superkingdom']
        phy = q_res.summarized_lineage_results['phylum']
        assert len(sk) == 2
        assert sk[0].lineage == RankLineageInfo(lineage_str="a")
        print(phy)
        if query_name == 'queryA':
            # check superkingdom results
            assert sk[0].fraction == approx(0.8)
            assert sk[0].f_weighted_at_rank == approx(0.9)
            assert sk[0].bp_match_at_rank == 80
            assert sk[1].fraction == approx(0.2)
            assert sk[1].f_weighted_at_rank == approx(0.1)
            assert sk[1].bp_match_at_rank == 20
            assert sk[1].lineage == RankLineageInfo()
            # check phylum results
            assert len(phy) == 3
            assert phy[0].fraction == approx(0.5)
            assert phy[0].f_weighted_at_rank == approx(0.5)
            assert phy[0].bp_match_at_rank == 50
            assert phy[0].lineage ==  RankLineageInfo(lineage_str="a;b")
            assert phy[1].fraction == approx(0.3)
            assert phy[1].f_weighted_at_rank == approx(0.4)
            assert phy[1].bp_match_at_rank == 30
            assert phy[1].lineage ==  RankLineageInfo(lineage_str="a;c")
            assert phy[2].fraction == approx(0.2)
            assert phy[2].f_weighted_at_rank == approx(0.1)
            assert phy[2].bp_match_at_rank == 20
            assert phy[2].lineage == RankLineageInfo()
        if query_name == 'queryB':
            # check superkingdom results
            assert sk[0].fraction == approx(0.3)
            assert sk[0].f_weighted_at_rank == approx(0.3)
            assert sk[0].bp_match_at_rank == 30
            assert sk[1].fraction == approx(0.7)
            assert sk[1].f_weighted_at_rank == approx(0.7)
            assert sk[1].bp_match_at_rank == 70
            assert sk[1].lineage == RankLineageInfo()
            # check phylum results
            assert len(phy) == 2
            assert phy[0].fraction == approx(0.3)
            assert phy[0].f_weighted_at_rank == approx(0.3)
            assert phy[0].bp_match_at_rank == 30
            assert phy[0].lineage ==  RankLineageInfo(lineage_str="a;c")
            assert phy[1].fraction == approx(0.7)
            assert phy[1].f_weighted_at_rank == approx(0.7)
            assert phy[1].bp_match_at_rank == 70
            assert phy[1].lineage == RankLineageInfo()


def test_QueryTaxResult_build_summarized_result_missing_lineage():
    "build summarized_result with missing lineage"
    taxD = make_mini_taxonomy([("gA", "a;b;c")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    print(q_res.summarized_lineage_results['superkingdom'])

    sk = [SummarizedGatherResult(rank='superkingdom', fraction=0.1, f_weighted_at_rank=0.2, 
                                 lineage=RankLineageInfo(lineage_str="a"), 
                                 bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)),
          SummarizedGatherResult(rank='superkingdom', fraction=0.9, lineage=RankLineageInfo(),f_weighted_at_rank=0.8,
                                 bp_match_at_rank=80, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(rank='phylum', fraction=0.1, f_weighted_at_rank=0.2,
                                  lineage=RankLineageInfo(lineage_str="a;b"), 
                                  bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)),
           SummarizedGatherResult(rank='phylum', fraction=0.9, lineage=RankLineageInfo(),f_weighted_at_rank=0.8,
                                  bp_match_at_rank=80, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(rank='class', fraction=0.1, lineage= RankLineageInfo(lineage_str="a;b;c"),
                                  f_weighted_at_rank=0.2, bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)),
          SummarizedGatherResult(rank='class', fraction=0.9, lineage=RankLineageInfo(), f_weighted_at_rank=0.8,
                                 bp_match_at_rank=80, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.2)
    assert q_res.total_f_classified['class'] == approx(0.1)
    assert q_res.total_bp_classified['superkingdom'] == 20


def test_QueryTaxResult_build_summarized_result_skipped_lineage():
    "build summarized_result with skipped lineage"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, skip_idents=['gB'])
    q_res.build_summarized_result()
    print(q_res.summarized_lineage_results.keys())
    print(q_res.summarized_lineage_results['superkingdom'])

    sk = [SummarizedGatherResult(rank='superkingdom', fraction=0.1, f_weighted_at_rank=0.2,  
                                 lineage=RankLineageInfo(lineage_str="a"), 
                                 bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)), 
          SummarizedGatherResult(rank='superkingdom', fraction=0.9, lineage=RankLineageInfo(),f_weighted_at_rank=0.8,
                                 bp_match_at_rank=80, query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['superkingdom'] == sk
    print(q_res.summarized_lineage_results['phylum'])
    phy = [SummarizedGatherResult(rank='phylum', fraction=0.1, lineage=RankLineageInfo(lineage_str="a;b"),
                                  f_weighted_at_rank=0.2, bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)),
           SummarizedGatherResult(rank='phylum', fraction=0.9, lineage=RankLineageInfo(), f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                  query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['phylum'] == phy
    print(q_res.summarized_lineage_results['class'])
    cl = [SummarizedGatherResult(rank='class', fraction=0.1,lineage=RankLineageInfo(lineage_str="a;b;c"),
                                  f_weighted_at_rank=0.2, bp_match_at_rank=20, query_ani_at_rank=approx(0.928, rel=1e-2)),
          SummarizedGatherResult(rank='class', fraction=0.9, lineage=RankLineageInfo(), f_weighted_at_rank=0.8, bp_match_at_rank=80,
                                 query_ani_at_rank=None)]
    assert q_res.summarized_lineage_results['class'] == cl

    assert q_res.total_f_weighted['phylum'] == approx(0.2)
    assert q_res.total_f_classified['class'] == approx(0.1)
    assert q_res.total_bp_classified['superkingdom'] == 20


def test_QueryTaxResult_build_summarized_result_over100percent():
    "summarize up ranks: different values"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB','f_unique_to_query': 0.95}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    # now summarize up the ranks
    assert len(q_res.raw_taxresults) == 2
    with pytest.raises(ValueError) as exc:
        q_res.build_summarized_result()
    print(str(exc))
    assert "Summarized fraction is > 100% of the query! This should not be possible" in str(exc)


def test_build_summarized_result_rank_fail_not_available_resummarize():
    "build classification result"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.summarize_up_ranks('superkingdom')
    with pytest.raises(ValueError) as exc:
        q_res.build_summarized_result(single_rank='order')
    print(str(exc))
    assert "Error: rank 'order' not in summarized rank(s), superkingdom" in str(exc)


def test_aggregate_by_lineage_at_rank():
    """test aggregate by lineage at rank"""
    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax, gB_tax])
    # make gather results
    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.5,'f_unique_to_query': 0.4,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    summarized, all_queries = aggregate_by_lineage_at_rank([q_res], rank='phylum', by_query=False)
    print(summarized)
    assert summarized == {'a;b': 0.4,
                          'a;c': 0.3,
                          'unclassified': approx(0.3, rel=1e-2)}
    assert all_queries == ['queryA']


def test_aggregate_by_lineage_at_rank_not_available():
    """test aggregate by lineage at rank"""
    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax, gB_tax])
    # make gather results
    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.5,'f_unique_to_query': 0.4,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    with pytest.raises(ValueError) as exc:
        aggregate_by_lineage_at_rank([q_res], rank='species', by_query=False)
    print(str(exc))
    assert "Error: rank 'species' not available for aggregation." in str(exc)


def test_aggregate_by_lineage_at_rank_by_query():
    """test two queries, aggregate by lineage at rank by query"""
    # make mini taxonomy
    gA_tax = ("gA", "a;b")
    gB_tax = ("gB", "a;c")
    taxD = make_mini_taxonomy([gA_tax, gB_tax])
    # make gather results
    gather_results = [{'query_name': 'queryA', 'name': 'gA', 'f_unique_weighted': 0.2,'f_unique_to_query': 0.2,'unique_intersect_bp': 50}, 
                      {'query_name': 'queryA', "name": 'gB', 'f_unique_weighted': 0.3,'f_unique_to_query': 0.3,'unique_intersect_bp': 30},
                      {'query_name': 'queryB', "name": 'gB', 'f_unique_weighted': 0.4,'f_unique_to_query': 0.4,'unique_intersect_bp': 30}]
    gres = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, summarize=True)
    # check by query
    summarized, all_queries = aggregate_by_lineage_at_rank(gres.values(), rank='superkingdom', by_query=True)
    print(summarized)
    assert summarized == {"a": {'queryA': 0.5, 'queryB': 0.4},
                          "unclassified": {'queryA': 0.5, 'queryB': 0.6}}
    #assert summarized == {'a': {'queryA': approx(0.1, rel=1e-2), 'queryB': 0.7}}
    assert all_queries == ['queryA', 'queryB']
    summarized, all_queries = aggregate_by_lineage_at_rank(gres.values(), rank='phylum', by_query=True)
    print(summarized)
    assert summarized == {'a;c': {'queryA': 0.3, 'queryB': 0.4}, 
                          'a;b': {'queryA': 0.2}, 
                          "unclassified": {'queryA': 0.5, 'queryB': 0.6}}
    

def test_build_classification_result_containment_threshold_fail():
    "classification result: improper containment threshold"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(containment_threshold=1.2)
    print(str(exc))
    assert "Containment threshold must be between 0 and 1 (input value: 1.2)." in str(exc)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(containment_threshold=-.1)
    print(str(exc))
    assert "Containment threshold must be between 0 and 1 (input value: -0.1)." in str(exc)


def test_build_classification_result_containment_threshold():
    "basic functionality: build classification result using containment threshold"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)

    q_res.build_classification_result(containment_threshold=0.1)
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)

    q_res.build_classification_result(containment_threshold=0.2)
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'phylum'
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b")
    assert q_res.classification_result.f_weighted_at_rank == 0.4
    assert q_res.classification_result.fraction == 0.2
    assert q_res.classification_result.bp_match_at_rank == 40
    assert q_res.classification_result.query_ani_at_rank == approx(0.95, rel=1e-2)

    q_res.build_classification_result(containment_threshold=1.0)
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'below_threshold'
    assert q_res.classification_result.rank == 'superkingdom'
    assert q_res.classification_result.fraction == 0.2
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a")
    assert q_res.classification_result.f_weighted_at_rank == 0.4
    assert q_res.classification_result.bp_match_at_rank == 40
    assert q_res.classification_result.query_ani_at_rank == approx(0.95, rel=1e-2)


def test_build_classification_result_ani_threshold():
    "basic functionality: build classification result"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)

    q_res.build_classification_result(ani_threshold=.92)
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)

    q_res.build_classification_result(ani_threshold=0.94) # should classify at phylum
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'phylum'
    assert q_res.classification_result.fraction == 0.2
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b")
    assert q_res.classification_result.f_weighted_at_rank == 0.4
    assert q_res.classification_result.bp_match_at_rank == 40
    assert q_res.classification_result.query_ani_at_rank == approx(0.95, rel=1e-2)

    # superk result, but doesn't meet ANI threshold
    q_res.build_classification_result(ani_threshold=0.96)
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'below_threshold'
    assert q_res.classification_result.rank == 'superkingdom'
    assert q_res.classification_result.fraction == 0.2
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a")
    assert q_res.classification_result.f_weighted_at_rank == 0.4
    assert q_res.classification_result.bp_match_at_rank == 40
    assert q_res.classification_result.query_ani_at_rank == approx(0.95, rel=1e-2)


def test_build_classification_result_ani_threshold_fail():
    "classification result: improper ANI threshold"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(ani_threshold=1.2)
    print(str(exc))
    assert "ANI threshold must be between 0 and 1 (input value: 1.2)." in str(exc)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(ani_threshold=-.1)
    print(str(exc))
    assert "ANI threshold must be between 0 and 1 (input value: -0.1)." in str(exc)


def test_build_classification_result_rank_fail_not_filled():
    "classification result: rank not available (wasn't filled in tax lineage matches)"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(rank='order')
    print(str(exc))
    assert "Error: rank 'order' was not available for any matching lineages." in str(exc)


def test_build_classification_result_rank_fail_not_available_resummarize():
    "classification result: rank not available (wasn't summarized)"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.summarize_up_ranks('superkingdom')
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(rank='order')
    print(str(exc))
    assert "Error: rank 'order' not in summarized rank(s), superkingdom" in str(exc)


def test_build_classification_result_rank_fail_not_available():
    "classification result: rank not available"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    with pytest.raises(ValueError) as exc:
        q_res.build_classification_result(rank='NotARank')
    print(str(exc))
    assert "Error: rank 'NotARank' not in available ranks (strain, species, genus, family, order, class, phylum, superkingdom)" in str(exc)


def test_build_classification_result_rank_containment_threshold():
    "classification result - rank and containment threshold (default)"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)

    q_res.build_classification_result(rank='class')
    print("classif: ", q_res.classification_result)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)

    q_res.build_classification_result(rank='class', containment_threshold=0.4)
    assert q_res.classification_result.status == 'below_threshold'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)


def test_build_classification_result_rank_ani_threshold():
    "classification result with rank and ANI threshold"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)

    q_res.build_classification_result(rank='class', ani_threshold=0.92)
    assert q_res.classification_result.status == 'match'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)

    q_res.build_classification_result(rank='class', ani_threshold=0.95)
    assert q_res.classification_result.status == 'below_threshold'
    assert q_res.classification_result.rank == 'class'
    assert q_res.classification_result.fraction == 0.1
    assert q_res.classification_result.lineage == RankLineageInfo(lineage_str="a;b;c")
    assert q_res.classification_result.f_weighted_at_rank == 0.2
    assert q_res.classification_result.bp_match_at_rank == 20
    assert q_res.classification_result.query_ani_at_rank == approx(0.928, rel=1e-2)


def test_krona_classified():
    "basic functionality: build classification result using containment threshold"
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_classification_result()
    assert q_res.krona_classified == None
    q_res.build_classification_result(rank='phylum')#, force_resummarize=True)
    print(q_res.krona_classified)
    assert q_res.krona_classified == (0.2, 'a', 'b')
    assert q_res.krona_unclassified == (0.8, 'unclassified', 'unclassified')
    q_res.build_classification_result(rank='superkingdom')
    print(q_res.krona_classified)
    assert q_res.krona_classified == (0.2, 'a')
    assert q_res.krona_unclassified == (0.8, 'unclassified')
    # make sure this goes back to None if we reclassify without rank
    q_res.build_classification_result()
    assert q_res.krona_classified == None
    assert q_res.krona_unclassified == None
    assert q_res.krona_header == []


def test_make_krona_header_basic():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    phy_header = ["fraction", "superkingdom", "phylum"]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_classification_result(rank='phylum')
    print(q_res.krona_classified)
    print(q_res.krona_header)
    assert q_res.krona_header == phy_header
    hd = q_res.make_krona_header('phylum')
    print("header: ", hd)
    assert hd == phy_header


def test_make_krona_header_basic_1():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    class_header = ["fraction", "superkingdom", "phylum", "class"]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True)
    q_res.build_classification_result(rank='class')
    assert q_res.krona_header == class_header
    hd = q_res.make_krona_header(min_rank='class')
    print("header: ", hd)
    assert hd == class_header


def test_make_krona_header_fail():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    with pytest.raises(ValueError) as exc:
        q_res.make_krona_header("order")
    assert "Rank 'order' not present in summarized ranks." in str(exc.value)
    with pytest.raises(ValueError) as exc:
        q_res.make_krona_header("NotARank")
    assert "Rank 'NotARank' not present in summarized ranks." in str(exc.value)


def test_make_human_summary():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    hs = q_res.make_human_summary(display_rank = "superkingdom")
    print(hs)
    assert hs == [{'rank': 'superkingdom', 'fraction': '0.800', 'lineage': 'unclassified',
                   'f_weighted_at_rank': '60.0%', 'bp_match_at_rank': "60", 'query_ani_at_rank': '-    ',
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn',
                   'total_weighted_hashes': "0"},
                  {'rank': 'superkingdom', 'fraction': '0.200', 'lineage': "a",
                  'f_weighted_at_rank': '40.0%', 'bp_match_at_rank': "40", 'query_ani_at_rank': '94.9%',
                  'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': "0"}]


def test_make_human_summary_2():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    hs = q_res.make_human_summary(display_rank = "phylum")
    print(hs)
    assert hs == [{'rank': 'phylum', 'fraction': '0.800', 'lineage': 'unclassified',
                   'f_weighted_at_rank': '60.0%', 'bp_match_at_rank': "60", 'query_ani_at_rank': '-    ',
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn',
                   'total_weighted_hashes': "0"},
                  {'rank': 'phylum', 'fraction': '0.200', 'lineage': 'a;b',
                  'f_weighted_at_rank': '40.0%', 'bp_match_at_rank': "40", 'query_ani_at_rank': '94.9%',
                  'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': "0"}]


def test_make_human_summary_classification():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, classify=True, classify_rank="superkingdom")
    hs = q_res.make_human_summary(display_rank = "superkingdom", classification=True)
    print(hs)
    assert hs == [{'rank': 'superkingdom', 'fraction': '0.200', 'lineage': 'a',
                  'f_weighted_at_rank': '40.0%', 'bp_match_at_rank': "40",
                  'query_ani_at_rank': '94.9%', 'status': 'match', 'query_name': 'q1',
                  'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': "0"}]


def test_make_human_summary_classification_2():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, classify=True, classify_rank="phylum")
    hs = q_res.make_human_summary(display_rank = "phylum", classification=True)
    print(hs)
    assert hs == [{'rank': 'phylum', 'fraction': '0.200', 'lineage': 'a;b',
                   'f_weighted_at_rank': '40.0%', 'bp_match_at_rank': "40",
                   'query_ani_at_rank': '94.9%', 'status': 'match',
                   'query_name': 'q1', 'query_md5': 'md5',
                   'query_filename': 'query_fn', 'total_weighted_hashes': "0"}]


def test_make_full_summary():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    header, fs = q_res.make_full_summary()
    assert header == ['query_name', 'rank', 'fraction', 'lineage', 'query_md5', 'query_filename', 
                   'f_weighted_at_rank', 'bp_match_at_rank', 'query_ani_at_rank', 'total_weighted_hashes']
    print(fs)
    assert fs == [{'rank': 'superkingdom', 'fraction': '0.2', 'lineage': 'a', 'f_weighted_at_rank': '0.4',
                   'bp_match_at_rank': '40', 'query_ani_at_rank': approx(0.949,rel=1e-3), 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'superkingdom', 'fraction': '0.8', 'lineage': 'unclassified', 'f_weighted_at_rank':
                   '0.6', 'bp_match_at_rank': '60', 'query_ani_at_rank': None,
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn',
                   'total_weighted_hashes': '0'},
                   {'rank': 'phylum', 'fraction': '0.2', 'lineage': 'a;b', 'f_weighted_at_rank': '0.4',
                   'bp_match_at_rank': '40', 'query_ani_at_rank': approx(0.949,rel=1e-3), 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'phylum', 'fraction': '0.8', 'lineage': 'unclassified', 'f_weighted_at_rank': '0.6',
                   'bp_match_at_rank': '60', 'query_ani_at_rank': None, 'query_name': 'q1', 'query_md5': 'md5',
                   'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.1', 'lineage': 'a;b;c', 'f_weighted_at_rank': '0.2',
                   'bp_match_at_rank': '20', 'query_ani_at_rank': approx(0.928, rel=1e-3),
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.1', 'lineage': 'a;b;d','f_weighted_at_rank': '0.2',
                   'bp_match_at_rank': '20', 'query_ani_at_rank': approx(0.928, rel=1e-3), 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.8', 'lineage': 'unclassified', 'f_weighted_at_rank': '0.6',
                   'bp_match_at_rank': '60', 'query_ani_at_rank': None, 'query_name': 'q1', 'query_md5': 'md5',
                   'query_filename': 'query_fn', 'total_weighted_hashes': '0'}]
    
    header, fs = q_res.make_full_summary(limit_float=True)
    assert header == ['query_name', 'rank', 'fraction', 'lineage', 'query_md5', 'query_filename', 
                   'f_weighted_at_rank', 'bp_match_at_rank', 'query_ani_at_rank', 'total_weighted_hashes']
    print(fs)
    assert fs == [{'rank': 'superkingdom', 'fraction': '0.200', 'lineage': 'a', 'f_weighted_at_rank': '0.400',
                   'bp_match_at_rank': '40', 'query_ani_at_rank': "0.949", 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                  {'rank': 'superkingdom', 'fraction': '0.800', 'lineage': 'unclassified', 'f_weighted_at_rank':
                   '0.600', 'bp_match_at_rank': '60', 'query_ani_at_rank': None,
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn',
                   'total_weighted_hashes': '0'},
                   {'rank': 'phylum', 'fraction': '0.200', 'lineage': 'a;b', 'f_weighted_at_rank': '0.400',
                   'bp_match_at_rank': '40', 'query_ani_at_rank': "0.949", 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'phylum', 'fraction': '0.800', 'lineage': 'unclassified', 'f_weighted_at_rank': '0.600',
                   'bp_match_at_rank': '60', 'query_ani_at_rank': None, 'query_name': 'q1', 'query_md5': 'md5',
                   'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.100', 'lineage': 'a;b;c', 'f_weighted_at_rank': '0.200',
                   'bp_match_at_rank': '20', 'query_ani_at_rank': "0.928",
                   'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.100', 'lineage': 'a;b;d','f_weighted_at_rank': '0.200',
                   'bp_match_at_rank': '20', 'query_ani_at_rank': "0.928", 'query_name': 'q1',
                   'query_md5': 'md5', 'query_filename': 'query_fn', 'total_weighted_hashes': '0'},
                   {'rank': 'class', 'fraction': '0.800', 'lineage': 'unclassified', 'f_weighted_at_rank': '0.600',
                   'bp_match_at_rank': '60', 'query_ani_at_rank': None, 'query_name': 'q1', 'query_md5': 'md5',
                   'query_filename': 'query_fn', 'total_weighted_hashes': '0'}]


def test_make_full_summary_summarization_fail():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=False)
    with pytest.raises(ValueError) as exc:
        q_res.make_full_summary()
    print(str(exc))
    assert 'not summarized yet' in str(exc)


def test_make_full_summary_classification():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, classify=True)
    header, fs = q_res.make_full_summary(classification=True)
    assert header == ["query_name", "status", "rank", "fraction", "lineage",
                     "query_md5", "query_filename", "f_weighted_at_rank",
                     "bp_match_at_rank", "query_ani_at_rank"]
    print(fs)
    assert fs == [{'rank': 'class', 'fraction': '0.1', 'lineage': 'a;b;c', 'f_weighted_at_rank': '0.2', 
                   'bp_match_at_rank': '20', 'query_ani_at_rank': approx(0.928, rel=1e-3),
                   'status': 'match', 'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn'}]

 
def test_make_full_summary_classification_limit_float():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, classify=True)
    header, fs = q_res.make_full_summary(classification=True, limit_float=True)
    assert header == ["query_name", "status", "rank", "fraction", "lineage",
                     "query_md5", "query_filename", "f_weighted_at_rank",
                     "bp_match_at_rank", "query_ani_at_rank"]
    print(fs)
    assert fs == [{'rank': 'class', 'fraction': '0.100', 'lineage': 'a;b;c', 'f_weighted_at_rank': '0.200', 
                   'bp_match_at_rank': '20', 'query_ani_at_rank': "0.928",
                   'status': 'match', 'query_name': 'q1', 'query_md5': 'md5', 'query_filename': 'query_fn'}]


def test_make_full_summary_classification_fail():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    with pytest.raises(ValueError) as exc:
        q_res.make_full_summary(classification=True)
    print(str(exc))
    assert 'not classified yet' in str(exc)


def test_make_kreport_results():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;c;d;e;f;g")])
    #need to go down to species to check that `num_bp_assigned` is happening correctly
    gather_results = [{"total_weighted_hashes":100}, {"name": 'gB', "total_weighted_hashes":100}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    header, krepD = q_res.make_kreport_results()
    print(krepD)
    assert krepD == [{'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'D', 'sci_name': 'a', 'ncbi_taxid': None},
                    {'num_bp_assigned': '60', 'percent_containment': '60.00', 'num_bp_contained': '60',
                    'sci_name': 'unclassified', 'rank_code': 'U', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'P', 'sci_name': 'b', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'C', 'sci_name': 'c', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'O', 'sci_name': 'd', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'F', 'sci_name': 'e', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'G', 'sci_name': 'f', 'ncbi_taxid': None},
                    {'num_bp_assigned': '20', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'S', 'sci_name': 'g', 'ncbi_taxid': None}]


def test_make_kreport_results_with_taxids():
    taxD = make_mini_taxonomy_with_taxids([("gA", "a;b;c", "1;2;3"), ("gB", "a;b;c;d;e;f;g", "1;2;3;4;5;6;7")])
    print(taxD)
    #need to go down to species to check that `num_bp_assigned` is happening correctly
    gather_results = [{"total_weighted_hashes":100}, {"name": 'gB', "total_weighted_hashes":100}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    header, krepD = q_res.make_kreport_results()
    print(krepD)
    assert krepD == [{'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'D', 'sci_name': 'a', 'ncbi_taxid': '1'},
                    {'num_bp_assigned': '60', 'percent_containment': '60.00', 'num_bp_contained': '60',
                    'sci_name': 'unclassified', 'rank_code': 'U', 'ncbi_taxid': None},
                    {'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'P', 'sci_name': 'b', 'ncbi_taxid': '2'},
                    {'num_bp_assigned': '0', 'percent_containment': '40.00', 'num_bp_contained': '40',
                    'rank_code': 'C', 'sci_name': 'c', 'ncbi_taxid': '3'},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'O', 'sci_name': 'd', 'ncbi_taxid': '4'},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'F', 'sci_name': 'e', 'ncbi_taxid': '5'},
                    {'num_bp_assigned': '0', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'G', 'sci_name': 'f', 'ncbi_taxid': '6'},
                    {'num_bp_assigned': '20', 'percent_containment': '20.00', 'num_bp_contained': '20',
                    'rank_code': 'S', 'sci_name': 'g', 'ncbi_taxid': '7'}]


def test_make_kreport_results_fail():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=False)
    with pytest.raises(ValueError) as exc:
        q_res.make_kreport_results()
    print(str(exc))
    assert 'not summarized yet' in str(exc)


def test_make_kreport_results_fail_pre_v450():
    taxD = make_mini_taxonomy([("gA", "a;b;c"), ("gB", "a;b;d")])
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    with pytest.raises(ValueError) as exc:
        q_res.make_kreport_results()
    print(str(exc))
    assert "cannot produce 'kreport' format from gather results before sourmash v4.5.0" in str(exc)


def test_make_cami_results_with_taxids():
    taxD = make_mini_taxonomy_with_taxids([("gA", "a;b;c", "1;2;3"), ("gB", "a;b;c;d;e;f;g", "1;2;3;4;5;6;7")])
    print(taxD)
    #need to go down to species to check that `num_bp_assigned` is happening correctly
    gather_results = [{"total_weighted_hashes":100}, {"name": 'gB', "total_weighted_hashes":100}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True)
    header, camires = q_res.make_cami_bioboxes()
    print(camires)
    assert camires == [['1', 'superkingdom', '1', 'a', '40.00'],
                       ['2', 'phylum', '1|2', 'a|b', '40.00'],
                       ['3', 'class', '1|2|3', 'a|b|c', '40.00'],
                       ['4', 'order', '1|2|3|4', 'a|b|c|d', '20.00'],
                       ['5', 'family', '1|2|3|4|5', 'a|b|c|d|e', '20.00'],
                       ['6', 'genus', '1|2|3|4|5|6', 'a|b|c|d|e|f', '20.00'],
                       ['7', 'species', '1|2|3|4|5|6|7', 'a|b|c|d|e|f|g', '20.00']]


def test_make_lingroup_results():
    taxD = make_mini_taxonomy([("gA", "1;0;0"), ("gB", "1;0;1"), ("gC", "1;1;0")], LIN=True)
    print(taxD)
    lingroupD = {"1":"lg1", "1;0":'lg2', '1;1': "lg3"}
    print(lingroupD)
    gather_results = [{"total_weighted_hashes":100},
                      {"name": 'gB', "total_weighted_hashes":100},
                      {"name": 'gC', "total_weighted_hashes":100}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True, LIN=True)
    print(q_res.summarized_lineage_results)

    header, lgD = q_res.make_lingroup_results(LINgroupsD = lingroupD)
    print(header)
    assert header == ['name', 'lin', 'percent_containment', 'num_bp_contained']
    # order may change, just check that each lg entry is present in list of results
    lg1 = {'percent_containment': '60.00', 'num_bp_contained': '60',
                    'lin': '1', 'name': 'lg1'}
    lg2 = {'percent_containment': '40.00', 'num_bp_contained': '40',
                    'lin': '1;0', 'name': 'lg2'}
    lg3 = {'percent_containment': '20.00', 'num_bp_contained': '20',
                    'lin': '1;1', 'name': 'lg3'}
    assert lg1 in lgD
    assert lg2 in lgD
    assert lg3 in lgD


def test_make_lingroup_results_fail_pre_v450():
    taxD = make_mini_taxonomy([("gA", "1;0;0"), ("gB", "1;0;1"), ("gC", "1;1;0")], LIN=True)
    gather_results = [{}, {"name": 'gB'}]
    q_res = make_QueryTaxResults(gather_info=gather_results, taxD=taxD, single_query=True, summarize=True, LIN=True)
    lingroupD = {"1":"lg1", "1;0":'lg2', '1;1': "lg3"}
    with pytest.raises(ValueError) as exc:
        q_res.make_lingroup_results(lingroupD)
    print(str(exc))
    assert "cannot produce 'lingroup' format from gather results before sourmash v4.5.0" in str(exc)


def test_read_lingroups(runtmp):
    lg_file = runtmp.output("test.lg.csv")
    with open(lg_file, 'w') as out:
        out.write('lin,name\n')
        out.write('1,lg1\n')
        out.write('1;0,lg2\n')
        out.write('1;1,lg3\n')
    lgD = read_lingroups(lg_file)

    assert lgD == {"1":"lg1", "1;0":'lg2', '1;1': "lg3"}

def test_read_lingroups_empty_file(runtmp):
    lg_file = runtmp.output("test.lg.csv")
    with open(lg_file, 'w') as out:
        out.write("")
    with pytest.raises(ValueError) as exc:
        read_lingroups(lg_file)
    print(str(exc))
    assert f"Cannot read lingroups from '{lg_file}'. Is file empty?" in str(exc)


def test_read_lingroups_only_header(runtmp):
    lg_file = runtmp.output("test.lg.csv")
    with open(lg_file, 'w') as out:
        out.write('lin,name\n')
    with pytest.raises(ValueError) as exc:
        read_lingroups(lg_file)
    print(str(exc))
    assert f"No lingroups loaded from {lg_file}" in str(exc)


def test_read_lingroups_bad_header(runtmp):
    lg_file = runtmp.output("test.lg.csv")
    with open(lg_file, 'w') as out:
        out.write('LINgroup_pfx,LINgroup_nm\n')
    with pytest.raises(ValueError) as exc:
        read_lingroups(lg_file)
    print(str(exc))
    assert f"'{lg_file}' must contain the following columns: 'name', 'lin'." in str(exc)


def test_LineageTree_init():
    x = "a;b"
    lin1 = RankLineageInfo(lineage_str=x)
    print(lin1)
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a'):
                         { LineagePair('phylum', 'b') : {}} }

def test_LineageTree_init_mult():
    x = "a;b"
    y = "a;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    print(lin1)
    from sourmash.tax.tax_utils import LineageTree
    tree = LineageTree([lin1, lin2])
    assert tree.tree == {LineagePair(rank='superkingdom', name='a', taxid=None): 
                          {LineagePair(rank='phylum', name='b', taxid=None): {},
                           LineagePair(rank='phylum', name='c', taxid=None): {}}}


def test_LineageTree_init_and_add_lineage():
    x = "a;b"
    y = "a;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    print(lin1)
    from sourmash.tax.tax_utils import LineageTree
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a'):
                         { LineagePair('phylum', 'b') : {}} }
    tree.add_lineage(lin2)
    assert tree.tree == {LineagePair(rank='superkingdom', name='a', taxid=None): 
                          {LineagePair(rank='phylum', name='b', taxid=None): {},
                           LineagePair(rank='phylum', name='c', taxid=None): {}}}


def test_LineageTree_init_and_add_lineages():
    x = "a;b"
    y = "a;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    print(lin1)
    from sourmash.tax.tax_utils import LineageTree
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a'):
                         { LineagePair('phylum', 'b') : {}} }
    tree.add_lineages([lin2])
    assert tree.tree == {LineagePair(rank='superkingdom', name='a', taxid=None): 
                          {LineagePair(rank='phylum', name='b', taxid=None): {},
                           LineagePair(rank='phylum', name='c', taxid=None): {}}}


def test_build_tree_RankLineageInfo():
    x = "a;b"
    lin1 = RankLineageInfo(lineage_str=x)
    print(lin1)
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a'):
                         { LineagePair('phylum', 'b') : {}} }


def test_build_tree_LINLineageInfo():
    x = "0;3"
    lin1 = LINLineageInfo(lineage_str=x)
    print(lin1)
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('0', '0'):
                         { LineagePair('1', '3') : {}} }


def test_build_tree_2():
    x = "a;b"
    y = "a;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    print(lin1)
    print(lin2)
    tree = LineageTree([lin1,lin2])

    assert tree.tree == { LineagePair('superkingdom', 'a'): { LineagePair('phylum', 'b') : {},
                                           LineagePair('phylum', 'c') : {}} }


def test_build_tree_2_LineagePairs():
    # build tree from LineagePairs
    tree = LineageTree([[LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b')],
                       [LineagePair('superkingdom', 'a'), LineagePair('phylum', 'c')],
                      ])

    assert tree.tree == { LineagePair('superkingdom', 'a'): { LineagePair('phylum', 'b') : {},
                                           LineagePair('phylum', 'c') : {}} }


def test_build_tree_3():
    # empty phylum name
    x='a;'
    lin1 = RankLineageInfo(lineage_str=x)
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a'): {} }


def test_build_tree_3_LineagePairs():
    # empty phylum name: LineagePair input
    lin1 = (LineagePair('superkingdom', "a", '3'),
            LineagePair('phylum', '', ''),)
    tree = LineageTree([lin1])
    assert tree.tree == { LineagePair('superkingdom', 'a', '3'): {} }


def test_build_tree_5():
    with pytest.raises(ValueError):
        tree = LineageTree([])


def test_build_tree_5b():
    with pytest.raises(ValueError):
        tree = LineageTree("")


def test_build_tree_iterable():
    with pytest.raises(ValueError) as exc:
        tree = LineageTree(RankLineageInfo())
    assert "Must pass in an iterable containing LineagePair or LineageInfo objects"  in str(exc)


def test_find_lca():
    x='a;b'
    lin1 = RankLineageInfo(lineage_str=x)
    tree = LineageTree([lin1])
    lca = tree.find_lca()

    assert lca == ((LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b'),), 0)


def test_find_lca_LineagePairs():
    tree = LineageTree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2')]])
    lca = tree.find_lca()

    assert lca == ((LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2'),), 0)


def test_find_lca_2():
    x = "a;b"
    y = "a;c"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)

    tree = LineageTree([lin1, lin2])
    lca = tree.find_lca()

    assert lca == ((LineagePair('superkingdom', 'a'),), 2)


def test_find_lca_LIN():
    x = "5;6"
    y = "5;10"
    lin1 = LINLineageInfo(lineage_str=x)
    lin2 = LINLineageInfo(lineage_str=y)

    tree = LineageTree([lin1, lin2])
    lca = tree.find_lca()

    assert lca == ((LineagePair('0', '5'),), 2)
    print(lca)


def test_find_lca_2_LineagePairs():
    tree = LineageTree([[LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2a')],
                       [LineagePair('rank1', 'name1'), LineagePair('rank2', 'name2b')],
                      ])
    lca = tree.find_lca()

    assert lca == ((LineagePair('rank1', 'name1'),), 2)


def test_find_lca_3():
    lin1 = RankLineageInfo(lineage_str="a;b;c")
    lin2 = RankLineageInfo(lineage_str="a;b")

    tree = LineageTree([lin1, lin2])
    lca, reason = tree.find_lca()
    assert lca == lin1.filled_lineage           # find most specific leaf node
    print(lca)


def test_build_tree_with_initial():
    x = "a;b;c"
    y = "a;b;d"
    z = "a;e"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    lin3 = RankLineageInfo(lineage_str=z)

    tree = LineageTree([lin1, lin2])
    lca = tree.find_lca()

    print(lca)
    assert lca == ((LineagePair(rank='superkingdom', name='a', taxid=None),
                    LineagePair(rank='phylum', name='b', taxid=None)), 2)
    tree.add_lineages([lin3])
    lca2 = tree.find_lca()
    print(lca2)
    assert lca2 == ((LineagePair('superkingdom', 'a'),), 2)


def test_LineageTree_find_ordered_paths():
    x = "a;b;c"
    y = "a;b;d"
    z = "a;e"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    lin3 = RankLineageInfo(lineage_str=z)

    tree = LineageTree([lin1, lin2, lin3])
    paths = tree.ordered_paths()

    print(paths)
    assert paths == [(LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='e', taxid=None)),
                     (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='b', taxid=None),
                        LineagePair(rank='class', name='c', taxid=None)),
                     (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='b', taxid=None),
                        LineagePair(rank='class', name='d', taxid=None))]


def test_LineageTree_find_ordered_paths_include_internal():
    x = "a;b;c"
    y = "a;b;d"
    z = "a;e"
    lin1 = RankLineageInfo(lineage_str=x)
    lin2 = RankLineageInfo(lineage_str=y)
    lin3 = RankLineageInfo(lineage_str=z)

    tree = LineageTree([lin1, lin2, lin3])
    paths = tree.ordered_paths(include_internal=True)

    print(paths)

    assert paths == [(LineagePair(rank='superkingdom', name='a', taxid=None),),
                     (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='e', taxid=None)),
                     (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='b', taxid=None)),
                     (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='b', taxid=None),
                        LineagePair(rank='class', name='c', taxid=None)),
                      (LineagePair(rank='superkingdom', name='a', taxid=None),
                        LineagePair(rank='phylum', name='b', taxid=None),
                        LineagePair(rank='class', name='d', taxid=None))]
