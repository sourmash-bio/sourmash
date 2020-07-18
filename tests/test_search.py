"""
Tests for the 'search' and 'gather' internal functions
"""
from __future__ import print_function, unicode_literals

import glob

import pytest


from . import sourmash_tst_utils as utils
from sourmash.search import search_databases, gather_databases
from sourmash import sourmash_args


@pytest.fixture
def ksize():
    return 21


@pytest.fixture
def moltype():
    return "DNA"


@pytest.fixture
def scaled():
    return 0


@pytest.fixture
def threshold():
    return 0.08


@pytest.fixture
def threshold_bp():
    return 5e4


@pytest.fixture
def best_only():
    return False


@pytest.fixture
def containment():
    return False


@pytest.fixture
def ignore_abundance():
    return False


@pytest.fixture
def index_sigs():
    testdata_glob = utils.get_test_data("gather/GCF*.sig")
    testdata_sigs = glob.glob(testdata_glob)
    return testdata_sigs


@pytest.fixture
def index_sbt(index_sigs, ksize, moltype, scaled):
    from sourmash.commands import load_matching_signatures_into_tree

    tree = load_matching_signatures_into_tree(
        index_sigs, ksize=ksize, moltype=moltype, scaled=scaled
    )
    return tree


@pytest.fixture
def query_sig():
    return utils.get_test_data("gather/combined.sig")


@pytest.fixture
def query(query_sig, ksize, moltype):
    q = sourmash_args.load_query_signature(
        query_sig, ksize=ksize, select_moltype=moltype
    )
    return q


@pytest.fixture
def databases(index_sbt):
    """Output list of databases like sourmash_args.load_dbs_and_sigs"""
    return [(index_sbt, "gcf_all", "SBT")]


def test_search_databases(
    query, databases, threshold, containment, best_only, ignore_abundance
):
    results = search_databases(
        query,
        databases,
        threshold,
        containment,
        best_only,
        ignore_abundance,
        unload_data=False,
    )
    results_displayable = [
        (sr.match._display_name(60), sr.similarity) for sr in results
    ]
    assert results_displayable == [
        (
            "NC_003198.1 Salmonella enterica subsp. enterica serovar T...",
            0.3321964529331514,
        ),
        (
            "NC_003197.2 Salmonella enterica subsp. enterica serovar T...",
            0.33083219645293316,
        ),
        (
            "NC_006905.1 Salmonella enterica subsp. enterica serovar C...",
            0.3219645293315143,
        ),
        (
            "NC_011294.1 Salmonella enterica subsp. enterica serovar E...",
            0.32128240109140516,
        ),
        (
            "NC_011080.1 Salmonella enterica subsp. enterica serovar N...",
            0.3199181446111869,
        ),
        (
            "NC_011274.1 Salmonella enterica subsp. enterica serovar G...",
            0.3117326057298772,
        ),
        (
            "NC_004631.1 Salmonella enterica subsp. enterica serovar T...",
            0.3035470668485675,
        ),
        (
            "NC_006511.1 Salmonella enterica subsp. enterica serovar P...",
            0.291268758526603,
        ),
        (
            "NC_000853.1 Thermotoga maritima MSB8 chromosome, complete...",
            0.13096862210095497,
        ),
        (
            "NC_009486.1 Thermotoga petrophila RKU-1, complete genome",
            0.1296043656207367,
        ),
        (
            "NC_011978.1 Thermotoga neapolitana DSM 4359, complete genome",
            0.12824010914051842,
        ),
        (
            "NC_002163.1 Campylobacter jejuni subsp. jejuni NCTC 11168...",
            0.10709413369713507,
        ),
    ]


def test_gather_databases(query, databases, threshold_bp, ignore_abundance):
    results = list(gather_databases(query, databases, threshold_bp, ignore_abundance,))
    assert True
    results_displayable = [
        (
            sr[0].match._display_name(40),
            sr[0].f_unique_weighted,
            sr[0].f_match,
            sr[0].average_abund,
        )
        for sr in results
    ]
    assert results_displayable == [
        ("NC_003198.1 Salmonella enterica subsp...", 0.3321964529331514, 1.0, 0),
        ("NC_000853.1 Thermotoga maritima MSB8 ...", 0.13096862210095497, 1.0, 0),
        (
            "NC_011978.1 Thermotoga neapolitana DS...",
            0.11527967257844475,
            0.898936170212766,
            0,
        ),
        ("NC_002163.1 Campylobacter jejuni subs...", 0.10709413369713507, 1.0, 0),
        (
            "NC_003197.2 Salmonella enterica subsp...",
            0.10368349249658936,
            0.3134020618556701,
            0,
        ),
        (
            "NC_009486.1 Thermotoga petrophila RKU...",
            0.06275579809004093,
            0.4842105263157895,
            0,
        ),
        (
            "NC_006905.1 Salmonella enterica subsp...",
            0.05184174624829468,
            0.16101694915254236,
            0,
        ),
        (
            "NC_011080.1 Salmonella enterica subsp...",
            0.04024556616643929,
            0.1257995735607676,
            0,
        ),
        (
            "NC_011274.1 Salmonella enterica subsp...",
            0.0286493860845839,
            0.09190371991247265,
            0,
        ),
        (
            "NC_006511.1 Salmonella enterica subsp...",
            0.021145975443383355,
            0.07259953161592506,
            0,
        ),
        (
            "NC_011294.1 Salmonella enterica subsp...",
            0.0047748976807639835,
            0.014861995753715499,
            0,
        ),
    ]
