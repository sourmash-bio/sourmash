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
            "NC_003197.2 Salmonella enterica subsp. enterica serovar T...",
            0.33083219645293316,
        ),
        (
            "NC_006905.1 Salmonella enterica subsp. enterica serovar C...",
            0.3219645293315143,
        ),
        (
            "NC_011080.1 Salmonella enterica subsp. enterica serovar N...",
            0.3199181446111869,
        ),
        (
            "NC_011978.1 Thermotoga neapolitana DSM 4359, complete genome",
            0.12824010914051842,
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
        ("NC_003197.2 Salmonella enterica subsp...", 0.33083219645293316, 1.0, 0),
        ("NC_011978.1 Thermotoga neapolitana DS...", 0.12824010914051842, 1.0, 0),
        (
            "NC_006905.1 Salmonella enterica subsp...",
            0.0791268758526603,
            0.2457627118644068,
            0,
        ),
        (
            "NC_011080.1 Salmonella enterica subsp...",
            0.04433833560709413,
            0.13859275053304904,
            0,
        ),
    ]
