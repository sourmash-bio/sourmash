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

#from sourmash.lca import lca_utils
#from sourmash.lca.lca_utils import LineagePair

## api tests
# def test_api_xxx


## command line tests
def test_run_sourmash_tax():
    status, out, err = utils.runscript('sourmash', ['tax'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)







## some test ideas to start with -- see test_lca.py for add'l ideas

def test_summarize_empty_gather_results():
    pass
def test_summarize_bad_gather_results():
    pass
def test_summarize_empty_lineage_input():
    pass
def test_summarize_bad_lineage_input():
    pass
def test_summarize_bad_rank():
    pass

def test_classify_empty_gather_results():
    pass
def test_classify_bad_gather_results():
    pass
def test_classify_empty_lineage_input():
    pass
def test_classify_bad_lineage_input():
    pass
def test_single_classify_empty():
    pass
def test_mult_classify_empty():
    pass
