"""
Tests for the 'sourmash database' command line.
"""
import csv
import shutil
import os
import glob

import pytest
import screed

import sourmash_tst_utils as utils
import sourmash
from sourmash.signature import load_signatures
from sourmash.manifest import CollectionManifest
from sourmash_tst_utils import SourmashCommandFailed


def test_database_extract_picklist_name(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    c=runtmp
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('index', '-k', '31', 'zzz.sbt.zip', sig47, sig63)
    assert c.last_result.status == 0

    test_sbt = c.output('zzz.sbt.zip')
    assert os.path.exists(test_sbt)
    name = "NC_009665.1 Shewanella baltica OS185, complete genome"
    md5 = "09a08691ce52952152f0e866a59f6261"
    md5prefix = "09a08691"

    c.sourmash('database', 'extract', test_sbt, '--name', name)

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


def test_database_extract_picklist_md5(runtmp):
    # extract 47 from 47, using a picklist w/full md5
    c=runtmp
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('index', '-k', '31', 'zzz.sbt.zip', sig47, sig63)
    assert c.last_result.status == 0

    test_sbt = c.output('zzz.sbt.zip')
    assert os.path.exists(test_sbt)
    name = "NC_009665.1 Shewanella baltica OS185, complete genome"
    md5 = "09a08691ce52952152f0e866a59f6261"
    md5prefix = "09a08691"

    c.sourmash('database', 'extract', test_sbt, '--md5', md5)

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig

    # select on any of these attributes
    #row = dict(exactName='NC_009665.1 Shewanella baltica OS185, complete genome',
    #           md5full='09a08691ce52952152f0e866a59f6261',
    #           md5short='09a08691ce5295215',
    #           fullIdent='NC_009665.1',
    #           nodotIdent='NC_009665')
