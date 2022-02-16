"""
Tests for the 'sourmash signature fileinfo' command line.
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

## command line tests


def test_fileinfo_1_sig(runtmp):
    c = runtmp

    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')

    shutil.copyfile(sig47, runtmp.output('sig47.sig'))
    c.run_sourmash('sig', 'fileinfo', 'sig47.sig')

    out = c.last_result.out
    print(c.last_result.out)

    expected_output = """\
path filetype: MultiIndex
location: sig47.sig
is database? no
has manifest? yes
is nonempty? yes
num signatures: 1
5177 total hashes
abundance information available: no
ksizes present: 31
moltypes present: DNA
scaled vals present: 1000
no num sketches present
""".splitlines()
    for line in expected_output:
        assert line.strip() in out
