"""
Tests for the 'sourmash signature' command line.
"""
from __future__ import print_function, unicode_literals

from . import sourmash_tst_utils as utils
import sourmash_lib

## command line tests


def test_run_sourmash_signature_cmd():
    status, out, err = utils.runscript('sourmash', ['signature'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_sig_cmd():
    status, out, err = utils.runscript('sourmash', ['sig'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)
