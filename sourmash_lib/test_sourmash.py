from __future__ import print_function
from . import sourmash_tst_utils as utils

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_plot_comparison():
    status, out, err = utils.runscript('plot-comparison.py', [], fail_ok=True)
    assert status != 0
