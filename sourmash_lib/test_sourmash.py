from __future__ import print_function
from . import sourmash_tst_utils as utils
import matplotlib
matplotlib.use('Agg')

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_plot_comparison():
    status, out, err = utils.runscript('plot-comparison.py', [], fail_ok=True)
    assert status != 0


def test_do_plot_comparison():
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    status, out, err = utils.runscript('sourmash',
                                       ['compute', testdata1, testdata2])

    
    status, out, err = utils.runscript('sourmash',
                                       ['compare', 'short.fa.sig',
                                        'short2.fa.sig', '-o', 'cmp'])

    status, out, err = utils.runscript('plot-comparison.py', ['cmp'])


def test_do_plot_comparison():
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    status, out, err = utils.runscript('sourmash',
                                       ['compute', testdata1, testdata2])

    
    status, out, err = utils.runscript('sourmash',
                                       ['compare', 'short.fa.sig',
                                        'short2.fa.sig', '-o', 'cmp'])

    status, out, err = utils.runscript('plot-comparison.py',
                                       ['cmp', '--pdf'])


def test_do_plot_comparison_3():
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    status, out, err = utils.runscript('sourmash',
                                       ['compute', testdata1, testdata2])

    
    status, out, err = utils.runscript('sourmash',
                                       ['compare', 'short.fa.sig',
                                        'short2.fa.sig', '-o', 'cmp'])

    status, out, err = utils.runscript('plot-comparison.py',
                                       ['cmp', '--labels'])

