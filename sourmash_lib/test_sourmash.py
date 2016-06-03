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
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

    
        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('plot-comparison.py',
                                           ['cmp'],
                                           in_directory=location)


def test_do_plot_comparison_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

    
        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('plot-comparison.py',
                                           ['cmp', '--pdf'],
                                           in_directory=location)


def test_do_plot_comparison_3():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

    
        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('plot-comparison.py',
                                           ['cmp', '--labels'],
                                           in_directory=location)
