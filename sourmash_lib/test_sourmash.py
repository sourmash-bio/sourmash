from __future__ import print_function
import os
import glob

from . import sourmash_tst_utils as utils
import matplotlib
matplotlib.use('Agg')

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


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

        status, out, err = utils.runscript('sourmash', ['plot', 'cmp'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


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

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--pdf'],
                                           in_directory=location)
        assert os.path.exists(os.path.join(location, "cmp.dendro.pdf"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.pdf"))


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

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--labels'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


def test_sourmash_search():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '0.958' in out


def test_mash_csv_to_sig():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa.msh.dump')
        testdata2 = utils.get_test_data('short.fa')

        status, out, err = utils.runscript('sourmash', ['import_csv',
                                                        testdata1,
                                                        '-o', 'xxx.sig'],
                                           in_directory=location)
        
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '31',
                                            'short.fa.sig', 'xxx.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches:' in out
