from __future__ import print_function, unicode_literals
from . import sourmash_tst_utils as utils

def test_bug_781():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('protein_781.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['compare',
                                            '--protein', '--no-dna', '--no-dayhoff',
                                            '--no-hp',  '-k', '11', '-o', 'testing',
                                            testdata1],
                                           in_directory=location)
        print(out)
        print(err)
        assert status == 0
