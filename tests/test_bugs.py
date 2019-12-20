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


@utils.in_tempdir
def test_bug_803(c):
    # can we do a 'sourmash search' on an LCA database and a query with abundance?
    query = utils.get_test_data('47.abunds.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    c.run_sourmash('search', query, lca_db)
    print(c)
    assert 'NC_009665.1 Shewanella baltica OS185, complete genome' in str(c)
