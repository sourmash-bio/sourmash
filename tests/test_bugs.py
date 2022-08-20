import sourmash_tst_utils as utils

@utils.in_tempdir
def test_bug_803(c):
    # can we do a 'sourmash search' on an LCA database and a query with abundance?
    query = utils.get_test_data('track_abund/47.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    c.run_sourmash('search', query, lca_db, '--ignore-abundance')
    print(c)
    assert 'NC_009665.1 Shewanella baltica OS185, complete genome' in str(c)
