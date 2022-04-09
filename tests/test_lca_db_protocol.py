"""
Test the behavior of LCA databases.
"""
import pytest
import sourmash_tst_utils as utils

import sourmash


def build_inmem_lca_db(runtmp):
    # test command-line creation of LCA database with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = runtmp.output('protein.lca.json')

    runtmp.sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '2', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--protein')

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]

    return db2
    

def build_sql_lca_db(runtmp):
    # test command-line creation of LCA database with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = runtmp.output('protein.lca.sqldb')

    runtmp.sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '2', '--split-identifiers', '--require-taxonomy',
                    '--scaled', '100', '-k', '19', '--protein', '-F', 'sql')

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]

    return db2


@pytest.fixture(params=[build_inmem_lca_db,
                        build_sql_lca_db])
def lca_db_obj(request, runtmp):
    build_fn = request.param

    return build_fn(runtmp)


def test_get_lineage_assignments(lca_db_obj):
    lineages = lca_db_obj.get_lineage_assignments(178936042868009693)

    assert len(lineages) == 1
    lineage = lineages[0]

    x = []
    for tup in lineage:
        x.append((tup[0], tup[1]))

    assert x == [('superkingdom', 'd__Archaea'),
                 ('phylum', 'p__Crenarchaeota'),
                 ('class', 'c__Bathyarchaeia'),
                 ('order', 'o__B26-1'),
                 ('family', 'f__B26-1'), 
                 ('genus', 'g__B26-1'),
                 ('species', 's__B26-1 sp001593925'),
                 ('strain', '')]


def test_hashvals(lca_db_obj):
    # test getting individual hashvals
    hashvals = set(lca_db_obj.hashvals)
    assert 178936042868009693 in hashvals


def test_get_identifiers_for_hashval(lca_db_obj):
    # test getting identifiers belonging to individual hashvals
    idents = lca_db_obj.get_identifiers_for_hashval(178936042868009693)
    idents = list(idents)
    assert len(idents) == 1

    ident = idents[0]
    assert ident == 'GCA_001593925'


def test_get_identifiers_for_hashval_2(lca_db_obj):
    # test systematic hashval => identifiers
    all_idents = set()

    for hashval in lca_db_obj.hashvals:
        idents = lca_db_obj.get_identifiers_for_hashval(hashval)
        #idents = list(idents)
        all_idents.update(idents)

    all_idents = list(all_idents)
    print(all_idents)
    assert len(all_idents) == 2

    assert 'GCA_001593925' in all_idents
    assert 'GCA_001593935' in all_idents
