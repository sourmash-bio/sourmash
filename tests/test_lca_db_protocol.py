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
    

@pytest.fixture(params=[build_inmem_lca_db,])
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


def test_hashval_to_ident(lca_db_obj):
    # @CTB: abstract me a bit. hashvals, hashval_to_ident, hashval_to_lineages?
    idxlist = lca_db_obj.hashval_to_idx[178936042868009693]

    assert len(idxlist) == 1
    idx = idxlist[0]

    ident = lca_db_obj._idx_to_ident[idx]
    assert ident == 'GCA_001593925'
