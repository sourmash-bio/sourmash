"""
Test the behavior of LCA databases. New LCA database classes should support
all of this functionality.
"""
import pytest
import sourmash_tst_utils as utils

import sourmash
from sourmash.tax.tax_utils import MultiLineageDB
from sourmash.lca.lca_db import (LCA_Database, load_single_database)


def build_inmem_lca_db(runtmp):
    # test in-memory LCA_Database
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    ss1 = sourmash.load_one_signature(sigfile1)
    ss2 = sourmash.load_one_signature(sigfile2)

    lineages_file = utils.get_test_data('prot/gtdb-subset-lineages.csv')
    lineages = MultiLineageDB.load([lineages_file])

    db = LCA_Database(ksize=19, scaled=100, moltype='protein')

    ident1 = ss1.name.split(' ')[0].split('.')[0]
    assert lineages[ident1]
    db.insert(ss1, ident=ident1, lineage=lineages[ident1])
    ident2 = ss2.name.split(' ')[0].split('.')[0]
    assert lineages[ident2]
    db.insert(ss2, ident=ident2, lineage=lineages[ident2])

    return db


def build_json_lca_db(runtmp):
    # test saved/loaded JSON database
    db = build_inmem_lca_db(runtmp)
    db_out = runtmp.output('protein.lca.json')

    db.save(db_out, format='json')

    x = load_single_database(db_out)
    db_load = x[0]

    return db_load


def build_sql_lca_db(runtmp):
    # test saved/loaded SQL database
    db = build_inmem_lca_db(runtmp)
    db_out = runtmp.output('protein.lca.json')

    db.save(db_out, format='sql')

    x = load_single_database(db_out)
    db_load = x[0]

    return db_load


@pytest.fixture(params=[build_inmem_lca_db,
                        build_json_lca_db,
                        build_sql_lca_db])
def lca_db_obj(request, runtmp):
    build_fn = request.param

    return build_fn(runtmp)


def test_get_lineage_assignments(lca_db_obj):
    # test get_lineage_assignments for a specific hash
    lineages = lca_db_obj.get_lineage_assignments(178936042868009693)

    assert len(lineages) == 1
    lineage = lineages[0]

    x = []
    for tup in lineage:
        if tup[0] != 'strain' or tup[1]: # ignore empty strain
            x.append((tup[0], tup[1]))

    assert x == [('superkingdom', 'd__Archaea'),
                 ('phylum', 'p__Crenarchaeota'),
                 ('class', 'c__Bathyarchaeia'),
                 ('order', 'o__B26-1'),
                 ('family', 'f__B26-1'), 
                 ('genus', 'g__B26-1'),
                 ('species', 's__B26-1 sp001593925'),]


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


def test_downsample_scaled(lca_db_obj):
    # check the downsample_scaled method
    assert lca_db_obj.scaled == 100
    lca_db_obj.downsample_scaled(500)
    assert lca_db_obj.scaled == 500


def test_downsample_scaled_fail(lca_db_obj):
    # check the downsample_scaled method - should fail if lower scaled.
    assert lca_db_obj.scaled == 100

    with pytest.raises(ValueError):
        lca_db_obj.downsample_scaled(50)
