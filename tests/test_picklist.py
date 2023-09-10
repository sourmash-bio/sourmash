"""
Tests for the picklist API.
"""
import pytest
import sourmash

import copy

import sourmash_tst_utils as utils
from sourmash import picklist
from sourmash.picklist import SignaturePicklist
from sourmash.index import LinearIndex, MultiIndex
from sourmash.index.sqlite_index import SqliteIndex


def test_load_empty_picklist_fail():
    empty = utils.get_test_data('picklist/empty.csv')

    pl = SignaturePicklist('manifest', pickfile=empty)
    with pytest.raises(ValueError):
        pl.load(allow_empty=False)


def test_load_empty_picklist_allow():
    empty = utils.get_test_data('picklist/empty.csv')

    pl = SignaturePicklist('manifest', pickfile=empty)
    pl.load(allow_empty=True)


def test_dup_md5_picked(runtmp):
    # load a sig, duplicate, and see if a picklist gets the right one
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_signatures(sig47)
    sig = list(ss)[0]

    # save a manifest with one entry
    xl = LinearIndex([sig])
    ml = MultiIndex.load([xl], [None], None)

    print(ml.manifest.rows)
    assert len(ml.manifest) == 1

    mf_csv = runtmp.output('select.csv')
    ml.manifest.write_to_filename(mf_csv)

    # now make an index to select against, with an identical signature
    # (but diff name)
    new_sig = sig.to_mutable()
    new_sig.name = 'foo'
    xl = LinearIndex([sig, new_sig])
    ml2 = MultiIndex.load([xl], [None], None)

    assert len(ml2) == 2

    # create a picklist...
    pl = SignaturePicklist('manifest', pickfile=mf_csv)
    print(pl.load())
    print('loaded:', len(pl.pickset))

    # use in select
    ml3 = ml2.select(picklist=pl)
    print('picked:', len(ml3))

    assert len(pl.pickset) == len(ml3)


def test_dup_md5_picked_mf_to_picklist(runtmp):
    # load a sig, duplicate, and see if a picklist gets the right one
    # uses an in memory picklist
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_signatures(sig47)
    sig = list(ss)[0]

    # save a manifest with one entry
    xl = LinearIndex([sig])
    ml = MultiIndex.load([xl], [None], None)

    print(ml.manifest.rows)
    assert len(ml.manifest) == 1

    pl = ml.manifest.to_picklist()

    # now make an index to select against, with an identical signature
    # (but diff name)
    new_sig = sig.to_mutable()
    new_sig.name = 'foo'
    xl = LinearIndex([sig, new_sig])
    ml2 = MultiIndex.load([xl], [None], None)

    assert len(ml2) == 2

    # use picklist in select
    ml3 = ml2.select(picklist=pl)
    print('picked:', len(ml3))

    assert len(pl.pickset) == len(ml3)


def test_dup_md5_picked_mf_to_picklist_sqlite(runtmp):
    # load a sig, duplicate, and see if a picklist gets the right one
    # use a sqlite db with its own to_picklist behavior.
    sig47 = utils.get_test_data('47.fa.sig')
    ss = sourmash.load_signatures(sig47)
    sig = list(ss)[0]

    # save a manifest with one entry
    xl = SqliteIndex.create(':memory:')
    xl.insert(sig)

    print(xl.manifest.rows)
    assert len(xl.manifest) == 1

    pl = xl.manifest.to_picklist()

    # now make an index to select against, with an identical signature
    # (but diff name)
    new_sig = sig.to_mutable()
    new_sig.name = 'foo'
    xl = LinearIndex([sig, new_sig])
    ml2 = MultiIndex.load([xl], [None], None)

    assert len(ml2) == 2

    # use picklist in select
    ml3 = ml2.select(picklist=pl)
    print('picked:', len(ml3))

    assert len(pl.pickset) == len(ml3)