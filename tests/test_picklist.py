"""
Tests for the picklist API.
"""
import pytest
import sourmash

import sourmash_tst_utils as utils
from sourmash.picklist import SignaturePicklist


def test_load_empty_picklist_fail():
    empty = utils.get_test_data('picklist/empty.csv')

    pl = SignaturePicklist('manifest')
    with pytest.raises(ValueError):
        pl.load(empty, 'foo', allow_empty=False)


def test_load_empty_picklist_allow():
    empty = utils.get_test_data('picklist/empty.csv')

    pl = SignaturePicklist('manifest')
    pl.load(empty, 'foo', allow_empty=True)
