from tempfile import NamedTemporaryFile

import pytest

from sourmash.nodegraph import Nodegraph

from . import sourmash_tst_utils as utils


def test_nodegraph_to_khmer_basic():
    pytest.importorskip('khmer')

    ng_file = utils.get_test_data('.sbt.v3/internal.0')

    sourmash_ng = Nodegraph.load(ng_file)
    khmer_sm_ng = sourmash_ng.to_khmer_nodegraph()

    assert sourmash_ng.ksize == khmer_sm_ng.ksize()


def test_nodegraph_khmer_compare():
    khmer = pytest.importorskip('khmer')

    khmer_ng = khmer.Nodegraph(3, 23, 6)
    khmer_ng.count("ACG")
    khmer_ng.count("TTA")
    khmer_ng.count("CGA")

    sm_ng = Nodegraph(3, 23, 6)
    sm_ng.count("ACG")
    sm_ng.count("TTA")
    sm_ng.count("CGA")

    assert sm_ng.ksize == khmer_ng.ksize()
    assert sm_ng.n_tables == len(khmer_ng.hashsizes())
    assert sm_ng.tablesize == sum(khmer_ng.hashsizes())
    assert sm_ng.get("ACG")
    assert sm_ng.get("TTA")
    assert sm_ng.get("CGA")

    assert khmer_ng.get("ACG")
    assert khmer_ng.get("TTA")
    assert khmer_ng.get("CGA")


def test_nodegraph_same_file():
    khmer = pytest.importorskip('khmer')

    ng_file = utils.get_test_data('.sbt.v3/internal.0')
    with open(ng_file, 'rb') as f:
        ng_data = f.read()

    sourmash_ng = Nodegraph.load(ng_file)
    khmer_sm_ng = sourmash_ng.to_khmer_nodegraph()

    khmer_ng = khmer.load_nodegraph(ng_file)

    with NamedTemporaryFile() as f1, NamedTemporaryFile() as f2, NamedTemporaryFile() as f3:
        sourmash_ng.save(f1.name)
        khmer_sm_ng.save(f2.name)
        khmer_ng.save(f3.name)

        f1.seek(0)
        sm_data = f1.read()

        f2.seek(0)
        kh_sm_data = f2.read()

        f3.seek(0)
        kh_data = f3.read()

        assert sm_data == kh_data
        assert sm_data == kh_sm_data

        assert ng_data == sm_data
        assert ng_data == kh_data
        assert ng_data == kh_sm_data
