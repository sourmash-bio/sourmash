from tempfile import NamedTemporaryFile

import pytest

from sourmash.nodegraph import Nodegraph

from . import sourmash_tst_utils as utils


def test_nodegraph_to_khmer_basic():
    khmer = pytest.importorskip('khmer')

    ng_file = utils.get_test_data('.sbt.v3/internal.0')

    sourmash_ng = Nodegraph.load(ng_file)
    khmer_sm_ng = sourmash_ng.to_khmer_nodegraph()

    assert sourmash_ng.ksize == khmer_sm_ng.ksize()


def test_nodegraph_same_file():
    khmer = pytest.importorskip('khmer')

    ng_file = utils.get_test_data('.sbt.v3/internal.0')
    with open(ng_file, 'rb') as f:
        ng_data = f.read()

    sourmash_ng = Nodegraph.load(ng_file)
    khmer_sm_ng = sourmash_ng.to_khmer_nodegraph()

    khmer_ng = khmer.load_nodegraph(ng_file)

    with NamedTemporaryFile() as f1, NamedTemporaryFile() as f2:
        sourmash_ng.save(f1.name)
        khmer_ng.save(f2.name)

        f1.seek(0)
        sm_data = f1.read()

        f2.seek(0)
        kh_data = f2.read()

        assert sm_data == kh_data
        assert ng_data == sm_data
        assert ng_data == kh_data
