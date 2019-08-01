import glob
import os

import numpy as np
import pytest

import sourmash
from sourmash.compare import compare_all_pairs, _compare_serial, memmap_siglist
from . import sourmash_tst_utils as utils


@pytest.fixture()
def siglist():
    demo_path = utils.get_test_data("demo")
    filenames = sorted(glob.glob(os.path.join(demo_path, "*.sig")))
    sigs = []
    for filename in filenames:
        sigs.extend(sourmash.load_signatures(filename))
    return sigs


@pytest.fixture()
def ignore_abundance(track_abundance):
    return not track_abundance


def test__compare_serial(siglist, ignore_abundance):
    similarities = _compare_serial(siglist, ignore_abundance)

    true_similarities = np.array(
        [[1., 0.356, 0.078, 0.086, 0., 0., 0.],
         [0.356, 1., 0.072, 0.078, 0., 0., 0.],
         [0.078, 0.072, 1., 0.074, 0., 0., 0.],
         [0.086, 0.078, 0.074, 1., 0., 0., 0.],
         [0., 0., 0., 0., 1., 0.382, 0.364],
         [0., 0., 0., 0., 0.382, 1., 0.386],
         [0., 0., 0., 0., 0.364, 0.386, 1.]])

    assert np.array_equal(similarities, true_similarities)


def test_memmap_siglist(siglist):
    memmapped = memmap_siglist(siglist)
    # Assert that the data didn't change as a result of memory-mapping
    np.testing.assert_array_equal(memmapped, siglist)


def test_compare_all_pairs(siglist, ignore_abundance):
    similarities_parallel = compare_all_pairs(siglist, ignore_abundance, n_jobs=2)
    similarities_serial = _compare_serial(siglist, ignore_abundance)
    assert np.array_equal(similarities_parallel, similarities_serial)
