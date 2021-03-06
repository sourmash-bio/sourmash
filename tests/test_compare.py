import glob
import os

import numpy as np
import pytest

import sourmash
from sourmash.compare import (compare_all_pairs, compare_parallel,
                              compare_serial)
import sourmash_tst_utils as utils


@pytest.fixture()
def siglist():
    demo_path = utils.get_test_data("demo")
    filenames = sorted(glob.glob(os.path.join(demo_path, "*.sig")))
    sigs = []
    for filename in filenames:
        sigs.extend(sourmash.load_file_as_signatures(filename))
    return sigs


@pytest.fixture()
def ignore_abundance(track_abundance):
    return not track_abundance


def test_compare_serial(siglist, ignore_abundance):
    similarities = compare_serial(siglist, ignore_abundance, downsample=False)

    true_similarities = np.array(
        [[1., 0.356, 0.078, 0.086, 0., 0., 0.],
         [0.356, 1., 0.072, 0.078, 0., 0., 0.],
         [0.078, 0.072, 1., 0.074, 0., 0., 0.],
         [0.086, 0.078, 0.074, 1., 0., 0., 0.],
         [0., 0., 0., 0., 1., 0.382, 0.364],
         [0., 0., 0., 0., 0.382, 1., 0.386],
         [0., 0., 0., 0., 0.364, 0.386, 1.]])

    np.testing.assert_array_equal(similarities, true_similarities)


def test_compare_parallel(siglist, ignore_abundance):
    similarities = compare_parallel(siglist, ignore_abundance, downsample=False, n_jobs=2)

    true_similarities = np.array(
        [[1., 0.356, 0.078, 0.086, 0., 0., 0.],
         [0.356, 1., 0.072, 0.078, 0., 0., 0.],
         [0.078, 0.072, 1., 0.074, 0., 0., 0.],
         [0.086, 0.078, 0.074, 1., 0., 0., 0.],
         [0., 0., 0., 0., 1., 0.382, 0.364],
         [0., 0., 0., 0., 0.382, 1., 0.386],
         [0., 0., 0., 0., 0.364, 0.386, 1.]])

    np.testing.assert_array_equal(similarities, true_similarities)


def test_compare_all_pairs(siglist, ignore_abundance):
    similarities_parallel = compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=2)
    similarities_serial = compare_serial(siglist, ignore_abundance, downsample=False)
    np.testing.assert_array_equal(similarities_parallel, similarities_serial)
