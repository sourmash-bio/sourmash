import glob
import os

import numpy as np
import pytest

import sourmash
from sourmash.compare import (compare_all_pairs, compare_parallel,
                              compare_serial, compare_serial_containment, compare_serial_max_containment)
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
def scaled_siglist():
    demo_path = utils.get_test_data("scaled")
    filenames = sorted(glob.glob(os.path.join(demo_path, "*.sig")))
    sigs = []
    for filename in filenames:
        these_sigs = sourmash.load_file_as_signatures(filename)
        scaled_sigs = [s for s in these_sigs if s.minhash.scaled != 0]
        sigs.extend(scaled_sigs)
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


def test_compare_serial_jaccardANI(scaled_siglist, ignore_abundance):
    similarities = compare_serial(scaled_siglist, ignore_abundance, downsample=False, return_ANI=True)

    true_similarities = np.array(
        [[1., 0.942, 0.988, 0.986, 0.],
        [0.942, 1., 0.960, 0., 0.],
        [0.988, 0.960, 1., 0., 0.],
        [0.986, 0., 0., 1., 0.],
        [0., 0., 0., 0., 1.]])

    np.testing.assert_array_almost_equal(similarities, true_similarities, decimal=3)


def test_compare_parallel_jaccardANI(scaled_siglist, ignore_abundance):
    similarities = compare_parallel(scaled_siglist, ignore_abundance, downsample=False, n_jobs=2, return_ANI=True)

    true_containment = np.array(
        [[1., 0.942, 0.988, 0.986, 0.],
        [0.942, 1., 0.960, 0., 0.],
        [0.988, 0.960, 1., 0., 0.],
        [0.986, 0., 0., 1., 0.],
        [0., 0., 0., 0., 1.]])

    np.testing.assert_array_almost_equal(similarities, true_containment, decimal=3)


def test_compare_all_pairs_jaccardANI(scaled_siglist, ignore_abundance):
    similarities_parallel = compare_all_pairs(scaled_siglist, ignore_abundance, downsample=False, n_jobs=2, return_ANI=True)
    similarities_serial = compare_serial(scaled_siglist, ignore_abundance, downsample=False, return_ANI=True)
    np.testing.assert_array_equal(similarities_parallel, similarities_serial)


def test_compare_serial_containmentANI(scaled_siglist, ignore_abundance):
    containment = compare_serial_containment(scaled_siglist, ignore_abundance, return_ANI=True)

    true_containment = np.array(
        [[1., 1., 1., 1., 0.],
        [0.92391599, 1., 0.94383993, 0., 0.],
        [0.97889056, 1., 1., 0., 0.],
        [0.97685474, 0., 0., 1., 0.],
        [0., 0., 0., 0., 1.]])

    np.testing.assert_array_almost_equal(containment, true_containment, decimal=3)

    # check max_containment ANI
    max_containment = compare_serial_max_containment(scaled_siglist, ignore_abundance, return_ANI=True)

    true_max_containment = np.array(
        [[1., 1., 1., 1., 0.],
        [1., 1., 1., 0., 0.],
        [1., 1., 1., 0., 0.],
        [1., 0., 0., 1., 0.],
        [0., 0., 0., 0., 1.,]])

    np.testing.assert_array_almost_equal(max_containment, true_max_containment, decimal=3)
