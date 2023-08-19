import glob
import os

import numpy as np
import pytest

import sourmash
from sourmash.compare import (compare_all_pairs, compare_parallel,
                              compare_serial, compare_serial_containment,
                              compare_serial_max_containment, compare_serial_avg_containment)
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
    sigfiles = ["2.fa.sig", "2+63.fa.sig", "47.fa.sig", "63.fa.sig"]
    filenames = [utils.get_test_data(c) for c in sigfiles]
    sigs = []
    for filename in filenames:
        these_sigs = sourmash.load_file_as_signatures(filename, ksize=31)
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
    jANI = compare_serial(scaled_siglist, ignore_abundance, downsample=False, return_ani=True)
    print(jANI)
    
    true_jaccard_ANI = np.array(
           [[1., 0.978, 0., 0.],
           [0.978, 1., 0.96973012, 0.99262776],
           [0., 0.96973012, 1., 0.97697011],
           [0., 0.99262776, 0.97697011, 1.]])

    np.testing.assert_array_almost_equal(jANI, true_jaccard_ANI, decimal=3)


def test_compare_parallel_jaccardANI(scaled_siglist, ignore_abundance):
    jANI = compare_parallel(scaled_siglist, ignore_abundance, downsample=False, n_jobs=2, return_ani=True)

    true_jaccard_ANI = np.array(
           [[1., 0.978, 0., 0.],
           [0.978, 1., 0.96973012, 0.99262776],
           [0., 0.96973012, 1., 0.97697011],
           [0., 0.99262776, 0.97697011, 1.]])

    np.testing.assert_array_almost_equal(jANI, true_jaccard_ANI, decimal=3)


def test_compare_all_pairs_jaccardANI(scaled_siglist, ignore_abundance):
    similarities_parallel = compare_all_pairs(scaled_siglist, ignore_abundance, downsample=False, n_jobs=2, return_ani=True)
    similarities_serial = compare_serial(scaled_siglist, ignore_abundance, downsample=False, return_ani=True)
    np.testing.assert_array_equal(similarities_parallel, similarities_serial)


def test_compare_serial_containmentANI(scaled_siglist):
    containment_ANI = compare_serial_containment(scaled_siglist, return_ani=True)
    print(containment_ANI)

    true_containment_ANI = np.array(
        [[1, 0.966, 0., 0.],
        [1, 1., 0.97715525, 1.],
        [0., 0.96377054, 1., 0.97678608],
        [0., 0.98667513, 0.97715525, 1.]])

    np.testing.assert_array_almost_equal(containment_ANI, true_containment_ANI, decimal=3)


def test_compare_serial_maxcontainmentANI(scaled_siglist):

    # check max_containment ANI
    max_containment_ANI = compare_serial_max_containment(scaled_siglist, return_ani=True)
    print(max_containment_ANI)

    true_max_containment_ANI = np.array(
        [[1., 1., 0., 0.],
        [1., 1., 0.97715525, 1.],
        [0., 0.97715525, 1., 0.97715525],
        [0., 1., 0.97715525, 1.]])

    np.testing.assert_array_almost_equal(max_containment_ANI, true_max_containment_ANI, decimal=3)


def test_compare_serial_avg_containmentANI(scaled_siglist):

    # check avg_containment ANI
    avg_containment_ANI = compare_serial_avg_containment(scaled_siglist, return_ani=True)
    print(avg_containment_ANI)

    true_avg_containment_ANI = np.array(
        [[1., 0.983, 0., 0.],
        [0.983, 1., 0.97046289, 0.99333757],
        [0., 0.97046289, 1., 0.97697067],
        [0., 0.99333757, 0.97697067, 1.]])

    np.testing.assert_array_almost_equal(avg_containment_ANI, true_avg_containment_ANI, decimal=3)
