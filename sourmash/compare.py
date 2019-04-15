from itertools import combinations
import os
import tempfile

from joblib import Parallel, delayed, load, dump
from scipy.spatial.distance import squareform

import numpy as np


def _compare_serial(siglist, ignore_abundance):
    n = len(siglist)

    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = combinations(range(n), 2)

    similarities = np.ones((n, n))

    for i, j in iterator:
        sig1 = siglist[i]
        similarity = sig1.similarity(siglist[j], ignore_abundance)
        similarities[i][j] = similarity
        similarities[j][i] = similarity

    return similarities


def memmap_siglist(siglist):
    """Write a memory-mapped array of signatures"""
    temp_folder = tempfile.mkdtemp()
    filename = os.path.join(temp_folder, 'siglist.mmap')
    if os.path.exists(filename): os.unlink(filename)
    _ = dump(siglist, filename)
    large_memmap = load(filename, mmap_mode='r+')
    return large_memmap


def compare_all_pairs(siglist, ignore_abundance, n_jobs=None):
    n = len(siglist)

    if n_jobs is None or n_jobs == 1:
        similarities = _compare_serial(siglist, ignore_abundance)
    else:
        # Create a memory-mapped array
        memmapped = memmap_siglist(siglist)
        sig_iterator = combinations(memmapped, 2)

        # This creates a condensed distance matrix
        condensed = Parallel(n_jobs=n_jobs, require='sharedmem',
                             backend='threading')(
            delayed(sig1.similarity)(sig2, ignore_abundance) for sig1, sig2 in sig_iterator)
        similarities = squareform(condensed)

        # 'squareform' was made for *distance* matrices not *similarity*
        # so need to replace diagonal values with 1.
        # np.fill_digonal modifies 'similarities' in-place
        np.fill_diagonal(similarities, 1)
    return similarities
