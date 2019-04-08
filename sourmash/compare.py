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

    values = np.ones((n, n))

    for i, j in iterator:
        sig1 = siglist[i]
        similarity = sig1.similarity(siglist[j], ignore_abundance)
        values[i][j] = similarity
        values[j][i] = similarity

    return values


def _memmap_siglist(siglist):
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
        values = _compare_serial(siglist, ignore_abundance)
    else:
        # Create a memory-mapped array
        memmapped = _memmap_siglist(siglist)
        sig_iterator = combinations(memmapped, 2)

        # This creates a condensed distance matrix
        condensed = Parallel(n_jobs=n_jobs, require='sharedmem',
                             backend='threading')(
            delayed(sig1.similarity)(sig2, ignore_abundance) for sig1, sig2 in sig_iterator)
        values = squareform(condensed)
    return values
