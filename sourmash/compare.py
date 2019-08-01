from itertools import combinations
from functools import partial
import os
import tempfile
import multiprocessing
from multiprocessing.pool import ThreadPool
from scipy.spatial.distance import squareform
import numpy as np

from .logging import notify, error, print_results


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
    if os.path.exists(filename):
        os.unlink(filename)
    nrows = len(siglist)
    f = np.memmap(filename, mode='w+', shape=(nrows), dtype=object)
    for i in range(nrows):
        f[i] = siglist[i]
    del f
    large_memmap = np.memmap(filename, dtype=object, shape=(nrows))
    return large_memmap


def similarity(ignore_abundance, downsample, sig1, sig2):
    "Compute similarity with the other MinHash signature."

    print("calculating similarity")
    try:
        return sig1.minhash.similarity(sig2.minhash, ignore_abundance)
    except ValueError as e:
        if 'mismatch in max_hash' in str(e) and downsample:
            xx = sig1.minhash.downsample_max_hash(sig2.minhash)
            yy = sig2.minhash.downsample_max_hash(sig1.minhash)
            return similarity(xx, yy, ignore_abundance)
        else:
            raise


def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None):
    print("in compare_all_pairs jobs {}", n_jobs)
    if n_jobs is None or n_jobs == 1:
        similarities = _compare_serial(siglist, ignore_abundance)
    else:
        # Create a memory-mapped array
        print("memory mapping initialized")
        memmapped = memmap_siglist(siglist)
        del siglist
        print("memory mapped complete")
        sig_iterator = combinations(memmapped, 2)
        print("combinations sig iterator created")

        # do some other stuff in the main process
        pool = ThreadPool(processes=n_jobs)
        print("ThreadPool initialized")
        func = partial(similarity, ignore_abundance, downsample)
        print("partial func initialized")

        condensed = pool.starmap(
            func,
            [(arg1, arg2) for arg1, arg2 in sig_iterator])
            
        pool.close()
        pool.join()
        similarities = squareform(condensed)

        # 'squareform' was made for *distance* matrices not *similarity*
        # so need to replace diagonal values with 1.
        # np.fill_digonal modifies 'similarities' in-place
        np.fill_diagonal(similarities, 1)
    return similarities
