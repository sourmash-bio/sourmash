import itertools
from functools import partial
import os
import tempfile
import time
import multiprocessing
from scipy.spatial.distance import squareform
import numpy as np

from .logging import notify, error, print_results


def _compare_serial(siglist, ignore_abundance):
    n = len(siglist)

    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = itertools.combinations(range(n), 2)

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


class MultiprocessCompare():
    sig_iterator = None

    @classmethod
    def set_sig_iterator(cls, sig_iterator):
        cls.sig_iterator = sig_iterator

    @classmethod
    def similarity(cls, ignore_abundance, downsample, index):
        "Compute similarity with the other MinHash signature."
        startt = time.time()
        sig1, sig2 = cls.sig_iterator[index]
        try:
            sig = sig1.minhash.similarity(sig2.minhash, ignore_abundance)
            notify("time taken in calculating similarity is {} seconds ", time.time() - startt)
            return sig
        except ValueError as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = sig1.minhash.downsample_max_hash(sig2.minhash)
                yy = sig2.minhash.downsample_max_hash(sig1.minhash)
                sig = similarity(xx, yy, ignore_abundance)
                notify("time taken in calculating similarity is {} seconds ", time.time() - startt)
                return sig
            else:
                raise


def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None):

    def nCr(n,r): 
        import math
        f = math.factorial 
        return f(n) // (f(r) * f(n-r))

    if n_jobs is None or n_jobs == 1:
        similarities = _compare_serial(siglist, ignore_abundance)
    else:
        length_combinations = nCr(len(siglist), 2)
        sig_iterator = list(itertools.combinations(siglist, 2))
        notify("Created combinations list")
        del siglist
        multiprocess_compare = MultiprocessCompare()
        multiprocess_compare.set_sig_iterator(sig_iterator)
        notify("Created MultiprocessCompare class and set sig iterator")
        del sig_iterator
        func = partial(multiprocess_compare.similarity, ignore_abundance, downsample)
        notify("similarity func initialized")
        with multiprocessing.Pool(n_jobs) as pool:
            condensed = list(pool.imap(func, range(length_combinations)))
        notify("condensed list done")
        similarities = squareform(condensed)
        # 'squareform' was made for *distance* matrices not *similarity*
        # so need to replace diagonal values with 1.
        # np.fill_digonal modifies 'similarities' in-place
        notify("squareformed")
        np.fill_diagonal(similarities, 1)
        notify("filled diagonal")
    return similarities
