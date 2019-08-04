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


def similarity(sig1, sig2, ignore_abundance, downsample):
    "Compute similarity with the other MinHash signature."
    try:
        sig = sig1.minhash.similarity(sig2.minhash, ignore_abundance)
        return sig
    except ValueError as e:
        if 'mismatch in max_hash' in str(e) and downsample:
            xx = sig1.minhash.downsample_max_hash(sig2.minhash)
            yy = sig2.minhash.downsample_max_hash(sig1.minhash)
            sig = similarity(xx, yy, ignore_abundance)
            return sig
        else:
            raise


def similarity_args_unpack(args, ignore_abundance, downsample):
    return similarity(*args, ignore_abundance, downsample)

def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None):

    def nCr(n,r): 
        import math
        f = math.factorial 
        return f(n) // (f(r) * f(n-r))

    if n_jobs is None or n_jobs == 1:
        similarities = _compare_serial(siglist, ignore_abundance)
    else:
        startt = time.time()
        sig_iterator = itertools.combinations(siglist, 2)
        length_combinations = nCr(len(siglist), 2)
        notify("Created combinations list")
        del siglist
        func = partial(
            similarity_args_unpack,
            ignore_abundance=ignore_abundance,
            downsample=downsample)
        notify("similarity func initialized")
        condensed = []
        chunksize, extra = divmod(length_combinations, n_jobs) 
        if extra: 
            chunksize += 1
        chunk = 0
        while True:
            pool = multiprocessing.Pool(n_jobs)
            imapchunksize, extra = divmod(length_combinations, n_jobs) 
            if extra: 
                imapchunksize += 1
            result = list(pool.imap(func, itertools.islice(sig_iterator, chunksize), chunksize=imapchunksize))
            chunk += 1
            if result:
                condensed.extend(result)
                notify("{} chunk done", chunk)
            else:
                break

        similarities = squareform(condensed)
        notify("condensed list done")
        del sig_iterator
        # 'squareform' was made for *distance* matrices not *similarity*
        # so need to replace diagonal values with 1.
        # np.fill_digonal modifies 'similarities' in-place
        notify("squareformed")
        np.fill_diagonal(similarities, 1)
        notify("filled diagonal")
        notify("time taken to compare all pairs parallely is {} seconds ", time.time() - startt)
    return similarities
