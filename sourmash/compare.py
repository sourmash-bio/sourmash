import itertools
import math
from functools import partial
import os
import tempfile
import time
import multiprocessing
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


def similarity_args(args, ignore_abundance, downsample):
    return similarity(*args, ignore_abundance, downsample)


def similarity_args_unpack(index, ignore_abundance, downsample, siglist):
    startt = time.time()
    sig_iterator = itertools.product([siglist[index]], siglist[index + 1:])
    func = partial(similarity_args,
                   ignore_abundance=ignore_abundance,
                   downsample=downsample)
    ret = list(map(func, sig_iterator))
    notify("comparison for index {} done in {:.5f} seconds", index, time.time() - startt)
    return ret


def nCr(n, r):
    f = math.factorial
    return f(n) // (f(r) * f(n - r))


def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None):

    if n_jobs is None or n_jobs == 1:
        similarities = _compare_serial(siglist, ignore_abundance)
    else:
        startt = time.time()
        length_siglist = len(siglist)
        siglist = memmap_siglist(siglist)
        notify("Created memmapped siglist")
        func = partial(
            similarity_args_unpack,
            siglist=siglist,
            ignore_abundance=ignore_abundance,
            downsample=downsample)
        notify("Created similarity func")

        length_combinations = nCr(length_siglist, 2)
        d = int(np.ceil(np.sqrt(length_combinations * 2)))

        # Check that v is of valid dimensions.
        if d * (d - 1) != length_combinations * 2:
            raise ValueError('Incompatible vector size. It must be a binomial '
                             'coefficient n choose 2 for some integer n >= 2.')
        similarities = np.eye(d, dtype=np.float64)
        temp_folder = tempfile.mkdtemp()
        filename = os.path.join(temp_folder, 'similarities.mmap')
        if os.path.exists(filename):
            os.unlink(filename)
        shape = (d, d)
        memmap_similarities = np.memmap(filename, mode='w+', shape=shape, dtype=np.float64)
        for i in range(d):
            memmap_similarities[i, i] = similarities[i, i]
        notify("Initialized memmapped similarities matrix")

        with multiprocessing.Pool(n_jobs) as pool:
            chunksize, extra = divmod(length_siglist, n_jobs)
            if extra:
                chunksize += 1
            notify("Calculated chunk size")
            result = pool.imap(func, range(length_siglist), chunksize=chunksize)
            notify("Pool.imap completed")
            for index, l in enumerate(result):
                startt = time.time()
                col_idx = index + 1
                for idx_condensed, item in enumerate(l):
                    memmap_similarities[index, col_idx + idx_condensed] = memmap_similarities[idx_condensed + col_idx, index] = item
                notify("setting similarities matrix for index {} done in {:.5f} seconds", index, time.time() - startt)
            notify("Setting similarities completed")

        del siglist
        del memmap_similarities
    large_memmap = np.memmap(filename, dtype=np.float64, shape=shape)
    notify("time taken to compare all pairs parallely is {:.5f} seconds ", time.time() - startt)
    return large_memmap
