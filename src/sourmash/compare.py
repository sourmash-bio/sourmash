"""Functionality for comparing many signatures, used in sourmash compare."""

import itertools
from functools import partial
import time
import multiprocessing

from sourmash.sketchcomparison import FracMinHashComparison

from .logging import notify
from sourmash.np_utils import to_memmap


def compare_serial(siglist, ignore_abundance, *, downsample=False, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations serially on a single
    process. Best to use when there is few signatures.

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)
    jaccard_ani_untrustworthy = False
    potential_false_negatives = False

    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = itertools.combinations(range(n), 2)

    similarities = np.ones((n, n))

    for i, j in iterator:
        if return_ani:
            ani_result = siglist[i].jaccard_ani(siglist[j],downsample=downsample)
            if not potential_false_negatives and ani_result.p_exceeds_threshold:
                potential_false_negatives = True
            if not jaccard_ani_untrustworthy and ani_result.je_exceeds_threshold:
                jaccard_ani_untrustworthy = True
            ani = ani_result.ani
            if ani == None:
                ani = 0.0
            similarities[i][j] = similarities[j][i] = ani
        else:
            similarities[i][j] = similarities[j][i] = siglist[i].similarity(siglist[j], ignore_abundance=ignore_abundance, downsample=downsample)

    if jaccard_ani_untrustworthy:
        notify("WARNING: Jaccard estimation for at least one of these comparisons is likely inaccurate. Could not estimate ANI for these comparisons.")
    if potential_false_negatives:
        notify("WARNING: Some of these sketches may have no hashes in common based on chance alone (false negatives). Consider decreasing your scaled value to prevent this.")
    return similarities


def compare_serial_containment(siglist, *, downsample=False, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of containments. Processes combinations serially on a single
    process. Best to only use when there are few signatures.

    :param list siglist: list of signatures to compare
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)
    potential_false_negatives = False

    containments = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                containments[i][j] = 1
            elif return_ani:
                ani_result = siglist[j].containment_ani(siglist[i], downsample=downsample)
                ani = ani_result.ani
                if not potential_false_negatives and ani_result.p_exceeds_threshold:
                    potential_false_negatives = True
                if ani == None:
                    ani = 0.0
                containments[i][j] = ani
            else:
                containments[i][j] = siglist[j].contained_by(siglist[i],
                                                         downsample=downsample)

    if potential_false_negatives:
        notify("WARNING: Some of these sketches may have no hashes in common based on chance alone (false negatives). Consider decreasing your scaled value to prevent this.")

    return containments


def compare_serial_max_containment(siglist, *, downsample=False, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of max_containments. Processes combinations serially on a single
    process. Best to only use when there are few signatures.

    :param list siglist: list of signatures to compare
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)
    potential_false_negatives = False
    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = itertools.combinations(range(n), 2)

    containments = np.ones((n, n))

    for i, j in iterator:
        if return_ani:
            ani_result = siglist[j].max_containment_ani(siglist[i], downsample=downsample)
            ani = ani_result.ani
            if not potential_false_negatives and ani_result.p_exceeds_threshold:
                potential_false_negatives = True
            if ani == None:
                ani = 0.0
            containments[i][j] = containments[j][i] = ani
        else:
            containments[i][j] = containments[j][i] = siglist[j].max_containment(siglist[i],
                                                        downsample=downsample)
    if potential_false_negatives:
        notify("WARNING: Some of these sketches may have no hashes in common based on chance alone (false negatives). Consider decreasing your scaled value to prevent this.")

    return containments


def compare_serial_avg_containment(siglist, *, downsample=False, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of avg_containments. Processes combinations serially on a single
    process. Best to only use when there are few signatures.

    :param list siglist: list of signatures to compare
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)
    potential_false_negatives = False
    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = itertools.combinations(range(n), 2)

    containments = np.ones((n, n))

    for i, j in iterator:
        if return_ani:
            cmp = FracMinHashComparison(siglist[j].minhash, siglist[i].minhash)
            ani = cmp.avg_containment_ani
            if ani == None:
                ani = 0.0
            if not potential_false_negatives and cmp.potential_false_negative:
                potential_false_negatives = True
            containments[i][j] = containments[j][i] = ani
        else:
            containments[i][j] = containments[j][i] = siglist[j].avg_containment(siglist[i],
                                                        downsample=downsample)

    if potential_false_negatives:
        notify("WARNING: Some of these sketches may have no hashes in common based on chance alone (false negatives). Consider decreasing your scaled value to prevent this.")

    return containments


def similarity_args_unpack(args, ignore_abundance, *, downsample, return_ani=False):
    """Helper function to unpack the arguments. Written to use in pool.imap
    as it can only be given one argument."""
    sig1, sig2 = args
    if return_ani:
        ani = sig1.jaccard_ani(sig2, downsample=downsample).ani
        if ani == None:
            ani = 0.0
        return ani
    else:
        return sig1.similarity(sig2,
                           ignore_abundance=ignore_abundance,
                           downsample=downsample)


def get_similarities_at_index(index, ignore_abundance, downsample, siglist, *, return_ani=False):
    """Returns similarities of all the combinations of signature at index in
    the siglist with the rest of the indices starting at index + 1. Doesn't
    redundantly calculate signatures with all the other indices prior to
    index - 1

    :param int index: generate masks from this image
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :param siglist list of signatures
    :return: list of similarities for the combinations of signature at index
        with rest of the signatures from index+1
    """
    startt = time.time()
    sig_iterator = itertools.product([siglist[index]], siglist[index + 1:])
    func = partial(similarity_args_unpack,
                   ignore_abundance=ignore_abundance,
                   downsample=downsample,
                   return_ani=return_ani)
    similarity_list = list(map(func, sig_iterator))
    notify(
        f"comparison for index {index} done in {time.time() - startt:.5f} seconds", end='\r')
    return similarity_list


def compare_parallel(siglist, ignore_abundance, downsample, n_jobs, *, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations parallely on number of processes
    given by n_jobs

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :param int n_jobs number of processes to run the similarity calculations on
    :return: np.array similarity matrix
    """
    import numpy as np

    # Starting time - calculate time to keep track in case of lengthy siglist
    start_initial = time.time()

    # Create a memory map of the siglist using numpy to avoid memory burden
    # while accessing small parts in it
    siglist, _ = to_memmap(np.array(siglist))
    notify("Created memmapped siglist")

    # Check that length of combinations can result in a square similarity matrix
    length_siglist = len(siglist)

    # Initialize with ones in the diagonal as the similarity of a signature with
    # itself is one
    similarities = np.eye(length_siglist, dtype=np.float64)
    memmap_similarities, filename = to_memmap(similarities)
    notify("Initialized memmapped similarities matrix")

    # Initialize the function using func.partial with the common arguments like
    # siglist, ignore_abundance, downsample, for computing all the signatures
    # The only changing parameter that will be mapped from the pool is the index
    func = partial(
        get_similarities_at_index,
        siglist=siglist,
        ignore_abundance=ignore_abundance,
        downsample=downsample,
        return_ani=return_ani)
    notify("Created similarity func")

    # Initialize multiprocess.pool
    pool = multiprocessing.Pool(processes=n_jobs)

    # Calculate chunk size, by default pool.imap chunk size is 1
    chunksize, extra = divmod(length_siglist, n_jobs)
    if extra:
        chunksize += 1
    notify("Calculated chunk size for multiprocessing")

    # This will not generate the results yet, since pool.imap returns a generator
    result = pool.imap(func, range(length_siglist), chunksize=chunksize)
    notify("Initialized multiprocessing pool.imap")

    # Enumerate and calculate similarities at each of the indices
    # and set the results at the appropriate combination coordinate
    # locations inside the similarity matrix
    for index, l in enumerate(result):
        startt = time.time()
        col_idx = index + 1
        for idx_condensed, item in enumerate(l):
            memmap_similarities[index, col_idx + idx_condensed] = memmap_similarities[idx_condensed + col_idx, index] = item
        notify(
            f"Setting similarities matrix for index {index} done in {time.time() - startt:.5f} seconds", end='\r')
    notify("Setting similarities completed")

    pool.close()
    pool.join()

    notify(f"Time taken to compare all pairs parallely is {time.time() - start_initial:.5f} seconds ")
    return np.memmap(filename, dtype=np.float64, shape=(length_siglist, length_siglist))


def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None, return_ani=False):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations either serially or
    based on parallely on number of processes given by n_jobs

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :param int n_jobs number of processes to run the similarity calculations on,
    if number of jobs is None or 1, compare serially, otherwise parallely.
    :return: np.array similarity matrix
    """
    if n_jobs is None or n_jobs == 1:
        similarities = compare_serial(siglist, ignore_abundance=ignore_abundance, downsample=downsample, return_ani=return_ani)
    else:
        similarities = compare_parallel(siglist, ignore_abundance, downsample, n_jobs, return_ani=return_ani)
    return similarities
