# This is a simple function designed to raise a warning whenever the sketch size is too small
# for the set being estimated. See the discussion in issue #1798
import numpy as np


def chernoff(set_size, scale, relative_error=0.05):
    """
    Computes the probability that the estimate: sketch_size * scaled deviates from the true
    set_size by more than relative_error. This relies on the fact that the sketch_size
    is binomially distributed with parameters sketch_size and 1/scale. The two-sided Chernoff
    bounds are used.
    @param set_size: The number of distinct k-mers in the given set
    @param relative_error: the desired relative error (defaults to 5%)
    @return: float (the upper bound probability)
    """
    upper_bound = 1 - 2 * np.exp(- relative_error**2*set_size/(scale * 3))
    return upper_bound


def get_set_size(scale, num_sketches):
    """
    This returns the expected number of distinct k-mers from the scale size and number of sketches.
    @param scale: the scale factor used
    @param num_sketches: the actual size of the FracMinHash sketch
    @return: int (the expected number of distinct k-mers in the original set of k-mers)
    """
    #  TODO: replace with HLL when that gets implemented
    return int(np.floor(scale * num_sketches))


def sketch_size_is_large_enough(scale, num_sketches, relative_error=0.05, confidence=0.95):
    set_size = get_set_size(scale, num_sketches)
    probability = chernoff(set_size, scale, relative_error)
    if probability >= confidence:
        return True
    else:
        return False
