# This is a simple function designed to raise a warning whenever the sketch size is too small
# for the set being estimated. See the discussion in issue #1798
import numpy as np


def set_size_chernoff(set_size, scale, relative_error=0.05):
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


def set_size_estimate_is_accurate(scale, num_sketches, relative_error=0.05, confidence=0.95):
    set_size = get_set_size(scale, num_sketches)
    probability = set_size_chernoff(set_size, scale, relative_error)
    if probability >= confidence:
        return True
    else:
        return False


def test_set_size_chernoff():
    eps = 10**(-6)
    rel_error = 0.01
    set_size = 1000000
    s = 1/0.1  # I'm used to using a scale value between 0 and 1
    value_from_mathematica = 0.928652
    assert np.abs(set_size_chernoff(set_size, s, rel_error) - value_from_mathematica) < eps

    rel_error = 0.05
    set_size = 10000
    s = 1
    value_from_mathematica = 0.999519
    assert np.abs(set_size_chernoff(set_size, s, rel_error) - value_from_mathematica) < eps

    rel_error = 0.001
    set_size = 10
    s = 1/.01
    value_from_mathematica = -1
    assert np.abs(set_size_chernoff(set_size, s, rel_error) - value_from_mathematica) < eps


def test_set_size_estimate_is_accurate():
    eps = 10 ** (-6)
    rel_error = 0.05
    set_size = 1000000
    s = 1 / 0.1  # I'm used to using a scale value between 0 and 1
    num_sketches = set_size / s  # idealized case
    confidence = 0.95
    assert set_size_estimate_is_accurate(scale=s, num_sketches=num_sketches, relative_error=rel_error, confidence=confidence) is True
    confidence = set_size_chernoff(set_size=set_size, scale=s, relative_error=rel_error)
    assert set_size_estimate_is_accurate(scale=s, num_sketches=num_sketches, relative_error=rel_error, confidence=confidence) is True
    # Horrible values
    assert set_size_estimate_is_accurate(scale=10000, num_sketches=num_sketches, relative_error=0, confidence=1) is False
    # Less horrible, but still bad values
    confidence = set_size_chernoff(set_size=set_size, scale=s, relative_error=rel_error)
    assert set_size_estimate_is_accurate(scale=s, num_sketches=num_sketches, relative_error=rel_error, confidence=confidence*2) is False
    # one where the confidence is negative
    rel_error = .001
    set_size = 10
    s = 100
    num_sketches = set_size/s
    assert set_size_estimate_is_accurate(scale=s, num_sketches=num_sketches, relative_error=rel_error, confidence=confidence) is False
    assert set_size_estimate_is_accurate(scale=s, num_sketches=0, relative_error=rel_error, confidence=confidence) is False


def run_tests():
    test_set_size_chernoff()
    test_set_size_estimate_is_accurate()


if __name__ == '__main__':
    print("Running tests")
    run_tests()
    print("Tests completed successfully")
