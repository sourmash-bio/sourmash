import numpy as np
from sourmash.signature import SourmashSignature
import sourmash
from sourmash.np_utils import to_memmap


def test_memmap():

    e1 = sourmash.MinHash(n=1, ksize=20)
    sig1 = SourmashSignature(e1)

    e2 = sourmash.MinHash(n=1, ksize=25)
    sig2 = SourmashSignature(e2)
    siglist = [sig1, sig2]
    memmapped, filename = to_memmap(np.array(siglist))
    # Assert that the data didn't change as a result of memory-mapping
    np.testing.assert_array_equal(memmapped, siglist)
    assert filename.endswith(".mmap")
