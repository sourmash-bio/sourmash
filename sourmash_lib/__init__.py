#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.

Legacy / deprecated; will be removed in sourmash 4.0.
"""
import sys
import warnings

warnings.warn("Please import sourmash, instead of sourmash_lib; sourmash_lib will be removed in 4.x", FutureWarning)

import sourmash

sys.modules[__name__] = sys.modules['sourmash']
#sys.modules[__name__] = __import__('sourmash')
