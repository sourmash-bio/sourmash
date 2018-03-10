#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
import sys

import sourmash

sys.modules[__name__] = sys.modules['sourmash']
#sys.modules[__name__] = __import__('sourmash')
