#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re
import math
from ._minhash import (MinHash, get_minhash_default_seed, get_minhash_max_hash)
import os

# retrieve VERSION from sourmash_lib/VERSION.
thisdir = os.path.dirname(__file__)
version_file = open(os.path.join(thisdir, 'VERSION'))
VERSION = version_file.read().strip()

DEFAULT_SEED = get_minhash_default_seed()
MAX_HASH = get_minhash_max_hash()

def scaled_to_max_hash(scaled):
    if scaled > 1:
        max_hash = MAX_HASH / float(scaled)
        max_hash = int(round(max_hash, 0))
        return max_hash
    else:
        return 0
