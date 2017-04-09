#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re
import math
from ._minhash import (MinHash, dotproduct, get_minhash_default_seed,
                       get_minhash_max_hash)

DEFAULT_SEED = get_minhash_default_seed()
MAX_HASH = get_minhash_max_hash()
