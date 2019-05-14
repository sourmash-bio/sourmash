#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function

import math
import os
import re

from pkg_resources import DistributionNotFound, get_distribution

from . import lca, sbt, sbt_storage, sbtmh, signature
from ._minhash import MinHash, get_minhash_default_seed, get_minhash_max_hash
from .sbtmh import create_sbt_index, load_sbt_index, search_sbt_index
from .signature import (SourmashSignature, load_one_signature, load_signatures,
                        save_signatures)

DEFAULT_SEED = get_minhash_default_seed()
MAX_HASH = get_minhash_max_hash()


try:
    VERSION = get_distribution(__name__).version
except DistributionNotFound:  # pragma: no cover
    try:
        from .version import version as VERSION  # noqa
    except ImportError:  # pragma: no cover
        raise ImportError(
            "Failed to find (autogenerated) version.py. "
            "This might be because you are installing from GitHub's tarballs, "
            "use the PyPI ones."
        )
