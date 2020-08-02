"Legacy / deprecated; will be removed in sourmash 4.0."
import warnings

warnings.warn("Please import from the top level sourmash module instead of using _minhash, which will be renamed in 4.x", FutureWarning)

from .minhash import *
