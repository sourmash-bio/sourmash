import pytest

from sourmash.utils import RustObject
from sourmash.minhash import to_bytes


def test_rustobj_init():
    with pytest.raises(TypeError):
        RustObject()


def test_to_bytes():
    with pytest.raises(TypeError):
        to_bytes([9882])

    assert to_bytes(98) == bytes([98])
    assert to_bytes("abc") == b"abc"
    assert to_bytes(b"abc") == b"abc"
