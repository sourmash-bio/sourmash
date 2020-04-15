# -*- coding: UTF-8 -*-
from __future__ import unicode_literals, division

from tempfile import NamedTemporaryFile

from ._compat import string_types, range_type
from ._lowlevel import ffi, lib
from ._minhash import to_bytes, MinHash
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError


class Nodegraph(RustObject):
    __dealloc_func__ = lib.nodegraph_free

    def __init__(self, ksize, starting_size, n_tables):
        self._objptr = lib.nodegraph_with_tables(ksize, int(starting_size), n_tables)

    @staticmethod
    def load(filename):
        ng_ptr = rustcall(lib.nodegraph_from_path, to_bytes(filename))
        return Nodegraph._from_objptr(ng_ptr)

    @staticmethod
    def from_buffer(buf):
        ng_ptr = rustcall(lib.nodegraph_from_buffer, buf, len(buf))
        return Nodegraph._from_objptr(ng_ptr)

    def save(self, filename):
        self._methodcall(lib.nodegraph_save, to_bytes(filename))

    def update(self, other):
        if isinstance(other, Nodegraph):
            return self._methodcall(lib.nodegraph_update, other._objptr)
        elif isinstance(other, MinHash):
            return self._methodcall(lib.nodegraph_update_mh, other._objptr)
        else:
            raise TypeError("Must be a Nodegraph or MinHash")

    def count(self, h):
        if isinstance(h, string_types):
            return self._methodcall(lib.nodegraph_count_kmer, to_bytes(h))
        return self._methodcall(lib.nodegraph_count, h)

    def get(self, h):
        if isinstance(h, string_types):
            return self._methodcall(lib.nodegraph_get_kmer, to_bytes(h))
        return self._methodcall(lib.nodegraph_get, h)

    @property
    def n_occupied(self):
        return self._methodcall(lib.nodegraph_noccupied)

    @property
    def ksize(self):
        return self._methodcall(lib.nodegraph_ksize)

    @property
    def tablesize(self):
        return self._methodcall(lib.nodegraph_tablesize)

    @property
    def n_tables(self):
        return self._methodcall(lib.nodegraph_ntables)

    @property
    def expected_collisions(self):
        return self._methodcall(lib.nodegraph_expected_collisions)

    def matches(self, mh):
        if not isinstance(mh, MinHash):
            raise ValueError("mh must be a MinHash")

        return self._methodcall(lib.nodegraph_matches, mh._objptr)

    def to_khmer_nodegraph(self):
        import khmer
        try:
            load_nodegraph = khmer.load_nodegraph
        except AttributeError:
            load_nodegraph = khmer.Nodegraph.load

        with NamedTemporaryFile() as f:
            self.save(f.name)
            f.file.flush()
            f.file.seek(0)
            return load_nodegraph(f.name)


def extract_nodegraph_info(filename):
    ng = Nodegraph.load(filename)
    return ng.ksize, ng.tablesize, ng.n_tables


def calc_expected_collisions(graph, force=False, max_false_pos=.2):
    fp_all = graph.expected_collisions()

    if fp_all > max_false_pos:
        print("**", file=sys.stderr)
        print("** ERROR: the graph structure is too small for ",
              file=sys.stderr)
        print("** this data set.  Increase data structure size.",
              file=sys.stderr)
        print("** Do not use these results!!", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** (estimated false positive rate of %.3f;" % fp_all,
              file=sys.stderr, end=' ')
        print("max recommended %.3f)" % max_false_pos, file=sys.stderr)
        print("**", file=sys.stderr)

        if not force:
            sys.exit(1)

    return fp_all
