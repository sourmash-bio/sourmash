"""
RevIndex - a rust-based reverse index by hashes.
"""

import weakref

from sourmash.index import Index, IndexSearchResult, _check_select_parameters
from sourmash.minhash import MinHash
from sourmash.signature import SourmashSignature
from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str, encode_str
import sourmash._lowlevel
from sourmash.minhash import flatten_and_intersect_scaled


class RevIndex(RustObject, Index):
    __dealloc_func__ = lib.revindex_free

    def __init__(
        self,
        *,
        signatures=None,
        signature_paths=None,
        template=None,
        threshold=0,
        queries=None,
        keep_sigs=False,
    ):
        self.template = template
        self.threshold = threshold
        self.queries = queries
        self.keep_sigs = keep_sigs
        self.signature_paths = signature_paths
        self._signatures = signatures

        if signature_paths is None or signatures is None:
            # delay initialization
            self._objptr = ffi.NULL
        else:
            self._init_inner()

    def _init_inner(self):
        if self._objptr != ffi.NULL:
            # Already initialized
            return

        if (
            self.signature_paths is None
            and not self._signatures
            and self._objptr == ffi.NULL
        ):
            raise ValueError("No signatures provided")
        elif (self.signature_paths or self._signatures) and self._objptr != ffi.NULL:
            raise NotImplementedError("Need to update RevIndex")

        attached_refs = weakref.WeakKeyDictionary()

        queries_ptr = ffi.NULL
        queries_size = 0
        if self.queries:
            # get list of rust objects
            collected = []
            for obj in queries:
                rv = obj._get_objptr()
                attached_refs[rv] = obj
                collected.append(rv)
            queries_ptr = ffi.new("SourmashSignature*[]", collected)
            queries_size = len(queries)

        template_ptr = ffi.NULL
        if self.template:
            if isinstance(self.template, MinHash):
                template_ptr = self.template._get_objptr()
            else:
                raise ValueError("Template must be a MinHash")

        search_sigs_ptr = ffi.NULL
        sigs_size = 0
        collected = []
        if self.signature_paths:
            for path in self.signature_paths:
                collected.append(encode_str(path))
            search_sigs_ptr = ffi.new("SourmashStr*[]", collected)
            sigs_size = len(signature_paths)

            self._objptr = rustcall(
                lib.revindex_new_with_paths,
                search_sigs_ptr,
                sigs_size,
                template_ptr,
                self.threshold,
                queries_ptr,
                queries_size,
                self.keep_sigs,
            )
        elif self._signatures:
            # force keep_sigs=True, and pass SourmashSignature directly to RevIndex.
            for sig in self._signatures:
                collected.append(sig._get_objptr())
            search_sigs_ptr = ffi.new("SourmashSignature*[]", collected)
            sigs_size = len(self._signatures)

            self._objptr = rustcall(
                lib.revindex_new_with_sigs,
                search_sigs_ptr,
                sigs_size,
                template_ptr,
                self.threshold,
                queries_ptr,
                queries_size,
            )

    def signatures(self):
        self._init_inner()

        size = ffi.new("uintptr_t *")
        sigs_ptr = self._methodcall(lib.revindex_signatures, size)
        size = size[0]

        sigs = []
        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            sigs.append(sig)

        for sig in sigs:
            yield sig

        # if self._signatures:
        #    yield from self._signatures
        # else:
        #    raise NotImplementedError("Call into Rust and retrieve sigs")

    def __len__(self):
        if self._objptr:
            return self._methodcall(lib.revindex_len)
        else:
            return len(self._signatures)

    def insert(self, node):
        if self._signatures is None:
            self._signatures = []
        self._signatures.append(node)

    def save(self, path):
        pass

    @classmethod
    def load(cls, location):
        pass

    def select(self, ksize=None, moltype=None, **kwargs):
        if self.template:
            if ksize:
                self.template.ksize = ksize
            if moltype:
                self.template.moltype = moltype
        else:
            # TODO: deal with None/default values
            self.template = MinHash(ksize=ksize, moltype=moltype)

    #    def search(self, query, *args, **kwargs):
    #        """Return set of matches with similarity above 'threshold'.
    #
    #        Results will be sorted by similarity, highest to lowest.
    #
    #        Optional arguments:
    #          * do_containment: default False. If True, use Jaccard containment.
    #          * ignore_abundance: default False. If True, and query signature
    #            and database support k-mer abundances, ignore those abundances.
    #
    #        Note, the "best only" hint is ignored by LCA_Database
    #        """
    #        if not query.minhash:
    #            return []
    #
    #        # check arguments
    #        if "threshold" not in kwargs:
    #            raise TypeError("'search' requires 'threshold'")
    #        threshold = kwargs["threshold"]
    #        do_containment = kwargs.get("do_containment", False)
    #        ignore_abundance = kwargs.get("ignore_abundance", False)
    #
    #        self._init_inner()
    #
    #        size = ffi.new("uintptr_t *")
    #        results_ptr = self._methodcall(
    #            lib.revindex_search,
    #            query._get_objptr(),
    #            threshold,
    #            do_containment,
    #            ignore_abundance,
    #            size,
    #        )
    #
    #        size = size[0]
    #        if size == 0:
    #            return []
    #
    #        results = []
    #        for i in range(size):
    #            match = SearchResult._from_objptr(results_ptr[i])
    #            if match.score >= threshold:
    #                results.append(IndexSearchResult(match.score, match.signature, match.filename))
    #
    #        return results
    #
    #    def gather(self, query, *args, **kwargs):
    #        "Return the match with the best Jaccard containment in the database."
    #        if not query.minhash:
    #            return []
    #
    #        self._init_inner()
    #
    #        threshold_bp = kwargs.get("threshold_bp", 0.0)
    #        threshold = threshold_bp / (len(query.minhash) * self.scaled)
    #
    #        results = []
    #        size = ffi.new("uintptr_t *")
    #        results_ptr = self._methodcall(
    #            lib.revindex_gather, query._get_objptr(), threshold, True, True, size
    #        )
    #        size = size[0]
    #        if size == 0:
    #            return []
    #
    #        results = []
    #        for i in range(size):
    #            match = SearchResult._from_objptr(results_ptr[i])
    #            if match.score >= threshold:
    #                results.append(IndexSearchResult(match.score, match.signature, match.filename))
    #
    #        results.sort(reverse=True,
    #                     key=lambda x: (x.score, x.signature.md5sum()))
    #
    #        return results[:1]

    @property
    def scaled(self):
        return self._methodcall(lib.revindex_scaled)


class SearchResult(RustObject):
    __dealloc_func__ = lib.searchresult_free

    def __repr__(self):
        return f"SearchResult({self.score}, {self.signature}, {self.location})"

    def __iter__(self):
        return iter((self.score, self.signature, self.location))

    @property
    def score(self):
        return self._methodcall(lib.searchresult_score)

    @property
    def signature(self):
        sig_ptr = self._methodcall(lib.searchresult_signature)
        return SourmashSignature._from_objptr(sig_ptr)

    @property
    def location(self):
        result = decode_str(self._methodcall(lib.searchresult_filename))
        if result == "":
            return None
        return result


class DiskRevIndex_DatasetPicklist(RustObject):
    __dealloc_func__ = lib.dataset_picklist_free

    def __init__(self, idxs):
        idx_list = list(idxs)
        idx_list_size = len(idx_list)

        self._objptr = rustcall(
            lib.dataset_picklist_new_from_list, idx_list, idx_list_size
        )


class DiskRevIndex(RustObject):
    __dealloc_func__ = lib.disk_revindex_free
    is_database = True
    manifest = None

    def __init__(self, path, *, ptr=None):
        path_b = path.encode("utf-8")
        if ptr is None:
            self._objptr = rustcall(lib.disk_revindex_new_from_rocksdb, path_b)
        self.location = path

    @classmethod
    def from_sigs(self, siglist, path):
        path_b = path.encode("utf-8")

        collected = []
        for ss in siglist:
            rv = ss._get_objptr()
            collected.append(rv)

        sigs_ptr = ffi.new("SourmashSignature*[]", collected)
        sig_size = len(collected)

        _ = rustcall(lib.disk_revindex_new_with_sigs, sigs_ptr, sig_size, path_b)

        return DiskRevIndex(path)

    def __len__(self):
        return self._methodcall(lib.disk_revindex_len)

    def select(
        self,
        ksize=None,
        moltype=None,
        scaled=None,
        num=None,
        abund=None,
        containment=None,
        picklist=None,
        **kwargs,
    ):
        _check_select_parameters(
            ksize=ksize,
            moltype=moltype,
            scaled=scaled,
            num=num,
            abund=abund,
            containment=containment,
            picklist=picklist,
            **kwargs,
        )

        assert not abund
        assert num is None or num == 0
        # ignore containment!

        my_ksize = self._methodcall(lib.disk_revindex_ksize)
        my_scaled = self._methodcall(lib.disk_revindex_scaled)
        my_moltype = self._methodcall(lib.disk_revindex_moltype)

        if ksize is not None:
            if ksize != my_ksize:
                raise ValueError(f"revindex ksize is {my_ksize}, not {ksize}")
        if scaled is not None and scaled < my_scaled:
            raise ValueError(f"revindex scaled is {my_scaled}, not {scaled}")
        if 0 and moltype is not None and moltype != my_moltype:
            raise ValueError(f"revindex moltype is {my_moltype}, not {moltype}")

        return self

    def signatures(self):
        size = ffi.new("uintptr_t *")
        sigs_ptr = self._methodcall(lib.disk_revindex_signatures, size)
        size = size[0]

        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            yield sig

    def signatures_with_location(self):
        for ss in self.signatures():
            yield ss, self.location

    def _signatures_with_internal(self):
        for n, ss in enumerate(self.signatures()):
            yield ss, n

    def prefetch(self, query_ss, threshold_bp=0, **kwargs):
        if not query_ss.minhash:
            raise ValueError("empty query")

        threshold_bp = int(threshold_bp)

        size = ffi.new("uintptr_t *")
        results_ptr = self._methodcall(
            lib.disk_revindex_prefetch,
            query_ss._get_objptr(),
            threshold_bp,
            size,
        )
        size = size[0]

        matches = []
        for i in range(size):
            match = SearchResult._from_objptr(results_ptr[i])
            matches.append(match)
        return matches

    def search(
        self,
        query_ss,
        *,
        threshold=None,
        do_containment=False,
        do_max_containment=False,
        best_only=False,
        picklist=None,
        **kwargs,
    ):
        # @CTB: best_only? sorting?
        if not query_ss.minhash:
            raise ValueError("empty query")

        if threshold is None:
            raise TypeError("'search' requires 'threshold'")

        size = ffi.new("uintptr_t *")
        if do_containment:
            # calculate threshold_bp from threshold
            query_mh = query_ss.minhash
            threshold_bp = int(round(threshold * len(query_mh) * query_mh.scaled))
            results_ptr = self._methodcall(
                lib.disk_revindex_prefetch,
                query_ss._get_objptr(),
                threshold_bp,
                size,
            )
        elif do_max_containment:
            raise NotImplementedError(
                "max_containment is not (yet) available on RocksDB"
            )
        else:  # jaccard
            if picklist is None:
                pl_ptr = ffi.NULL
            else:
                pl_ptr = picklist._objptr
            results_ptr = self._methodcall(
                lib.disk_revindex_search_jaccard,
                query_ss._get_objptr(),
                threshold,
                size,
                pl_ptr,
            )

        size = size[0]

        matches = []
        for i in range(size):
            match = SearchResult._from_objptr(results_ptr[i])
            matches.append(match)

        return matches

    def best_containment(self, query_ss, *, threshold_bp=0, **kwargs):
        if not query_ss.minhash:
            raise ValueError("empty query")

        threshold_bp = int(threshold_bp)

        try:
            ss_ptr = self._methodcall(
                lib.disk_revindex_best_containment, query_ss._get_objptr(), threshold_bp
            )
            match_ss = SourmashSignature._from_objptr(ss_ptr)
            if not match_ss.minhash:
                raise ValueError("no results")
        except:
            raise ValueError("no results")
        containment = query_ss.contained_by(match_ss)

        return IndexSearchResult(containment, match_ss, self.location)

    def peek(self, query_mh, *, threshold_bp=0):
        ss_ptr = self._methodcall(
            lib.disk_revindex_peek, query_mh._get_objptr(), int(threshold_bp)
        )

        match_ss = SourmashSignature._from_objptr(ss_ptr)
        if not match_ss:
            return []

        intersect_mh = flatten_and_intersect_scaled(match_ss.minhash, query_mh)
        containment = intersect_mh.contained_by(query_mh)

        return (IndexSearchResult(containment, match_ss, self.location), intersect_mh)

    def consume(self, intersect_mh):
        pass

    def counter_gather(self, query, threshold_bp, **kwargs):
        counter = DiskRevIndex_CounterGather(query, self, threshold_bp)
        for result in self.prefetch(query, threshold_bp=threshold_bp):
            counter.add(result.signature)

        return counter


class DiskRevIndex_CounterGather:
    def __init__(self, query, db, threshold_bp):
        self.query = query
        self.orig_query_mh = query.minhash.copy().flatten()
        self.found_mh = query.minhash.copy_and_clear().to_mutable()
        self.db = db
        self.threshold_bp = threshold_bp

    def add(self, match):
        query_mh = self.orig_query_mh
        match_mh = match.minhash.downsample(scaled=query_mh.scaled)
        intersect_mh = query_mh.intersection(match_mh)
        self.found_mh += intersect_mh

    def peek(self, query_mh, *, threshold_bp=None):
        if threshold_bp is None:
            threshold_bp = self.threshold_bp
            assert 0  # @CTB
        return self.db.peek(query_mh, threshold_bp=threshold_bp)

    def consume(self, intersect_mh):
        self.found_mh += intersect_mh

    @property
    def union_found(self):
        return self.found_mh

    def signatures(self):
        for sr in self.db.prefetch(self.query):
            yield sr.signature
