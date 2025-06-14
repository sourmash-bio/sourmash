"""
RevIndex and DiskRevIndex - Rust-based reverse indexes by hashes.
"""

import os
import weakref

from sourmash.index import Index, IndexSearchResult, _check_select_parameters
from sourmash.minhash import MinHash, flatten_and_intersect_scaled
from sourmash.signature import SourmashSignature
from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str, encode_str
import sourmash._lowlevel
from sourmash.minhash import flatten_and_intersect_scaled
from sourmash.manifest import CollectionManifest


class RevIndex(RustObject, Index):
    """
    Base class for both MemRevIndex and DiskRevIndex.

    Provides core FFI functionality to connect to Rust code, and implements
    basic RevIndex functionality based on RevIndexOps trait.
    """

    __dealloc_func__ = lib.revindex_free
    is_database = True
    location = None

    def __init__(self):
        self._objptr = ffi.NULL
        self._idx_picklist = None

    @property
    def _ffi_idx_picklist(self):
        if self._idx_picklist is None:
            return ffi.NULL
        return self._idx_picklist._objptr

    @property
    def manifest(self):
        self._init_inner()
        mf_objptr = rustcall(lib.revindex_manifest, self._objptr)
        return CollectionManifest._from_rust(mf_objptr)

    def _generate_idx_picklist_from_manifest(self, mf):
        # grab internal indices
        idx_list = [int(row["internal_location"]) for row in mf.rows]
        self._idx_picklist = RevIndex_DatasetPicklist(idx_list)

    def signatures(self, *, use_picklist=True):
        self._init_inner()

        size = ffi.new("uintptr_t *")
        if use_picklist:
            picklist_objptr = self._ffi_idx_picklist
        else:
            picklist_objptr = ffi.NULL

        sigs_ptr = self._methodcall(lib.revindex_signatures, size, picklist_objptr)

        size = size[0]

        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            yield sig

    def signatures_with_location(self):
        for ss in self.signatures():
            yield ss, self.location

    def _signatures_with_internal(self):
        # CTB note: this should return _all_ signatures, independent of
        # picklist.
        for n, ss in enumerate(self.signatures(use_picklist=False)):
            yield ss, n

    def __len__(self):
        self._init_inner()
        return self._methodcall(lib.revindex_len, self._ffi_idx_picklist)

    @property
    def scaled(self):
        scaled = self._methodcall(lib.revindex_scaled)
        return scaled

    def save(self, path):
        raise NotImplementedException

    @classmethod
    def load(cls, location):
        raise NotImplementedException

    def search(
        self,
        query_ss,
        *,
        threshold=None,
        do_containment=False,
        do_max_containment=False,
        best_only=False,
        ignore_abundance=False,
        **kwargs,
    ):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments:
          * do_containment: default False. If True, use Jaccard containment
            instead of Jaccard similarity.
          * ignore_abundance: default False. If True, and query signature
            and database support k-mer abundances, ignore those abundances.
        """
        # CTB note: could optimize for best_only, I 'spose. But only makes a
        # difference for RevIndex when searching with containment.

        if not query_ss.minhash:
            raise ValueError("empty query")

        if threshold is None:
            raise TypeError("'search' requires 'threshold'")

        self._init_inner()

        size = ffi.new("uintptr_t *")
        if do_containment:
            # calculate threshold_bp from threshold
            query_mh = query_ss.minhash
            threshold_bp = int(round(threshold * len(query_mh) * query_mh.scaled))
            results_ptr = self._methodcall(
                lib.revindex_prefetch,
                query_ss._get_objptr(),
                threshold_bp,
                size,
                self._ffi_idx_picklist,
            )
        elif do_max_containment:
            raise NotImplementedError("max_containment is not available on RevIndex")
        else:  # jaccard similarity
            results_ptr = self._methodcall(
                lib.revindex_search_jaccard,
                query_ss._get_objptr(),
                threshold,
                size,
                self._ffi_idx_picklist,
            )

        # retrieve results => SearchResult
        size = size[0]

        matches = []
        for i in range(size):
            match = SearchResult._from_objptr(results_ptr[i])
            matches.append(match)

        return matches

    def best_containment(self, query_ss, *, threshold_bp=0, **kwargs):
        """
        Return a SearchResult tuple for the sketch with the highest
        containment of the query.
        """
        if not query_ss.minhash:
            raise ValueError("empty query")

        self._init_inner()
        threshold_bp = int(threshold_bp)
        query_mh = query_ss.minhash

        ss_ptr = self._methodcall(
            lib.revindex_best_containment,
            query_mh._get_objptr(),
            threshold_bp,
            self._ffi_idx_picklist,
        )

        match_ss = SourmashSignature._from_objptr(ss_ptr)
        if not match_ss:
            raise ValueError("no match")

        containment = query_ss.contained_by(match_ss)

        return IndexSearchResult(containment, match_ss, self.location)

    def peek(self, query_mh, *, threshold_bp=0):
        """
        Return the best containment match for the query MinHash,
        plus the intersected hashes between the query and the match.

        Used for CounterGather functionality.
        """
        self._init_inner()
        ss_ptr = self._methodcall(
            lib.revindex_best_containment,
            query_mh._get_objptr(),
            int(threshold_bp),
            self._ffi_idx_picklist,
        )

        match_ss = SourmashSignature._from_objptr(ss_ptr)
        if not match_ss:
            return []

        # calculate the intersection
        match_mh = match_ss.minhash
        common_scaled = max(match_mh.scaled, query_mh.scaled)
        query_mh = query_mh.flatten().downsample(scaled=common_scaled)
        match_mh = match_mh.flatten().downsample(scaled=common_scaled)
        intersect_mh = query_mh & match_mh
        containment = query_mh.contained_by(intersect_mh)

        return (IndexSearchResult(containment, match_ss, self.location), intersect_mh)

    def consume(self, intersect_mh):
        """
        Provide CounterGather API - does nothing, in this case.
        """
        pass

    def counter_gather(self, query_ss, threshold_bp=0, **kwargs):
        """
        Return a CounterGather object that holds interim results for a
        'gather', and can be used to get iterative results.
        """
        if not query_ss.minhash:
            raise ValueError("empty query")

        self._init_inner()
        cg_ptr = self._methodcall(
            lib.revindex_prefetch_to_countergather,
            query_ss._get_objptr(),
            self._ffi_idx_picklist,
        )

        return RevIndex_CounterGather_Colors(cg_ptr, query_ss, self)

    def prefetch(self, query_ss, threshold_bp=0, **kwargs):
        """
        Return all containment matches above threshold for the query.
        """
        if not query_ss.minhash:
            raise ValueError("empty query")

        self._init_inner()
        threshold_bp = int(threshold_bp)

        size = ffi.new("uintptr_t *")
        results_ptr = self._methodcall(
            lib.revindex_prefetch,
            query_ss._get_objptr(),
            threshold_bp,
            size,
            self._ffi_idx_picklist,
        )
        size = size[0]

        matches = []
        for i in range(size):
            match = SearchResult._from_objptr(results_ptr[i])
            matches.append(match)
        return matches

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
        """
        Implement selection protocol, including picklists.

        Since RevIndex have fixed ksize and moltype and scaled, this
        mostly just checks to see that the desired ksize and moltype match,
        and the desired scaled is greater than the RevIndex scaled.

        Also does picklist matching using RevIndex-specific Idx.
        """
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
        self._init_inner()

        assert not abund
        assert num is None or num == 0
        # ignore containment!

        my_ksize = self._methodcall(lib.revindex_ksize)
        my_scaled = self._methodcall(lib.revindex_scaled)
        my_moltype = decode_str(self._methodcall(lib.revindex_moltype))

        if ksize is not None:
            if ksize != my_ksize:
                raise ValueError(f"revindex ksize is {my_ksize}, not {ksize}")
        if scaled is not None and scaled < my_scaled:
            raise ValueError(f"revindex scaled is {my_scaled}, not {scaled}")
        if moltype is not None and moltype != my_moltype:
            raise ValueError(f"revindex moltype is {my_moltype}, not {moltype}")

        if picklist is not None:
            if self._idx_picklist is not None:
                raise Exception("cannot use picklists multiple times, sorry")

            # select matching entries from our manifest:
            m = self.manifest.select_to_manifest(picklist=picklist)

            # build the internal picklist sing the internal Idx identifiers.
            self._idx_picklist = RevIndex_DatasetPicklist.from_manifest(m)

        return self


class SearchResult(RustObject):
    """
    Hold SearchResults from Rust.
    """

    __dealloc_func__ = lib.searchresult_free

    def __repr__(self):
        return f"SearchResult({self.score}, {self.signature}, {self.location})"

    def __iter__(self):
        return iter((self.score, self.signature, self.location))

    def __getitem__(self, i):
        return list(self)[i]

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


class MemRevIndex(RevIndex):
    """
    Memory-based RevIndex.
    """

    def __init__(self, *, template=None):
        """
        Create an empty MemRevIndex, holding sketches that match the template.
        """
        super().__init__()
        assert template is not None
        assert isinstance(template, MinHash)
        if template.num != 0:
            raise ValueError("must use scaled sketches")
        self.template = template.copy_and_clear()
        self._scaled = template.scaled
        self._signatures = []  # hold sketches _prior_ to construction
        self._orig_signatures = {}

    def _check_not_init(self, *, do_raise=True):
        "Confirm that this object is not initialized, optionally raising exc."
        if self._objptr != ffi.NULL:
            if do_raise:
                raise Exception("already initialized")
            return False
        return True

    def _init_inner(self):
        """
        Initialize the MemRevIndex from all signatures in self._signatures.
        """
        if self._objptr != ffi.NULL:
            # Already initialized
            return

        if not self._signatures and self._objptr == ffi.NULL:
            raise ValueError("No signatures provided")

        if self.template.scaled < self._scaled:
            self.template = self.template.downsample(scaled=self._scaled)

        # prepare FFI call
        template_ptr = self.template._get_objptr()

        search_sigs_ptr = ffi.NULL
        sigs_size = 0
        collected = []
        for sig in self._signatures:
            collected.append(sig._get_objptr())
            search_sigs_ptr = ffi.new("SourmashSignature*[]", collected)
            sigs_size = len(self._signatures)

        self._objptr = rustcall(
            lib.revindex_mem_new_with_sigs,
            search_sigs_ptr,
            sigs_size,
            template_ptr,
        )

        # provide a mapping between original signatures, and potentially
        # downsampled signatures stored in this object, based on md5sum
        # See https://github.com/sourmash-bio/sourmash/issues/3601.
        for n, (orig_ss, stored_ss) in enumerate(
            zip(self._signatures, self.signatures())
        ):
            self._orig_signatures[stored_ss.md5sum()] = orig_ss

    def insert(self, sig):
        "Add signature to internal list, tracking max scaled along way."
        self._check_not_init()

        if sig.minhash.scaled > self._scaled:
            self._scaled = sig.minhash.scaled
        self._signatures.append(sig)

    def search(self, *args, **kwargs):
        """
        Implement a search that returns the original signatures, not
        just the indexed ones.

        See also https://github.com/sourmash-bio/sourmash/issues/3601.
        """
        results = super().search(*args, **kwargs)
        results2 = []
        for match in results:
            match_md5 = match.signature.md5sum()
            orig_ss = self._orig_signatures[match_md5]
            results2.append(IndexSearchResult(match.score, orig_ss, match.location))
        return results2


class DiskRevIndex(RevIndex):
    """
    RocksDB-based low-memory on disk inverted index, implemented in Rust.
    """

    __dealloc_func__ = lib.revindex_free
    is_database = True

    def __init__(self, path):
        """
        Initialize based on a pre-existing RocksDB.
        """
        super().__init__()
        check_file = os.path.join(path, "CURRENT")
        if not os.path.exists(check_file):
            raise ValueError("not a RocksDB")

        # create via FFI
        path_b = path.encode("utf-8")
        self._objptr = rustcall(lib.revindex_new_from_rocksdb, path_b)

        # store location
        self._path = path

    def _init_inner(self):
        pass

    @property
    def location(self):
        return self._path

    def insert(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def create_from_sigs(self, siglist, path):
        """
        Create a _new_ DiskRevIndex based on a list of signatures.
        """
        path_b = path.encode("utf-8")

        collected = []
        for ss in siglist:
            rv = ss._get_objptr()
            collected.append(rv)

        sigs_ptr = ffi.new("SourmashSignature*[]", collected)
        sig_size = len(collected)

        _ = rustcall(lib.revindex_disk_create, sigs_ptr, sig_size, path_b)

        return DiskRevIndex(path)


class RevIndex_CounterGather:
    """
    Simple implementation of CounterGather API that tracks matches
    while passing most calls back to the parent RevIndex.

    CTB note: This is not used in the code base currently, but _is_
    tested in test_index_protocol, so I'm leaving it in for now.
    """

    def __init__(self, query, db, threshold_bp, *, allow_insert=False):
        """
        Initialize a CounterGather obj.

        Here, 'db' can be either a RevIndex or a DiskRevIndex.
        """
        self.query = query
        self.orig_query_mh = query.minhash.copy().flatten()
        self.found_mh = query.minhash.copy_and_clear().to_mutable()
        self.db = db
        self.threshold_bp = threshold_bp
        self.allow_insert = allow_insert
        self.locations = {}

    def __len__(self):
        return len(self.db)

    @property
    def scaled(self):
        return self.db.scaled

    def add(self, match_ss, *, location=None, require_overlap=True):
        if self.allow_insert:
            if self.db._check_not_init(do_raise=False):
                self.db.insert(match_ss)
            else:
                raise ValueError

        self.locations[match_ss.md5sum()] = location

        query_mh = self.orig_query_mh
        match_mh = match_ss.minhash
        intersect_mh = flatten_and_intersect_scaled(query_mh, match_mh)
        if require_overlap and not intersect_mh:
            raise ValueError("require overlap")

        if self.found_mh.scaled < intersect_mh.scaled:
            self.found_mh = self.found_mh.downsample(scaled=intersect_mh.scaled)

        self.found_mh += intersect_mh

    def peek(self, query_mh, *, threshold_bp=0):
        """
        CounterGather API: return best matching minhash, plus intersection.
        """
        if not query_mh:
            return []

        if query_mh.contained_by(self.orig_query_mh, True) != 1.0:
            raise ValueError
        # assert threshold_bp is not None

        res = self.db.peek(query_mh, threshold_bp=threshold_bp)
        if not res:
            return []

        sr, intersect_mh = res
        sr_ss = sr.signature
        sr_score = sr.score
        new_sr = IndexSearchResult(sr_score, sr_ss, self.locations.get(sr_ss.md5sum()))
        return new_sr, intersect_mh

    def consume(self, intersect_mh):
        "CounterGather API: track found sketches."
        self.db._init_inner()

        if self.found_mh.scaled < intersect_mh.scaled:
            self.found_mh = self.found_mh.downsample(scaled=intersect_mh.scaled)
        elif self.found_mh.scaled > intersect_mh.scaled:
            intersect_mh = intersect_mh.downsample(scaled=self.found_mh.scaled)

        # CTB: is this right?
        self.found_mh += intersect_mh

    @property
    def union_found(self):
        "Return all hashes found."
        return self.found_mh

    def signatures(self):
        # don't track actual signatures - go back to RevIndex
        # this is probably overkill.
        for sr in self.db.prefetch(self.query):
            yield sr.signature


class RevIndex_CounterGather_Colors(RustObject):
    """
    Implementation of CounterGather using RevIndex color-based counters.
    """

    __dealloc_func__ = lib.revindex_countergather_free

    def __init__(self, objptr, query_ss, db):
        """
        Track originating query & RevIndex.
        """
        self._objptr = objptr
        assert isinstance(db, RevIndex)
        self.db = db
        query_mh = query_ss.minhash
        self._scaled = query_mh.scaled

        if len(self) == 0:
            raise ValueError("no matches found")

        empty_mh = query_mh.copy_and_clear()

        # track found hashes:
        found_mh_ptr = self._methodcall(
            lib.revindex_countergather_found_hashes,
            empty_mh._objptr,
        )
        self.found_mh = MinHash._from_objptr(found_mh_ptr)

    @property
    def scaled(self):
        return self._scaled

    def add(self, match_ss, *, location=None, require_overlap=True):
        raise NotImplementedError

    def peek(self, query_mh, *, threshold_bp=0):
        """
        CounterGather API: return matching SearchResult + intersecting minhash.
        """
        threshold_hashes = int(threshold_bp / query_mh.scaled)
        match_ss_ptr = self._methodcall(
            lib.revindex_countergather_peek,
            self.db._objptr,
            threshold_hashes,
        )

        # empty SourmashSignature => nothing found.
        match_ss = SourmashSignature._from_objptr(match_ss_ptr)
        if not match_ss:
            return []

        match_mh = match_ss.minhash
        intersect_mh = flatten_and_intersect_scaled(query_mh, match_mh)
        containment = len(intersect_mh) / len(query_mh)

        return (
            IndexSearchResult(containment, match_ss, self.db.location),
            intersect_mh,
        )

    def consume(self, intersect_mh):
        "CounterGather API: decrement counters."
        _ = self._methodcall(lib.revindex_countergather_consume, intersect_mh._objptr)

    @property
    def union_found(self):
        "Return all found hashes."
        return self.found_mh

    def __len__(self):
        "Return number of distinct matching signatures remaining."
        return self._methodcall(lib.revindex_countergather_len)

    def signatures(self):
        "Return all signatures found."
        size = ffi.new("uintptr_t *")
        sigs_ptr = self._methodcall(
            lib.revindex_countergather_signatures, self.db._objptr, size
        )
        size = size[0]

        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            yield sig


class RevIndex_DatasetPicklist(RustObject):
    """
    Intermediate class for holding lists of Idx, internal identifiers
    used by RevIndex structs in Rust.
    """

    __dealloc_func__ = lib.dataset_picklist_free

    def __init__(self, idxs):
        idx_list = list(idxs)
        idx_list_size = len(idx_list)

        # CTB do some validation on Idx?
        self._objptr = rustcall(
            lib.dataset_picklist_new_from_list, idx_list, idx_list_size
        )

    @classmethod
    def from_manifest(cls, mf):
        idx_list = [int(row["internal_location"]) for row in mf.rows]
        return cls(idx_list)
