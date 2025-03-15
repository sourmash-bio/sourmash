"""
RevIndex and DiskRevIndex - a Rust-based reverse indexes by hashes.
"""

import os
import weakref

from sourmash.index import Index, IndexSearchResult, _check_select_parameters
from sourmash.minhash import MinHash
from sourmash.signature import SourmashSignature
from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str, encode_str
import sourmash._lowlevel
from sourmash.minhash import flatten_and_intersect_scaled
from sourmash.manifest import CollectionManifest


class RevIndex(RustObject):  # , Index):
    __dealloc_func__ = lib.revindex_free
    manifest = None
    is_database = True
    location = None

    def __init__(self, *, template=None):
        assert template is not None
        assert isinstance(template, MinHash)
        if template.num != 0:
            raise ValueError("must use scaled sketches")
        self.template = template.copy_and_clear().to_mutable()
        self._scaled = template.scaled
        self._signatures = []
        self._objptr = ffi.NULL

    def _check_not_init(self, *, do_raise=True):
        if self._objptr != ffi.NULL:
            if do_raise:
                raise Exception("already initialized")
            return False
        return True

    def _init_inner(self):
        if self._objptr != ffi.NULL:
            # Already initialized
            return

        if not self._signatures and self._objptr == ffi.NULL:
            raise ValueError("No signatures provided")

        if self.template.scaled != self._scaled:
            print(f'XXX downsampling: {self.template.scaled}, {self._scaled}')
            self.template = self.template.downsample(scaled=self._scaled)

        template_ptr = self.template._get_objptr()

        search_sigs_ptr = ffi.NULL
        sigs_size = 0
        collected = []
        for sig in self._signatures:
            collected.append(sig._get_objptr())
            search_sigs_ptr = ffi.new("SourmashSignature*[]", collected)
            sigs_size = len(self._signatures)

        self._objptr = rustcall(
            lib.revindex_new_with_sigs,
            search_sigs_ptr,
            sigs_size,
            template_ptr,
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

    def signatures_with_location(self):
        for ss in self.signatures():
            yield ss, self.location  # @CTB

    def __len__(self):
        self._init_inner()
        return self._methodcall(lib.revindex_len)

    def insert(self, sig):
        if sig.minhash.scaled > self._scaled:
            self._scaled = sig.minhash.scaled

        self._check_not_init()
        self._signatures.append(sig)

    def save(self, path):
        pass

    @classmethod
    def load(cls, location):
        pass

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

        my_ksize = self.template.ksize
        my_scaled = self.template.scaled
        my_moltype = self.template.moltype

        if ksize is not None:
            if ksize != my_ksize:
                raise ValueError(f"revindex ksize is {my_ksize}, not {ksize}")
        if scaled is not None and scaled < my_scaled:
            raise ValueError(f"revindex scaled is {my_scaled}, not {scaled}")
        if moltype is not None and moltype != my_moltype:
            raise ValueError(f"revindex moltype is {my_moltype}, not {moltype}")

        if picklist is not None:
            raise Exception("cannot use picklists, sry")

        return self

    def search(self, query, *args, **kwargs):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments:
          * do_containment: default False. If True, use Jaccard containment.
          * ignore_abundance: default False. If True, and query signature
            and database support k-mer abundances, ignore those abundances.

        Note, the "best only" hint is ignored by LCA_Database
        """
        if not query.minhash:
            return []

        print('XYY', query.minhash.scaled, self._scaled)

        # check arguments
        if "threshold" not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs["threshold"]
        do_containment = kwargs.get("do_containment", False)
        ignore_abundance = kwargs.get("ignore_abundance", False)

        self._init_inner()

        size = ffi.new("uintptr_t *")
        results_ptr = self._methodcall(
            lib.revindex_search,
            query._get_objptr(),
            threshold,
            do_containment,
            ignore_abundance,
            size,
        )

        size = size[0]
        if size == 0:
            return []

        results = []
        for i in range(size):
            match = SearchResult._from_objptr(results_ptr[i])
            if match.score >= threshold:
                results.append(
                    IndexSearchResult(match.score, match.signature, match.location)
                )

        return results

    @property
    def scaled(self):
        self._init_inner()
        return self._methodcall(lib.revindex_scaled)

    def prefetch(self, query_ss, threshold_bp=0, **kwargs):
        query_mh = query_ss.minhash
        if not query_mh:
            raise ValueError
        query_mh = query_mh.downsample(scaled=self._scaled)
        ss = SourmashSignature(query_mh)
        threshold = threshold_bp / query_mh.scaled / len(query_mh)
        print('XZZY prefetch', query_mh.scaled, self._scaled, threshold)
        sr = self.search(ss, threshold=threshold, do_containment=True)
        print(f'found: {len(sr)}')
        return sr

    def best_containment(self, query_ss, *, threshold_bp=0, **kwargs):
        query_mh = query_ss.minhash
        if not query_mh:
            raise ValueError("empty query")
        threshold = threshold_bp / query_mh.scaled / len(query_mh)
        print('XZX', query_mh.scaled, self._scaled, threshold)
        results = self.search(query_ss, threshold=threshold, do_containment=True)

        if results:
            results.sort(key=lambda x: -x.score)
            return results[0]
        raise ValueError("no results")

    def peek(self, query_mh, *, threshold_bp=0):
        if not len(query_mh):
            raise ValueError
        threshold = threshold_bp / query_mh.scaled / len(query_mh)
        print('XZZ peek', query_mh.scaled, self._scaled, threshold)
        query_ss = sourmash.SourmashSignature(query_mh)
        found = self.search(query_ss, threshold=threshold, do_containment=True)

        if found:
            match_mh = found[0].signature.minhash.flatten()
            intersect_mh = match_mh.intersection(query_mh.flatten())
            return found[0], intersect_mh
        return []

    def consume(self, intersect_mh):
        pass

    def counter_gather(self, query, threshold_bp, **kwargs):
        # will raise ValueError if empty:
        self._init_inner()

        counter = RevIndex_CounterGather(query, self, threshold_bp)
        #for result in self.prefetch(query, threshold_bp=threshold_bp):
        #    counter.add(result.signature)

        return counter


class SearchResult(RustObject):
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


class DiskRevIndex_DatasetPicklist(RustObject):
    __dealloc_func__ = lib.dataset_picklist_free

    def __init__(self, idxs):
        idx_list = list(idxs)
        idx_list_size = len(idx_list)

        self._objptr = rustcall(
            lib.dataset_picklist_new_from_list, idx_list, idx_list_size
        )


class DiskRevIndex(RustObject, Index):
    """
    RocksDB-based low-memory on disk inverted index, implemented in Rust.
    """

    __dealloc_func__ = lib.disk_revindex_free
    is_database = True
    manifest = None

    def __init__(self, path):
        check_file = os.path.join(path, "CURRENT")
        if not os.path.exists(check_file):
            raise ValueError("not a RocksDB")

        # create via FFI
        path_b = path.encode("utf-8")
        self._objptr = rustcall(lib.disk_revindex_new_from_rocksdb, path_b)

        # store location
        self._path = path
        self._idx_picklist = None

    @property
    def location(self):
        return self._path

    def insert(self, *args, **kwargs):
        raise NotImplementedError

    def load(self, *args, **kwargs):  # @CTB
        raise NotImplementedError

    def save(self, *args, **kwargs):
        raise NotImplementedError

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
        if 0 and moltype is not None and moltype != my_moltype:  #  @CTB
            raise ValueError(f"revindex moltype is {my_moltype}, not {moltype}")

        if picklist is not None:
            # @CTB building manifest this way is expensive!!
            m = CollectionManifest.create_manifest(
                self._signatures_with_internal(), include_signature=False
            )
            m = m.select_to_manifest(picklist=picklist)
            self._generate_idx_picklist_from_manifest(m)

        return self

    def _generate_idx_picklist_from_manifest(self, mf):
        if self._idx_picklist is not None:
            raise Exception("cannot use picklists multiple times, sorry")

        # grab internal indices
        idx_list = [int(row["internal_location"]) for row in mf.rows]
        self._idx_picklist = DiskRevIndex_DatasetPicklist(idx_list)

    @property
    def _ffi_idx_picklist(self):
        if self._idx_picklist is None:
            return ffi.NULL
        return self._idx_picklist._objptr

    def signatures(self):  # @CTB add picklist
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
        # CTB fix: don't use signatures() once we start paying attention
        # to picklists.
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
            self._ffi_idx_picklist,
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
                self._ffi_idx_picklist,
            )
        elif do_max_containment:
            raise NotImplementedError(
                "max_containment is not (yet) available on RocksDB"
            )
        else:  # jaccard
            results_ptr = self._methodcall(
                lib.disk_revindex_search_jaccard,
                query_ss._get_objptr(),
                threshold,
                size,
                self._ffi_idx_picklist,
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
                lib.disk_revindex_best_containment,
                query_ss._get_objptr(),
                threshold_bp,
                self._ffi_idx_picklist,
            )
            match_ss = SourmashSignature._from_objptr(ss_ptr)
            if not match_ss.minhash:
                raise ValueError("no results")
        except:
            raise ValueError("no results")
        containment = query_ss.contained_by(match_ss)

        return IndexSearchResult(containment, match_ss, self.location)

    #
    # implement CounterGather API
    #

    def peek(self, query_mh, *, threshold_bp=0):
        ss_ptr = self._methodcall(
            lib.disk_revindex_peek,
            query_mh._get_objptr(),
            int(threshold_bp),
            self._ffi_idx_picklist,
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
        counter = RevIndex_CounterGather(query, self, threshold_bp)
        for result in self.prefetch(query, threshold_bp=threshold_bp):
            counter.add(result.signature)

        return counter


class RevIndex_CounterGather:
    """
    Simple implementation of CounterGather API that tracks matches
    while passing most calls back to the parent RevIndex.
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
        self.locations = dict()

    def add(self, match_ss, *, location=None, require_overlap=True): # @CTB location
        if self.allow_insert:
            x = self.db._check_not_init(do_raise=False)
            print('checking', x)
            if self.db._check_not_init(do_raise=False):
                print('not init')
                self.db.insert(match_ss)
            else:
                raise ValueError

            self.locations[match_ss.md5sum()] = location

        query_mh = self.orig_query_mh
        match_mh = match_ss.minhash.downsample(scaled=query_mh.scaled)
        intersect_mh = query_mh.intersection(match_mh.flatten())
        if require_overlap and not intersect_mh:
            raise ValueError("require overlap")

        self.found_mh += intersect_mh

    def peek(self, query_mh, *, threshold_bp=0): # threshold_bp default?? @CTB
        if not query_mh:
            return []

        if query_mh.contained_by(self.orig_query_mh) != 1.0:
            raise ValueError
        #assert threshold_bp is not None
        print('BBB peek')

        res = self.db.peek(query_mh, threshold_bp=threshold_bp)
        if not res:
            return []

        sr, intersect_mh = res
        sr_ss = sr.signature
        sr_score = sr.score
        new_sr = IndexSearchResult(sr_score, sr_ss,
                                   self.locations[sr_ss.md5sum()])
        return new_sr, intersect_mh

    def consume(self, intersect_mh):
        self.db._init_inner()
        self.found_mh += intersect_mh

    @property
    def union_found(self):
        return self.found_mh

    def signatures(self):
        # don't track actual signatures - go back to RevIndex
        for sr in self.db.prefetch(self.query):
            print('FOO')
            yield sr.signature
