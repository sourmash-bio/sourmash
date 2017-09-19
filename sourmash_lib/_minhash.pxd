# -*- coding: UTF-8 -*-
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii

from __future__ import unicode_literals

from libcpp cimport bool
from libcpp.map cimport map
from libcpp.memory cimport unique_ptr
from libcpp.set cimport set as cppset
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t
from libcpp.vector cimport vector


cdef extern from "kmer_min_hash.hh":
    ctypedef uint64_t HashIntoType;
    ctypedef vector[HashIntoType] CMinHashType;


    cdef uint64_t _hash_murmur(const string, uint32_t seed)


    cdef cppclass KmerMinHash:
        const uint32_t seed;
        const unsigned int num;
        const unsigned int ksize;
        const bool is_protein;
        const HashIntoType max_hash;
        CMinHashType mins;

        KmerMinHash(unsigned int, unsigned int, bool, uint32_t, HashIntoType)
        void add_hash(HashIntoType) except +ValueError
        void add_word(string word) except +ValueError
        void add_sequence(const char *, bool) except +ValueError
        void merge(const KmerMinHash&) except +ValueError
        unsigned int count_common(const KmerMinHash&) except +ValueError
        unsigned long size()


    cdef cppclass KmerMinAbundance(KmerMinHash):
        CMinHashType abunds;

        KmerMinAbundance(unsigned int, unsigned int, bool, uint32_t, HashIntoType)
        void add_hash(HashIntoType) except +ValueError
        void add_word(string word) except +ValueError
        void add_sequence(const char *, bool) except +ValueError
        void merge(const KmerMinAbundance&) except +ValueError
        void merge(const KmerMinHash&) except +ValueError
        unsigned int count_common(const KmerMinAbundance&) except +ValueError
        unsigned long size()


cdef class MinHash(object):
    cdef unique_ptr[KmerMinHash] _this
    cdef public bool track_abundance

    cpdef get_mins(self, bool with_abundance=*)
    cpdef set_abundances(self, dict)
