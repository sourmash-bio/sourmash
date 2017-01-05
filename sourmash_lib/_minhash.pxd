# -*- coding: UTF-8 -*-
# cython: c_string_type=str, c_string_encoding=ascii

from __future__ import unicode_literals

from libcpp cimport bool
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t


cdef extern from "kmer_min_hash.hh":
    ctypedef uint64_t HashIntoType;
    ctypedef set[HashIntoType] CMinHashType;
    ctypedef map[HashIntoType, uint64_t] CMinAbundanceType;


    cdef uint64_t _hash_murmur(const string, uint32_t seed)


    cdef cppclass KmerMinHash:
        #const uint32_t seed;
        const unsigned int num;
        const unsigned int ksize;
        const bool is_protein;
        CMinHashType mins;

        KmerMinHash(unsigned int, unsigned int, bool, uint32_t)
        void _shrink()
        void add_hash(HashIntoType)
        void add_word(string word)
        void add_sequence(const char *, bool)
        void merge(const KmerMinHash&)
        unsigned int count_common(const KmerMinHash&) except +


    cdef cppclass KmerMinAbundance(KmerMinHash):
        CMinAbundanceType mins;

        KmerMinAbundance(unsigned int, unsigned int, bool, uint32_t)
        void add_hash(HashIntoType)
        void add_word(string word)
        void add_sequence(const char *, bool)
        void merge(const KmerMinAbundance&)
        void merge(const KmerMinHash&)
        unsigned int count_common(const KmerMinAbundance&) except +
        #unsigned int count_common(const KmerMinHash&) except +


cdef class MinHash(object):
    cdef KmerMinHash *_this
    cdef bool track_abundance

    cpdef add_sequence(self, str sequence, bool force=*)
    cpdef get_mins(self, bool with_abundance=*)
    cpdef add_hash(self, uint64_t h)
    cpdef set_abundances(self, dict)
    cpdef uint64_t count_common(self, MinHash other)
    cpdef add_protein(self, str sequence)
    #cdef __copy__()
    #cdef merge()
