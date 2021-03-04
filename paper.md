---
title: 'sourmash: a tool to quickly search, compare, and analyze genomic and metagenomic data sets'
tags:
  - MinHash
  - k-mers
  - Python
  - Rust
authors:
 - name: Luiz Irber
   orcid: 0000-0003-4371-9659
   affiliation: University of California, Davis
 - name: C. Titus Brown
   orcid: 0000-0001-6001-2677
   affiliation: University of California, Davis
date: 3 Mar 2021
bibliography: paper.bib
---

# Summary

sourmash is a command line tool and Python library for sketching
collections of DNA, RNA, and amino acid k-mers for biological sequence
search, comparison, and analysis [@Pierce2019]. The most recent
release, sourmash v4, is built on top of Rust and provides an
experimental Rust interface.

MinHash sketches provide a lightweight way to store "signatures" of
large DNA or RNA sequence collections, and then compare or search them
using a Jaccard index [@ondov2015fast].  MinHash sketches can be used
to identify samples, find similar samples, identify data sets with
shared sequences, and build phylogenetic trees

Since sourmash v1 was released in 2016 [@Brown2016], sourmash has expanded
to support new database types and many more command line functions.
In particular, sourmash now has robust support for both Jaccard similarity
and containment calculations, which supports robust analysis and comparison
of large metagenomic samples.

# References
