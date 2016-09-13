---
title: 'sourmash: a library for MinHash sketching of DNA'
tags:
  - MinHash
  - k-mers
  - Python
authors:
 - name: C. Titus Brown
   orcid: 0000-0001-6001-2677
   affiliation: University of California, Davis
 - name: Luiz Irber
   orcid: 0000-0003-4371-9659
   affiliation: University of California, Davis
date: 13 Sep 2016
bibliography: paper.bib
---

# Summary

sourmash is a toolbox for creating, comparing, and manipulating MinHash
sketches of genomic data.

MinHash sketches provide a lightweight way to store "signatures" of
large DNA or RNA sequence collections, and then compare or search them
using a Jaccard index.  MinHash sketches can be used to identify samples,
find similar samples, identify data sets with shared sequences, and
build phylogenetic trees [@ondov2015fast].

sourmash provides a command line script, a Python library, and a CPython
module for MinHash sketches.

# References
