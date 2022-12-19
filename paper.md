---
title: 'sourmash: a tool to quickly search, compare, and analyze genomic and metagenomic data sets'
tags:
  - FracMinHash
  - MinHash
  - k-mers
  - Python
  - Rust
authors:
 - name: Luiz Irber
   orcid: 0000-0003-4371-9659
   affiliation: University of California, Davis
 - name: N. Tessa Pierce-Ward
   orcid: 0000-0002-2942-5331
   affiliation: University of California, Davis
 - name: Tim Head
   orcid: 0000-0003-0931-3698
   affiliation: @@@
 - name: Harriet Alexander
   orcid: 0000-0003-1308-8008
   affiliation: Woods Hole Oceanographic Institute
 - name: Olga Botvinnik
   orcid: 0000-0003-4412-7970
   affiliation: Chan-Zuckerberg Biohub
 - name: Phillip Brooks
   orcid: 0000-0003-3987-244X
   affiliation: 10x Genomics
 - name: Laurent Gautier
   orcid: 0000-0003-0638-3391
   affiliation: @@@
 - name: Lisa Johnson
   orcid: 0000-0002-3600-7218
   affiliation: 10x Genomics
 - name: Fabian Klötzl
   orcid: @@@
   affiliation: @@@
 - name: David Koslicki
   orcid: 0000-0002-0640-954X
   affiliation: Pennsylvania State University
 - name: Mahmudur Rahman Hera
   orcid: 0000-0002-5992-9012
   affiliation: Pennsylvania State University
 - name: Katrin Leinweber
   orcid: 0000-0001-5135-5758
   affiliation: @@@
 - name: Ivan Ogasawara
   orcid: @@@
   affiliation: @@@
 - name: Taylor Reiter
   orcid: 0000-0002-7388-421X
   affiliation: University of California, Davis
 - name: Camille Scott
   orcid: 0000-0001-8822-8779
   affiliation: University of California, Davis
 - name: Connor Skennerton
   orcid: 0000-0003-1320-4873
   affiliation: @@@
 - name: Daniel Standage
   orcid: 0000-0003-0342-8531.
   affiliation: @@@
 - name: Joshua Swamidass
   orcid: 0000-0003-2191-0778
   affiliation: @@@
 - name: Connor Tiffany
   orcid: @@@
   affiliation: @@@
 - name: Erik Young
   orcid: @@@
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
search, comparison, and analysis [@Pierce:2019]. sourmash's FracMinHash sketching supports fast and accurate sequence comparisons between datasets of different sizes [@gather], including petabase-scale database search [@branchwater]. From release 4.x, sourmash is built on top of Rust and provides an experimental Rust interface.

FracMinHash sketching is lossy compression approach that represents
data sets using a "fractional" sketch containing $1/S$ of the original 
k-mers. Like other sequence sketching techniques (e.g. MinHash, [@Ondov:2015]), FracMinHash provides a lightweight way to store representations of large DNA or RNA sequence collections for comparison and search. Sketches can be used to identify samples, find similar samples, identify data sets with shared sequences, and build phylogenetic trees. FracMinHash sketches support estimation of overlap, bidirectional containment, and Jaccard similarity between data sets and is accurate even for data sets of different sizes.

Since sourmash v1 was released in 2016 [@Brown:2016], sourmash has expanded
to support new database types and many more command line functions.
In particular, sourmash now has robust support for both Jaccard similarity
and containment calculations, which enables analysis and comparison of data sets
of different sizes, including large metagenomic samples. As of v4.4,
sourmash can convert these to estimated Average Nucleotide Identity (ANI)
values, which can provide improved biological context to sketch comparisons [@hera2022debiasing].

# Statement of Need

Large collections of genomes, transcriptomes, and raw sequencing data
sets are readily available in biology, and the field needs lightweight
computational methods for searching and summarizing the content of
both public and private collections. sourmash provides a flexible set
of programmatic functionality for this purpose, together with a robust
and well-tested command-line interface. It has been used in well over 200
publications (based on citations of @Brown:2016 and @Pierce:2019) and it continues
to expand in functionality.

# Acknowledgements

This work is funded in part by the Gordon and Betty Moore Foundation’s
Data-Driven Discovery Initiative [GBMF4551 to CTB].

# References
