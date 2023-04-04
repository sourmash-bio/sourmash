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
   equal-contrib: true
   affiliation: 1
 - name: N. Tessa Pierce-Ward
   orcid: 0000-0002-2942-5331
   equal-contrib: true
   affiliation: 1
 - name: Mohamed Abuelanin
   orcid: 0000-0002-3419-4785
   affiliation: 1
 - name: Harriet Alexander
   orcid: 0000-0003-1308-8008
   affiliation: 2
 - name: Abhishek Anant
   orcid:  0000-0002-5751-2010
   affiliation: 9
 - name: Keya Barve
   orcid: 0000-0003-3241-2117
   affiliation: 1
 - name: Colton Baumler
   orcid: 0000-0002-5926-7792
   affiliation: 1
 - name: Olga Botvinnik
   orcid: 0000-0003-4412-7970
   affiliation: 3
 - name: Phillip Brooks
   orcid: 0000-0003-3987-244X
   affiliation: 1
 - name: Daniel Dsouza
   orcid: 0000-0001-7843-8596
   affiliation: 9
 - name: Laurent Gautier
   orcid: 0000-0003-0638-3391
   affiliation: 9
 - name: Mahmudur Rahman Hera
   orcid: 0000-0002-5992-9012
   affiliation: 4
 - name: Hannah Eve Houts
   orcid: 0000-0002-7954-4793
   affiliation: 1
 - name: Lisa K. Johnson
   orcid: 0000-0002-3600-7218
   affiliation: 1
 - name: Fabian Klötzl
   orcid: 0000-0002-6930-0592
   affiliation: 5
 - name: David Koslicki
   orcid: 0000-0002-0640-954X
   affiliation: 4
 - name: Marisa Lim
   orcid: 0000-0003-2097-8818
   affiliation: 1
 - name: Ricky Lim
   orcid: 0000-0003-1313-7076
   affiliation: 9
 - name: Ivan Ogasawara
   orcid: 0000-0001-5049-4289
   affiliation: 9
 - name: Taylor Reiter
   orcid: 0000-0002-7388-421X
   affiliation: 1
 - name: Camille Scott
   orcid: 0000-0001-8822-8779
   affiliation: 1
 - name: Andreas Sjödin
   orcid: 0000-0001-5350-4219
   affiliation: 6
 - name: Daniel Standage
   orcid: 0000-0003-0342-8531
   affiliation: 7
 - name: S. Joshua Swamidass
   orcid: 0000-0003-2191-0778
   affiliation: 8
 - name: Connor Tiffany
   orcid: 0000-0001-8188-7720
   affiliation: 9
 - name: Pranathi Vemuri
   orcid: 0000-0002-5748-9594
   affiliation: 3
 - name: Erik Young
   orcid: 0000-0002-9195-9801
   affiliation: 1
 - name: C. Titus Brown
   orcid: 0000-0001-6001-2677
   corresponding: true
   affiliation: 1
affiliations:
 - name: University of California, Davis
   index: 1
 - name: Woods Hole Oceanographic Institute
   index: 2
 - name: Chan-Zuckerberg Biohub
   index: 3
 - name:  Pennsylvania State University
   index: 4
 - name:  MPI for Evolutionary Biology
   index: 5
 - name: Swedish Defence Research Agency (FOI)
   index: 6
 - name: National Bioforensic Analysis Center
   index: 7
 - name: Washington University in St Louis
   index: 8
 - name: No affiliation
   index: 9

date: 27 Mar 2023
bibliography: paper.bib
---

# Summary

sourmash is a command line tool and Python library for sketching
collections of DNA, RNA, and amino acid k-mers for biological sequence
search, comparison, and analysis [@Pierce:2019]. sourmash's FracMinHash sketching supports fast and accurate sequence comparisons between datasets of different sizes [@gather], including petabase-scale database search [@branchwater]. From release 4.x, sourmash is built on top of Rust and provides an experimental Rust interface.

FracMinHash sketching is a lossy compression approach that represents
data sets using a "fractional" sketch containing $1/S$ of the original 
k-mers. Like other sequence sketching techniques (e.g. MinHash, [@Ondov:2015]), FracMinHash provides a lightweight way to store representations of large DNA or RNA sequence collections for comparison and search. Sketches can be used to identify samples, find similar samples, identify data sets with shared sequences, and build phylogenetic trees. FracMinHash sketching supports estimation of overlap, bidirectional containment, and Jaccard similarity between data sets and is accurate even for data sets of very different sizes.

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
