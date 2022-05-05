# Misc utility scripts

## Misc scripts

* trim-noV.sh - a script to do trimming of short reads. requires khmer >= 2.0.

## Debugging and testing scripts

* check-tree.py - exhaustively confirm the results of a search on an SBT.
* compute-dna-mh-another-way.py - a separate implementation of MinHash signature calculation for DNA.
* compute-input-prot-another-way.py - a separate implementation of MinHash signature calculation for proteins.
* compute-prot-mh-another-way.py - a separate implementation of MinHash signature computing for 6-frame translations of DNA into amino acid space.

CTB 1/2019

## Formula implementations

* cardinality_estimate_confidence.py - a function that will tell you if the sketch size is too small to trust the estimate `sketch_size * scale` as an estimate of the number of distinct k-mers.

DMK 5/2022
