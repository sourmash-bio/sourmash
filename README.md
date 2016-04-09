# sourmash

Compute MinHash signatures for RNAseq reads.

Usage::

   ./sourmash compute *.fq.gz
   ./sourmash compare *.sig -o distances
   ./plot-comparison.py distances
