# sourmash

Compute MinHash signatures for RNAseq reads.

Usage::

   ./sourmash compute *.fq.gz
   ./sourmash compare *.sig -o distances
   ./plot-comparison.py distances

----

The name comes from riffing off of https://github.com/marbl/Mash,
combined with my love of whiskey.  (Sour mash is used in making
whiskey.)
