# sourmash

Compute MinHash signatures for RNAseq reads.

Usage:

    ./sourmash compute *.fq.gz
    ./sourmash compare *.sig -o distances
    ./plot-comparison.py distances

This is an exploratory project; it's not really ready for people to use yet.
Buyers Beware.

----

The name comes from riffing off of https://github.com/marbl/Mash,
combined with my love of whiskey.  (Sour mash is used in making
whiskey.)

## Installation

You can do:

    pip install sourmash

It currently requires khmer and PyYAML, and runs under both
Python 2.7.11 and Python 3.5.

CTB

10.apr.2016
