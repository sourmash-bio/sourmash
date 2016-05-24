# sourmash


[![Build Status](https://drone.io/github.com/dib-lab/sourmash/status.png)](https://drone.io/github.com/dib-lab/sourmash/latest)
[![Coverage Status](https://coveralls.io/repos/github/dib-lab/sourmash/badge.svg?branch=master)](https://coveralls.io/github/dib-lab/sourmash?branch=master)

Compute MinHash signatures for RNAseq reads.

Usage:

    ./sourmash compute *.fq.gz
    ./sourmash compare *.sig -o distances
    ./plot-comparison.py distances

This is an exploratory project; it's not really ready for people to use yet.
Buyers Beware.

We have demo notebooks on binder:

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/dib-lab/sourmash)

----

The name comes from riffing off of https://github.com/marbl/Mash,
combined with my love of whiskey.  (Sour mash is used in making
whiskey.)

## Installation

You can do:

    pip install sourmash

It currently requires khmer and PyYAML, and runs under both
Python 2.7.11 and Python 3.5.

## Development

`sourmash` is the main command-line entry point; run it for help.

`sourmash_lib.py` contains the MinHash sketch implementation.

`sourmash_signature.py` contains the YAML sig formatting and I/O functions.

`sourmash_fig.py` contains some plotting functionality.

----

CTB

12.apr.2016
