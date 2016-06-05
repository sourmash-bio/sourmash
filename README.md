# sourmash

[![Documentation](https://readthedocs.org/projects/sourmash/badge/?version=latest)](http://sourmash.readthedocs.io/en/latest/)
[![Build Status](https://drone.io/github.com/dib-lab/sourmash/status.png)](https://drone.io/github.com/dib-lab/sourmash/latest)
[![codecov](https://codecov.io/gh/dib-lab/sourmash/branch/master/graph/badge.svg)](https://codecov.io/gh/dib-lab/sourmash)

Compute MinHash signatures for DNA sequences.

Usage:

    ./sourmash compute *.fq.gz
    ./sourmash compare *.sig -o distances
    ./plot-comparison.py distances

We have demo notebooks on binder:

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/dib-lab/sourmash)

----

The name is a riff off of [Mash](https://github.com/marbl/Mash), combined with
my love of whiskey.  (Sour mash is used in making whiskey.)

## Installation

You can do:

    pip install sourmash

sourmash runs under both Python 2.7.x and Python 3.5.  The base
requirements are screed and PyYAML, together with a C++ development
environment and the CPython development headers and libraries (for the
C++ extension).

The comparison code (`sourmash compare`) uses numpy, and the plotting
code uses matplotlib and scipy, but most of the code is usable without
these.

## Development

`sourmash` is the main command-line entry point; run it for help.

`sourmash_lib/` contains the library code.

Tests require py.test and can be run with `make test`.

----

CTB

6.jun.2016
