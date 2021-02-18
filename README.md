# sourmash

Quickly search, compare, and analyze genomic and metagenomic data sets.

[![Documentation](https://readthedocs.org/projects/sourmash/badge/?version=latest)](http://sourmash.readthedocs.io/en/latest/)
[![Build Status](https://github.com/dib-lab/sourmash/workflows/Python%20tests/badge.svg)](https://github.com/dib-lab/sourmash/actions/)
[![Bioconda install](https://img.shields.io/conda/dn/bioconda/sourmash.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/sourmash)
<a href="https://pypi.org/project/sourmash/"><img alt="PyPI" src="https://badge.fury.io/py/sourmash.svg"></a>
[![codecov](https://codecov.io/gh/dib-lab/sourmash/branch/latest/graph/badge.svg)](https://codecov.io/gh/dib-lab/sourmash)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00027/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00027)
<a href="https://github.com/dib-lab/sourmash/blob/latest/LICENSE"><img alt="License: 3-Clause BSD" src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg"></a>

<p align="center"><img src="https://raw.githubusercontent.com/dib-lab/sourmash/latest/doc/_static/logo.png" height="256" /></p>

Usage:

    sourmash sketch dna *.fq.gz
    sourmash compare *.sig -o distances.cmp -k 31
    sourmash plot distances.cmp

sourmash 1.0 is [published on JOSS](https://doi.org/10.21105/joss.00027); please cite that paper if you use sourmash (`doi: 10.21105/joss.00027`):.

The latest major release is sourmash v4, which has several
command-line and Python incompatibilities with previous
versions. Please
[visit our migration guide](https://sourmash.readthedocs.io/en/latest/support.html#migrating-from-sourmash-v3-x-to-sourmash-4-x)
to upgrade!

----

The name is a riff off of [Mash](https://github.com/marbl/Mash),
combined with @ctb's love of whiskey.
([Sour mash](https://en.wikipedia.org/wiki/Sour_mash) is used in
making whiskey.)

Primary authors: [C. Titus Brown](mailto:titus@idyll.org) ([@ctb](http://github.com/ctb)) and [Luiz C. Irber, Jr](mailto:sourmash@luizirber.org) ([@luizirber](http://github.com/luizirber)).

sourmash was initially developed by the
[Lab for Data-Intensive Biology](http://ivory.idyll.org/lab/) at the
[UC Davis School of Veterinary Medicine](http://www.vetmed.ucdavis.edu),
and now includes contributions from the global research and developer
community.

## Installation

We recommend using bioconda to install sourmash:

```
conda install -c conda-forge -c bioconda sourmash
```
This will install the latest stable version of sourmash 4.

You can also use pip to install sourmash:

```
pip install sourmash
```

A quickstart tutorial [is available](https://sourmash.readthedocs.io/en/latest/tutorials.html).

### Requirements

sourmash runs under Python 3.7 and later.  The base
requirements are screed, cffi, numpy, matplotlib, and scipy.  Conda
(see below) will install everything necessary, and is our recommended
installation method.

### Installation with conda

Bioconda is a channel for the
[conda](http://conda.pydata.org/docs/intro.html) package manager with
a focus on bioinformatics software. After
[installing conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/),
you can install sourmash by running:

```bash
$ conda create -n sourmash_env -c conda-forge -c bioconda sourmash python=3.7
$ source activate sourmash_env
$ sourmash --help
```

which will install
[the latest released version](https://github.com/dib-lab/sourmash/releases).

## Support

Please ask questions and files issues
[on Github](https://github.com/dib-lab/sourmash/issues).

## Development

Development happens on github at
[dib-lab/sourmash](https://github.com/dib-lab/sourmash).

sourmash is developed in Python and Rust, and you will need a Rust
environment to build it; see [the developer notes](doc/developer.md)
for our suggested development setup.

After installation, `sourmash` is the main command-line entry point;
run it with `python -m sourmash`, or do `pip install -e /path/to/repo` to
do a developer install in a virtual environment.

The `sourmash/` directory contains the Python library and command-line interface code.

The `src/core/` directory contains the Rust library implementing core
functionality.

Tests require py.test and can be run with `make test`.

Please see [the developer notes](doc/developer.md) for more information
on getting set up with a development environment.

----

CTB
Feb 2021
