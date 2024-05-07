# sourmash

Quickly search, compare, and analyze genomic and metagenomic data sets.

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<a href="https://github.com/sourmash-bio/sourmash/blob/latest/LICENSE"><img alt="License: 3-Clause BSD" src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg"></a>
[![Documentation](https://readthedocs.org/projects/sourmash/badge/?version=latest)](http://sourmash.readthedocs.io/en/latest/)
[![Gitter](https://badges.gitter.im/sourmash-bio/community.svg)](https://gitter.im/sourmash-bio/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00027/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00027)
[![pyOpenSci](https://tinyurl.com/y22nb8up)](https://github.com/pyOpenSci/software-submission/issues/129)

[![Bioconda install](https://img.shields.io/conda/dn/bioconda/sourmash.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/sourmash)
<a href="https://pypi.org/project/sourmash/"><img alt="PyPI" src="https://badge.fury.io/py/sourmash.svg"></a>
[![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/sourmash-minimal.svg)](https://anaconda.org/conda-forge/sourmash-minimal)

![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)
[![Build Status](https://github.com/sourmash-bio/sourmash/workflows/Python%20tests/badge.svg)](https://github.com/sourmash-bio/sourmash/actions/)
[![codecov](https://codecov.io/gh/sourmash-bio/sourmash/branch/latest/graph/badge.svg)](https://codecov.io/gh/sourmash-bio/sourmash)

<p align="center"><img src="https://raw.githubusercontent.com/sourmash-bio/sourmash/latest/doc/_static/logo.png" height="256" /></p>

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

sourmash is a k-mer analysis multitool, and we aim to provide stable, robust programmatic and command-line APIs for a variety of sequence comparisons. Some of our special sauce includes:
- `FracMinHash` sketching, which enables accurate comparisons (including ANI) between data sets of different sizes
- `sourmash gather`, a combinatorial k-mer approach for more accurate metagenomic profiling

Please see the [sourmash publications](https://sourmash.readthedocs.io/en/latest/publications.html#sourmash-fundamentals) for details.

The name is a riff off of [Mash](https://github.com/marbl/Mash),
combined with @ctb's love of whiskey.
([Sour mash](https://en.wikipedia.org/wiki/Sour_mash) is used in
making whiskey.)

Maintainers: [C. Titus Brown](mailto:titus@idyll.org) ([@ctb](http://github.com/ctb)), [Luiz C. Irber, Jr](mailto:luiz@sourmash.bio) ([@luizirber](http://github.com/luizirber)), and [N. Tessa Pierce-Ward](mailto:tessa@sourmash.bio) ([@bluegenes](http://github.com/bluegenes)).

sourmash was initially developed by the
[Lab for Data-Intensive Biology](http://ivory.idyll.org/lab/) at the
[UC Davis School of Veterinary Medicine](http://www.vetmed.ucdavis.edu),
and now includes contributions from the global research and developer
community.

## Installation

We recommend using conda-forge to install sourmash:

```
conda install -c conda-forge sourmash-minimal
```
This will install the latest stable version of sourmash 4.

You can also use pip to install sourmash:

```
pip install sourmash
```

A quickstart tutorial [is available](https://sourmash.readthedocs.io/en/latest/tutorials.html).

### Requirements

sourmash runs under Python 3.10 and later on Windows, Mac OS X, and
Linux.  The base requirements are screed, cffi, numpy, matplotlib, and
scipy.  Conda will install everything necessary, and is
our recommended installation method (see below).

### Installation with conda

conda-forge is a community maintained channel for the
[conda](http://conda.pydata.org/docs/intro.html) package manager.
[installing conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/),
you can install sourmash by running:

```bash
$ conda create -n sourmash_env -c conda-forge sourmash-minimal
$ conda activate sourmash_env
$ sourmash --help
```

which will install
[the latest released version](https://github.com/sourmash-bio/sourmash/releases).

## Support

For questions, please open an issue [on Github](https://github.com/sourmash-bio/sourmash/issues), or ask in our [chat](https://gitter.im/sourmash-bio/community?utm_source=share-link&utm_medium=link&utm_campaign=share-link).

## Development

Development happens on github at
[sourmash-bio/sourmash](https://github.com/sourmash-bio/sourmash).

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

CTB
Jan 2024
