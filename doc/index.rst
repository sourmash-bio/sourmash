Welcome to sourmash!
====================

sourmash is a command-line tool and Python library for computing
`hash sketches <https://en.wikipedia.org/wiki/MinHash>`__ from DNA
sequences, comparing them to each other, and plotting the results.
This allows you to estimate sequence similarity between even very
large data sets quickly and accurately.

sourmash can be used to quickly search large databases of genomes
for matches to query genomes and metagenomes; see `our list of
available databases <databases.html>`__.

sourmash also includes k-mer based taxonomic exploration and
classification routines for genome and metagenome analysis. These
routines can use the NCBI taxonomy but do not depend on it in any way.

We have `several tutorials <tutorials.html>`__ available! Start with
`Making signatures, comparing, and searching <tutorial-basic.html>`__.

The paper `Large-scale sequence comparisons with sourmash (Pierce et
al., 2019) <https://f1000research.com/articles/8-1006>`__ gives an
overview of how sourmash works and what its major use cases
are. Please also see the `mash software
<http://mash.readthedocs.io/en/latest/>`__ and the `mash paper (Ondov
et al., 2016) <http://dx.doi.org/10.1186/s13059-016-0997-x>`__ for
background information on how and why MinHash works.

**Questions? Thoughts?** Ask us on the `sourmash issue tracker <https://github.com/dib-lab/sourmash/issues/>`__!

----

To use sourmash, you must be comfortable with the UNIX command line;
programmers may find the `Python library and API <api.html>`__ useful as well.

If you use sourmash, please cite us!

   Brown and Irber (2016),
   **sourmash: a library for MinHash sketching of DNA**
   Journal of Open Source Software, 1(5), 27, `doi:10.21105/joss.00027 <https://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9>`__

sourmash in brief
-----------------

sourmash uses MinHash-style sketching to create "signatures", compressed
representations of DNA/RNA sequence.  These signatures can then
be stored, searched, explored, and taxonomically annotated.

* ``sourmash`` provides command line utilities for creating, comparing,
  and searching signatures, as well as plotting and clustering
  signatures by similarity (see `the command-line docs <command-line.html>`__).

* ``sourmash`` can **search very large collections of signatures** to find matches
  to a query.

* ``sourmash`` can also **identify parts of metagenomes that match known genomes**,
  and can **taxonomically classify genomes and metagenomes** against databases
  of known species.

* ``sourmash`` can be used to **search databases of public sequences**
  (e.g. all of GenBank) and can also be used to create and search databases
  of **private sequencing data**.

* ``sourmash`` supports saving, loading, and communication of signatures
  via `JSON <http://www.json.org/>`__, a ~human-readable &
  editable format.

* ``sourmash`` also has a simple Python API for interacting with signatures,
  including support for online updating and querying of signatures
  (see `the API docs <api.html>`__).

* ``sourmash`` isn't terribly slow, and relies on an underlying Cython
  module.

* ``sourmash`` is developed `on GitHub
  <https://github.com/dib-lab/sourmash>`__ and is **freely and openly
  available** under the BSD 3-clause license.  Please see `the README
  <https://github.com/dib-lab/sourmash/blob/master/README.md>`__ for
  more information on development, support, and contributing.

You can take a look at sourmash analyses on real data `in a saved
Jupyter notebook
<https://github.com/dib-lab/sourmash/blob/master/doc/sourmash-examples.ipynb>`__,
and experiment with it yourself `interactively in a Jupyter Notebook
<https://mybinder.org/v2/gh/dib-lab/sourmash/master?filepath=doc%2Fsourmash-examples.ipynb>`__ at `mybinder.org
<http://mybinder.org>`__.

Installing sourmash
-------------------

We currently suggest installing the latest pre-release in the sourmash
2.0 series; please see `the README file in
github.com/dib-lab/sourmash <https://github.com/dib-lab/sourmash/blob/master/README.md>`__
for information.  You can use pip or conda equally well.

Memory and speed
----------------

sourmash has relatively small disk and memory requirements compared to
many other software programs used for genome search and taxonomic
classification.

First, ``mash`` beats sourmash in speed and memory, so if you can use mash,
more power to you :)

``sourmash search`` and ``sourmash gather`` can be used to search all
genbank microbial genomes (`using our prepared
databases <databases.html>`__) with about 20 GB of disk and in under 1 GB
of RAM.  Typically a search for a single genome takes about 30 seconds
on a laptop.

``sourmash lca`` can be used to search/classify against all genbank
microbial genomes with about 200 MB of disk space and about 10 GB of
RAM. Typically a metagenome classification takes about 1 minute on a
laptop.

Limitations
-----------

**sourmash cannot find matches across large evolutionary distances.**

sourmash seems to work well to search and compare data sets for
matches at the species and genus level, but does not have much
sensitivity beyond that.  (It seems to be particularly good at
strain-level analysis.)  You should use protein-based analyses
to do searches across larger evolutionary distances.

**sourmash signatures can be very large.**

We use a modification of the MinHash sketch approach that allows us
to search the contents of metagenomes and large genomes with no loss
of sensitivity, but there is a tradeoff: there is no guaranteed limit
to signature size when using 'scaled' signatures.

Contents:
---------

.. toctree::
   :maxdepth: 2

   command-line
   tutorials
   using-sourmash-a-guide
   classifying-signatures
   api
   more-info
   support


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
