.. note::
    We are working on releasing sourmash 2.0,
    and this documentation reflects changes not released officially yet.
    If you want to see the documentation for sourmash 1.0 check
    the `stable version <https://sourmash.readthedocs.io/en/stable/>`__.

Welcome to sourmash!
====================

sourmash is a command-line tool and Python library for computing
`MinHash sketches <https://en.wikipedia.org/wiki/MinHash>`__ from DNA
sequences, comparing them to each other, and plotting the results.
This allows you to estimate sequence similarity between even very
large data sets quickly and accurately.

sourmash can also be used to quickly search large databases of genomes
for matches to query genomes and metagenomes; see `our list of
available databases <databases.html>`__.

sourmash also includes k-mer based taxonomic exploration and
classification routines for genome and metagenome analysis. These
routines can use the NCBI taxonomy but do not depend on it in any way.

Please see the `mash <http://mash.readthedocs.io/en/latest/>`__
software and the `mash paper (Ondov et al., 2016)
<http://biorxiv.org/content/early/2015/10/26/029827>`__ for background
information on how and why MinHash sketches work.

We have two tutorials available: one on `genome and metagenome searching <tutorial.html>`__, and one on `taxonomic analysis <tutorial-lca.html>`__.

----

To use sourmash, you must be comfortable with the UNIX command line;
programmers may find the Python library and API useful as well.

sourmash in brief
-----------------

sourmash uses MinHash sketching to create "signatures", compressed
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
<https://github.com/dib-lab/sourmash/blob/master/demo/00-demo.ipynb>`__,
and experiment with it yourself `interactively with a binder
<http://mybinder.org/repo/dib-lab/sourmash>`__ at `mybinder.org
<http://mybinder.org>`__.

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
   tutorials-lca
   classifying-signatures
   databases
   api
   api-example
   requirements
   more-info
   support


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
