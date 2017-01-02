Additional information on sourmash
==================================

Other MinHash implementations for DNA
-------------------------------------

In addition to `mash <https://github.com/marbl/Mash>`__, also see
`RKMH: Read Classification by Kmers
<https://github.com/edawson/rkmh>`__.

Blog posts
----------

We have a number of blog posts on sourmash and MinHash more generally:

* `Applying MinHash to cluster RNAseq samples
  <http://ivory.idyll.org/blog/2016-sourmash.html>`__

* `MinHash signatures as ways to find samples, and collaborators?
  <http://ivory.idyll.org/blog/2016-sourmash-signatures.html>`__

* `Efficiently searching MinHash Sketch collections
  <http://ivory.idyll.org/blog/2016-sourmash-sbt.html>`__ - indexing and
  search 42,000 bacterial genomes with Sequence Bloom Trees.

* `Quickly searching all the microbial genomes, mark 2 - now with
  archaea, phage, fungi, and protists!
  <http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html>`__ - indexing
  and searching 50,000 microbial genomes, round 2.

* `What metadata should we put in MinHash Sketch signatures?
  <http://ivory.idyll.org/blog/2016-sourmash-signatures-metadata.html>`__ -
  crowdsourcing ideas for what metadata belongs in a signature file.

* `Minhashing all the things (part 1): microbial genomes
  <http://blog.luizirber.org/2016/12/28/soursigs-arch-1/>`__ - on
  approaches to computing MinHashes for large collections of public data.

JSON format for the signature
-----------------------------

The JSON format is not necessarily final; this is a TODO item for future
releases.  In particular, we'd like to update it to store more metadata
for samples.

Interoperability with mash
--------------------------

The default sketches computed by sourmash and mash are comparable, but
we are still `working on ways to convert the file formats <https://github.com/marbl/Mash/issues/27>`__.

Developing sourmash
-------------------

Please see:

.. toctree::
   :maxdepth: 2

   developer
   release
