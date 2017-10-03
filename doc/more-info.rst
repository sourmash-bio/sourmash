Additional information on sourmash
==================================

Other MinHash implementations for DNA
-------------------------------------

In addition to `mash <https://github.com/marbl/Mash>`__, also see:

* `RKMH: Read Classification by Kmers <https://github.com/edawson/rkmh>`__.
* `mashtree <https://github.com/lskatz/mashtree/blob/master/README.md>`__ for building trees using Mash distances.
* `Finch: a Mash implementation in Rust
  <https://github.com/onecodex/finch-rs>`__. Quote, "Fast sketches,
  count histograms, better filtering."
* `BBMap and SendSketch <http://seqanswers.com/forums/showthread.php?t=74019>`__ - part of Brian Bushnell's tool collection.

If you are interested in exactly how these MinHash approaches
calculate the hashes of DNA sequences, please see some simple Python
code in sourmash, `utils/compute-dna-mh-another-way.py
<https://github.com/dib-lab/sourmash/blob/master/utils/compute-dna-mh-another-way.py>`__.

Papers and references
---------------------

`On the resemblance and containment of documents <http://ieeexplore.ieee.org/document/666900/?reload=true>`__, Broder, 1997. The original MinHash paper!

`Mash: fast genome and metagenome distance estimation using MinHash. <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`__, Ondov et al. 2016.

`sourmash: a library for MinHash sketching of DNA. <http://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9>`__, Brown and Irber, 2017.

`Improving MinHash via the Containment Index with Applications to Metagenomic Analysis <https://www.biorxiv.org/content/early/2017/09/04/184150>`__, Koslicki and Zabeti, 2017.

Presentations and posters
-------------------------

`Taxonomic classification of microbial metagenomes using MinHash signatures <https://osf.io/mu4gk/>`__, Brooks et al., 2017. Presented at ASM.

Blog posts
----------

We (and others) have a number of blog posts on sourmash and MinHash
more generally:

* `Some background on k-mers, and their use in taxonomy <http://ivory.idyll.org/blog/2017-something-about-kmers.html>`__

* From the Phillippy lab: `mash screen: what's in my sequencing run? <https://genomeinformatics.github.io/mash-screen/>`__

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

* `Comparing genome sets extracted from metagenomes <http://ivory.idyll.org/blog/2017-comparing-genomes-from-metagenomes.html>`__.

* `Taxonomic examinations of genome bins from Tara Oceans <http://ivory.idyll.org/blog/2017-taxonomy-of-tara-ocean-genomes.html>`__.

* `Classifying genome bins using a custom reference database, part I <http://ivory.idyll.org/blog/2017-classify-genome-bins-with-custom-db-part-1.html>`__.

* `Classifying genome bins using a custom reference database, part II <http://ivory.idyll.org/blog/2017-classify-genome-bins-with-custom-db-part-2.html>`__.

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
