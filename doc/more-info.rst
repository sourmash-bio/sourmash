Additional information on sourmash
==================================

Titus wrote a few blog posts on sourmash:

* `Applying MinHash to cluster RNAseq samples <http://ivory.idyll.org/blog/2016-sourmash.html>`__

* `MinHash signatures as ways to find samples, and collaborators? <http://ivory.idyll.org/blog/2016-sourmash-signatures.html>`__

YAML format for the signature
-----------------------------

The YAML format is not necessarily final; this is a TODO item for future
releases.  In particular, we'd like to update it to store more metadata
for samples.

Interoperability with mash
--------------------------

The default sketches computed by sourmash and mash are comparable, but
we are still [working on ways to convert the file
formats](https://github.com/marbl/Mash/issues/27).
