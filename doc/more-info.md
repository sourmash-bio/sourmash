# Additional information on sourmash


## Other MinHash implementations for DNA


In addition to [mash][0], also see:

* [RKMH][1]: Read Classification by Kmers
* [mashtree][2]:For building trees using Mash
* [Finch][3]:"Fast sketches,
  count histograms, better filtering."

If you are interested in exactly how these MinHash approaches
calculate the hashes of DNA sequences, please see some simple Python
code in sourmash, [utils/compute-dna-mh-another-way.py][4]


## Blog posts


We have a number of blog posts on sourmash and MinHash more generally:

* [Applying MinHash to cluster RNAseq samples][5]


* [MinHash signatures as ways to find samples, and collaborators?][6]


* [Efficiently searching MinHash Sketch collections][7]:
   indexing and
  search 42,000 bacterial genomes with Sequence Bloom Trees.

* [Quickly searching all the microbial genomes, mark 2 - now with archaea, phage, fungi, and protists!][8]:
 indexing
  and searching 50,000 microbial genomes, round 2.

* [What metadata should we put in MinHash Sketch signatures?][9]:
  crowdsourcing ideas for what metadata belongs in a signature file.

* [Minhashing all the things (part 1): microbial genomes][10]:
   on
  approaches to computing MinHashes for large collections of public data.

## JSON format for the signature


The JSON format is not necessarily final; this is a TODO item for future
releases.  In particular, we'd like to update it to store more metadata
for samples.

## Interoperability with mash


The default sketches computed by sourmash and mash are comparable, but
we are still [working on ways to convert the file formats][11]

## Developing sourmash


Please see:

.. toctree::
   :maxdepth: 2

   developer
   release

[0]:https://github.com/marbl/Mash
[1]:https://github.com/edawson/rkmh
[2]:https://github.com/lskatz/mashtree/blob/master/README.md
[3]:https://github.com/onecodex/finch-rs
[4]:https://github.com/dib-lab/sourmash/blob/master/utils/compute-dna-mh-another-way.py
[5]:http://ivory.idyll.org/blog/2016-sourmash.html
[6]:http://ivory.idyll.org/blog/2016-sourmash-signatures.html
[7]:http://ivory.idyll.org/blog/2016-sourmash-sbt.html
[8]:http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html
[9]:http://ivory.idyll.org/blog/2016-sourmash-signatures-metadata.html
[10]:http://blog.luizirber.org/2016/12/28/soursigs-arch-1/
[11]:https://github.com/marbl/Mash/issues/27
