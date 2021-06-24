# Additional information on sourmash

## Computational requirements

Read more about the [computational requirements, here.](requirements.md)

## Prepared search database

We offer a number of [prepared search databases.](databases.md)

## Other MinHash implementations for DNA

In addition to [mash][0], also see:

* [RKMH][1]: Read Classification by Kmers
* [mashtree][2]: For building trees using Mash
* [Finch][3]: "Fast sketches,
  count histograms, better filtering."
* [BBMap and SendSketch][13]: part of Brian Bushnell's tool collection.
* [PATRIC](https://patricbrc.org/) uses MinHash for genome search.

If you are interested in exactly how these MinHash approaches
calculate the hashes of DNA sequences, please see some simple Python
code in sourmash, [utils/compute-dna-mh-another-way.py][4]

## Papers and references

[On the resemblance and containment of documents][20],  Broder, 1997. The original MinHash paper!

[Mash: fast genome and metagenome distance estimation using MinHash.][21], Ondov et al. 2016.

[sourmash: a library for MinHash sketching of DNA.][22], Brown and Irber, 2017.

[Improving MinHash via the Containment Index with Applications to Metagenomic Analysis][23], Koslicki and Zabeti, 2017.

[Ultra-fast search of all deposited bacterial and viral genomic data](http://dx.doi.org/10.1038/s41587-018-0010-1), Bradley et al., 2019.
 
[Streaming histogram sketching for rapid microbiome analytics](https://www.biorxiv.org/content/10.1101/408070v1), Rowe et al., 2019.

[Dashing: Fast and Accurate Genomic Distances with HyperLogLog](https://www.biorxiv.org/content/10.1101/501726v2), Baker and Langmead, 2019.

## Presentations and posters

[Taxonomic classification of microbial metagenomes using MinHash signatures][12], Brooks et al., 2017. Presented at ASM.

## Blog posts


 We (and others) have a number of blog posts on sourmash and MinHash
 more generally:

 * [Some background on k-mers, and their use in taxonomy][15]

 * From the Phillippy lab: [mash screen: what's in my sequencing run?][14]

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
* [Comparing genome sets extracted from metagenomes][16]

* [Taxonomic examinations of genome bins from Tara Oceans][17]

* [Classifying genome bins using a custom reference database, part I][18]

* [Classifying genome bins using a custom reference database, part II][19]

## JSON format for the signature


The JSON format is not necessarily final; this is a TODO item for future
releases.  In particular, we'd like to update it to store more metadata
for samples.

## Interoperability with mash


The hashing functions used by sourmash and mash are the same, but we
are still [working on ways to convert the file formats][11]. Please
keep an eye on `sourmash signature import` and `sourmash signature
export`!

## Developing sourmash


Please see:
 * [Developer information][24]
 * [Release workflow][25]

## Known issues

For at least some versions of matplotlib, users may encounter an
error "Failed to connect to server socket:" or "RuntimeError:
Invalid DISPLAY variable".  This is because by default matplotlib
tries to connect to X11 to use the Tkinter backend.

The solution is to force the use of the 'Agg' backend in matplotlib;
see [this stackoverflow answer](https://stackoverflow.com/a/34294056)
or [this sourmash issue comment](https://github.com/sourmash-bio/sourmash/issues/254#issuecomment-304274590).

Newer versions of matplotlib do not seem to have this problem.

[0]:https://github.com/marbl/Mash
[1]:https://github.com/edawson/rkmh
[2]:https://github.com/lskatz/mashtree/blob/master/README.md
[3]:https://github.com/onecodex/finch-rs
[4]:https://github.com/sourmash-bio/sourmash/blob/latest/utils/compute-dna-mh-another-way.py
[5]:http://ivory.idyll.org/blog/2016-sourmash.html
[6]:http://ivory.idyll.org/blog/2016-sourmash-signatures.html
[7]:http://ivory.idyll.org/blog/2016-sourmash-sbt.html
[8]:http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html
[9]:http://ivory.idyll.org/blog/2016-sourmash-signatures-metadata.html
[10]:http://blog.luizirber.org/2016/12/28/soursigs-arch-1/
[11]:https://github.com/marbl/Mash/issues/27
[12]:https://osf.io/mu4gk/
[13]:http://seqanswers.com/forums/showthread.php?t=74019
[14]:https://genomeinformatics.github.io/mash-screen/
[15]:http://ivory.idyll.org/blog/2017-something-about-kmers.html
[16]:http://ivory.idyll.org/blog/2017-comparing-genomes-from-metagenomes.html
[17]:http://ivory.idyll.org/blog/2017-taxonomy-of-tara-ocean-genomes.html
[18]:http://ivory.idyll.org/blog/2017-classify-genome-bins-with-custom-db-part-1.html
[19]:http://ivory.idyll.org/blog/2017-classify-genome-bins-with-custom-db-part-2.html
[20]:http://ieeexplore.ieee.org/document/666900/?reload=true
[21]:https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x
[22]:http://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9
[23]:https://www.biorxiv.org/content/early/2017/09/04/184150
[24]:developer.md
[25]:release.md
