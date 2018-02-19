# Using sourmash LCA to do taxonomic classification

The `sourmash lca` sub-commands do k-mer classification using an
"lowest common ancestor" approach.  See "Some discussion" below for
links and details.

## Quickstart

(These `sourmash lca classify` and `sourmash lca summarize` steps require
about 4 GB of RAM when using the genbank database, as below.)

First, install sourmash 2.0a4 or later.
```
pip install -U https://github.com/dib-lab/sourmash/archive/master.zip
```

Next, download a genbank LCA database for k=31:

```
curl -L -o genbank-k31.lca.json.gz https://osf.io/zskb9/download
```

Download a random genome from genbank:
```
curl -L -o some-genome.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/178/875/GCF_000178875.2_ASM17887v2/GCF_000178875.2_ASM17887v2_genomic.fna.gz
```

Compute a signature for this genome:
```
sourmash compute -k 31 --scaled=1000 some-genome.fa.gz
```

Now, classify the signature with sourmash `lca classify`,
```
sourmash lca classify --db genbank-k31.lca.json.gz \
    --query some-genome.fa.gz.sig
```
and this will give you a taxonomic identification of your genome bin,
classify using all of the genbank microbial genomes.

You can also summarize the taxonomic distribution of the content with
`lca summarize`:
```
sourmash lca summarize --db genbank-k31.lca.json.gz \
    --query some-genome.fa.gz.sig
```

To apply this to your own genome(s), replace `some-genome.fa.gz` above
with your own filename(s).

You can also specify multiple databases and multiple query signatures
on the command line; separate them with `--db` or `--query`.

## A longer tutorial

Install sourmash as above; see Appendix (below) for dependencies.

Let's start by building your own LCA database, using your own taxonomy.

### Building your own LCA database

(This is an abbreviated version of [this blog post](http://ivory.idyll.org/blog/2017-classify-genome-bins-with-custom-db-try-again.html), updated to use the `sourmash lca` commands.)

Download some pre-computed signatures:

```
curl -L https://osf.io/bw8d7/download?version=1 -o delmont-subsample-sigs.tar.gz
tar xzf delmont-subsample-sigs.tar.gz
```

Next, grab the associated taxonomy spreadsheet

```
curl -O -L https://github.com/ctb/2017-sourmash-lca/raw/master/tara-delmont-SuppTable3.csv
```

Build a sourmash LCA database named `delmont.lca.json`:

```
sourmash lca index tara-delmont-SuppTable3.csv delmont.lca.json delmont-subsample-sigs/*.sig
```

### Using the LCA database to classify signatures

We can now use `delmont.lca.json` to classify signatures with k-mers
according to the database we just created.  (Note, the database is
completely self-contained at this point.)

Let's classify a single signature:
```
sourmash lca classify --db delmont.lca.json \
    --query delmont-subsample-sigs/TARA_RED_MAG_00003.fa.gz.sig
```

and you should see:

```
loaded 1 databases for LCA use.
ksize=31 scaled=10000
outputting classifications to stdout
ID,status,superkingdom,phylum,class,order,family,genus,species
TARA_RED_MAG_00003,found,Bacteria,Proteobacteria,Gammaproteobacteria,,,,
classified 1 signatures total
```

You can classify a bunch of signatures and also specify an output
location for the CSV:

```
sourmash lca classify --db delmont.lca.json \
    --query delmont-subsample-sigs/*.sig \
    -o out.csv
```

The `lca classify` command supports multiple databases as well as
multiple queries; e.g. `sourmash lca classify --db delmont.lca.json
other.lca.json` will classify based on the combination of taxonomies
in the two databases.

## Some discussion

Sourmash LCA is using k-mers to do taxonomic classification, using the
"lowest common ancestor" approach (pioneered by
[Kraken](http://ccb.jhu.edu/software/kraken/MANUAL.html), and described
[here](http://ivory.idyll.org/blog/2017-something-about-kmers.html)),
to identify each k-mer.  From this it can either find a consensus
taxonomy between all the k-mers (`sourmash classify`) or it can summarize
the mixture of k-mers present in one or more signatures (`sourmash summarize`).

The `sourmash lca index` command can be used to prepare custom taxonomy
databases; sourmash will happily ingest any taxonomy, whether or not
it matches NCBI. See
[the spreadsheet from Delmont et al., 2017](https://github.com/ctb/2017-sourmash-lca/blob/master/tara-delmont-SuppTable3.csv)
for an example format.

## Appendix: Installing sourmash from scratch

To install sourmash on an Ubuntu or Debian system, run:

```
sudo apt-get -y update && \
sudo apt-get install -y python3.5-dev python3.5-venv make \
    libc6-dev g++ zlib1g-dev
```

this installs Python 3.5.

Now, create a local software install and populate it with Jupyter and
other dependencies:

```
python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip
pip install -U Cython
pip install -U jupyter jupyter_client ipython pandas matplotlib scipy scikit-learn khmer
```

Last but not least, install `sourmash` from the LCA branch:

```
pip install -U https://github.com/dib-lab/sourmash/archive/add/lca.zip
```


[Return to index][3]

[0]:http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html
[1]:databases.html
[2]:https://www.ncbi.nlm.nih.gov/pubmed/233877
[3]:index.html
