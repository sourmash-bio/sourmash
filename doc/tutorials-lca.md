# Using sourmash LCA to do taxonomic classification

## Installing sourmash

To install sourmash, run:

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

pip install -U https://github.com/dib-lab/sourmash/archive/add/lca.zip

```

## Building an LCA database

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

## Using the LCA database to classify signatures

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

[Return to index][3]

[0]:http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html
[1]:databases.html
[2]:https://www.ncbi.nlm.nih.gov/pubmed/233877
[3]:index.html
