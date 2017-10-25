# Welcome to sourmash!


sourmash is a command-line tool and Python library for computing
[MinHash sketches][0] from DNA
sequences, comparing them to each other, and plotting the results.
This allows you to estimate sequence similarity between even very
large data sets quickly and accurately.

sourmash can also be used to quickly search large databases of genomes
for matches to query genomes and metagenomes; see [our list of available databases][9]

Please see the [mash][1] software and the [mash paper (Ondov et al., 2016)][2] for background
information on how and why MinHash sketches work.

-----

To use sourmash, you must be comfortable with the UNIX command line;
programmers may find the Python library and API useful as well.

In brief,

* `sourmash` provides command line utilities for creating, comparing,
  and searching MinHash sketches, as well as plotting and clustering
  sketches by distance (see [the command-line][10])

* `sourmash` supports saving, loading, and communication of MinHash
  sketches via [JSON][3], a ~human-readable &
  editable format.

* `sourmash` also has a simple Python API for interacting with sketches,
  including support for online updating and querying of sketches
  (see [the API docs][11])

* `sourmash` isn't terribly slow, and relies on an underlying CPython
  module.

* `sourmash` is developed [on GitHub][4]
  and is freely and openly
  available under the BSD 3-clause license.  Please see [the README][5] for
  more information on development, support, and contributing.

You can take a look at sourmash analyses on real data [in a saved Jupyter notebook][6],
and experiment with it yourself [interactively with a binder][7] at [mybinder.org][8]


[0]:https://en.wikipedia.org/wiki/MinHash
[1]:http://mash.readthedocs.io/en/latest/
[2]:http://biorxiv.org/content/early/2015/10/26/029827
[3]:http://www.json.org/
[4]:https://github.com/dib-lab/sourmash
[5]:https://github.com/dib-lab/sourmash/blob/master/README.md
[6]:https://github.com/dib-lab/sourmash/blob/master/demo/00-demo.ipynb
[7]:http://mybinder.org/repo/dib-lab/sourmash
[8]:http://mybinder.org
[9]:[databases.html]
[10]:[the-command-line.html]
[11]:[api.html]
