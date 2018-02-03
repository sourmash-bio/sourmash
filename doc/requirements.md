# Computational requirements


sourmash has no particular memory requirements; it will need to hold
the largest single sequence you have in memory, but the individual
signatures are quite small and we do no special buffer allocation.

sourmash's intensive computation is almost entirely computing k-mers
and hashes of k-mers.  It will create a sketch of several megabases
in a second or so on a rather slow 2016 Mac laptop.

MinHash sketches and signatures are quite small on disk.

sourmash should run with little modification on Linux and Mac OS X,
under Python 2.7.11 and Python 3.5.  Please see [the development repository README][0]
for
information on source code, tests, and continuous integration.
[0]:https://github.com/dib-lab/sourmash/blob/master/README.md
