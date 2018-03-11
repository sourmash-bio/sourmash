# Developer information

## Development environment


You can get the latest development master branch with:
```
git clone https://github.com/dib-lab/sourmash.git
```
To install all of the necessary dependencies, do:
```
pip install -r requirements.txt
```
Briefly, we use `py.test` for testing, and `coverage` for code
coverage analysis.

We suggest working on sourmash in a virtualenv; e.g. from within the
sourmash clone directory, you can do:
```
python -m virtualenv dev
. dev/bin/activate
pip install -e .
```
## Automated tests and code coverage calculation

We use [Travis][0] for continuous integration.

Code coverage calculation is enabled (on Linux only) by running
`make coverage`.  This recompiles the C++ extension without
optimization and with coverage configured.  See `setup.py` for
more information on this; the environment variable
`SOURMASH_COVERAGE` controls whether the C++ extension is
compiled with code coverage analysis enabled.

Code coverage can be viewed interactively at [codecov.io][1].

[0]:https://travis-ci.org/dib-lab/sourmash
[1]:https://codecov.io/gh/dib-lab/sourmash/
