# Developer information

## Development environment


You can get the latest development master branch with:
```
git clone https://github.com/dib-lab/sourmash.git
```
sourmash runs under both Python 2.7.x and Python 3.5+.  The base
requirements are screed and ijson, together with a Rust environment (for the
extension code). We suggest using `rustup` to install the Rust environment:

    curl https://sh.rustup.rs -sSf | sh

To install all of the necessary Python dependencies, do:
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

You can run tests by invoking `make test` or `python -m pytest` in the sourmash
directory.

## Automated tests and code coverage calculation

We use [Travis][0] and [GitHub Actions][2] for continuous integration.

Code coverage can be viewed interactively at [codecov.io][1].

[0]: https://travis-ci.com/dib-lab/sourmash
[1]: https://codecov.io/gh/dib-lab/sourmash/
[2]: https://github.com/dib-lab/sourmash/actions
