# Developer quickstart with conda

The instructions below are for people interested in developing sourmash.
They should get you up and running with a sourmash development environment
in about 15 minutes.

## Install conda (or mamba)

(You don't need to do this if you already have conda and mamba installed!)

Follow the
[installation instructions](https://github.com/conda-forge/miniforge#install)
for installing `mambaforge` (a conda distribution that uses
[`mamba`](https://github.com/TheSnakePit/mamba) and the
[`conda-forge`](https://conda-forge.org/) channel by default).

## Create a conda environment with basic requirements:

First, install all of the necessary base packages:
```
mamba create -y -n sourmash-dev python=3.10 pip \
    cxx-compiler make rust
```

## Clone sourmash repo

Get a copy of the sourmash repository from GitHub:
```
git clone https://github.com/dib-lab/sourmash.git
cd sourmash
```

## Activate conda environment and install the development version of sourmash

Activate your new conda environment, created above:
```
conda activate sourmash-dev
```
and then install sourmash from within the git working directory:
```
python -m pip install -e ".[all]"
```

The `-e` option to `pip install` installs sourmash in "developer"
mode, such that the installed package is simply a link to the current
directory. Directing pip to install `".[all]"` results in pip
installing the package in the current directory (the `.` argument)
with all of the documentation and test packages (`[all]`) as specified
in `pyproject.toml`.

## Run tests

You can run all the Python tests like so:
```
python -m pytest 
```

## Develop!

You are now ready to develop on sourmash! Any code you edit or update will
be run by the tests and/or by executing `sourmash` at the command line.

If you change any sourmash Rust code, you will need run `make` before
running sourmash or the tests. This recompiles the Rust extension
code.

You don't need to run `make` if you're just changing Python code.

For example,
```
make
python -m pytest
```

Once you've run `make`, you can run `sourmash` in that conda
environment and it will run the version of sourmash in the git working
directory, and any changes to the Rust code will be included.

## Additional information

### Running subsets of tests

You can run specific subsets of tests using `pytest -k`. For example, this:
```
python -m pytest -k sbt
```
will run any test that has `sbt` in the name.

### Building documentation

You can build the docs like so:
```
cd doc
make
```
and then the built docs will be in `_build/html/`.

### Changing/creating new branches

As long as you're in the directory that you installed, you can change
branches and/or create new branches, and the `sourmash` command will
run the code on the branch.

So, for example,
```
git switch branch_from_github
```
will change to the `branch_from_github`, and then `git commit` and
`git push` will make commits and send those commits back to github.
