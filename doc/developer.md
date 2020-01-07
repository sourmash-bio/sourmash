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

## Code organization

There are three main components in the sourmash repo:
- Python module (in `sourmash/`)
- The command-line interface (in `sourmash/cli`)
- The Rust core library  (in `src/core`)

`setup.py` has all the configuration to prepare a Python package containing these three components.
First it compiles the Rust core component into a shared library,
which is wrapped by CFFI and exposed to the Python module.

A short description of the high-level files and dirs in the sourmash repo:
```
.
├── benchmarks/         | Benchmarks for the Python module
├── binder/             | mybinder.org configuration
├── data/               | data used for demos
├── doc/                | the documentation rendered in sourmash.bio
├── include/            | C/C++ header files for using core library
├── sourmash/           | The Python module and CLI code
├── sourmash_lib/       | DEPRECATED: previous name of the Python module
├── src/                |
│   └── core            | Code for the core library (Rust)
├── tests/              | Tests for the Python module and CLI
├── utils/              |
├── asv.conf.json       | benchmarking config file (for ASV)
├── Cargo.toml          | Rust definition for a workspace
├── CITATION.cff        | Citation info
├── codemeta.json       | Metadata for software discovery
├── CODE_OF_CONDUCT.rst | Code of conduct
├── CONTRIBUTING.md     | Instruction for contributing to development
├── LICENSE             | License for the repo
├── Makefile            | Entry point for most development tasks
├── MANIFEST.in         | Describes what files to add to the Python package
├── matplotlibrc        | Configuration for matplotlib
├── netlify.toml        | Configuration for netlify (build docs for preview)
├── paper.bib           | References in the JOSS paper
├── paper.md            | JOSS paper content
├── pytest.ini          | pytest configuration
├── README.md           | Info to get started
├── requirements.txt    | Python dependencies for development
├── setup.py            | Python package definition
└── tox.ini             | Configuration for test automation
```

### The Python module

```
sourmash
├── cli/
├── command_compute.py
├── commands.py
├── compare.py
├── _compat.py
├── exceptions.py
├── fig.py
├── index.py
├── __init__.py
├── lca
│   ├── command_classify.py
│   ├── command_compare_csv.py
│   ├── command_gather.py
│   ├── command_index.py
│   ├── command_rankinfo.py
│   ├── command_summarize.py
│   ├── __init__.py
│   ├── lca_utils.py
│   └── __main__.py
├── logging.py
├── __main__.py
├── _minhash.py
├── np_utils.py
├── sbtmh.py
├── sbt.py
├── sbt_storage.py
├── search.py
├── sig
│   ├── __init__.py
│   └── __main__.py
├── signature_json.py
├── signature.py
├── sourmash_args.py
├── utils.py
└── version.py
```

### The Python CLI

```
sourmash
└── cli/
```

### The Rust core library

This is completely defined in `src/core` to avoid mixing with the code of other components
(and trying to make it easier to reason about changes).
If you're only working on the core,
you don't need to change any files outside this directory.

This is also published to [crates.io] (the Rust package repository) and [NPM],
after it is compiled to Webassembly.
The GitHub Actions workflow publishes new versions automatically to these
repositories.

[crates.io]: https://crates.io/crates/sourmash
[NPM]: https://www.npmjs.com/package/sourmash

```
src/core
├── benches              |
│   └── index.rs
├── Cargo.toml
├── cbindgen.toml
├── examples/
├── README.md
├── src
│   ├── cmd.rs
│   ├── errors.rs
│   ├── ffi
│   │   ├── minhash.rs
│   │   ├── mod.rs
│   │   ├── nodegraph.rs |
│   │   ├── signature.rs
│   │   └── utils.rs
│   ├── from.rs
│   ├── index
│   │   ├── bigsi.rs
│   │   ├── linear.rs
│   │   ├── mod.rs
│   │   ├── sbt
│   │   │   ├── mhbt.rs
│   │   │   ├── mhmt.rs
│   │   │   └── mod.rs
│   │   ├── search.rs
│   │   └── storage.rs
│   ├── lib.rs
│   ├── signature.rs
│   ├── sketch
│   │   ├── minhash.rs
│   │   ├── mod.rs
│   │   ├── nodegraph.rs
│   │   └── ukhs.rs
│   └── wasm.rs
└── tests
    ├── finch.rs
    ├── minhash.rs
    ├── signature.rs
    └── test.rs
```

## Versioning

We use [`setuptools_scm`] to generate versions based on git tags.
Versions are tagged in a `vMAJOR.MINOR.PATH` format,
following the [Semantic Versioning] convention.
From their definition:

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> MAJOR version when you make incompatible API changes,
> MINOR version when you add functionality in a backwards compatible manner, and
> PATCH version when you make backwards compatible bug fixes.

[`setuptools_scm`]: https://github.com/pypa/setuptools_scm
[Semantic Versioning]: https://semver.org/

For the Rust core library we use `rMAJOR.MINOR.PATH`
(note it starts with `r`, and not `v`).
The Rust version is not automated,
and must be bumped in `src/core/Cargo.toml`.
