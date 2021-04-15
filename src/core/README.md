<meta charset="utf-8"/>

# sourmash

🦀
[![](http://meritbadge.herokuapp.com/sourmash)](https://crates.io/crates/sourmash)
[![Rust API Documentation on docs.rs](https://docs.rs/sourmash/badge.svg)](https://docs.rs/sourmash)
[![build-status]][github-actions]
[![codecov](https://codecov.io/gh/dib-lab/sourmash/branch/latest/graph/badge.svg)](https://codecov.io/gh/dib-lab/sourmash)
<a href="https://github.com/dib-lab/sourmash/blob/latest/LICENSE"><img alt="License: 3-Clause BSD" src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg"></a>

[build-status]: https://github.com/dib-lab/sourmash/workflows/Rust%20checks/badge.svg
[github-actions]: https://github.com/dib-lab/sourmash/actions?query=workflow%3A%22Rust+checks%22

----

Compute MinHash signatures for nucleotide (DNA/RNA) and protein sequences.

This is the core library used by sourmash. It exposes a C API that can be
called from FFI in other languages, and it is how we use it in Python for
building the sourmash application (CLI and Python API).

----

sourmash is a product of the
[Lab for Data-Intensive Biology](http://ivory.idyll.org/lab/) at the
[UC Davis School of Veterinary Medicine](http://www.vetmed.ucdavis.edu).

## Support

Please ask questions and files issues
[on Github](https://github.com/dib-lab/sourmash/issues).

## Development

Development happens on github at
[dib-lab/sourmash](https://github.com/dib-lab/sourmash).

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.42.0.
