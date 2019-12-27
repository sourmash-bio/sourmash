PYTHON ?= python

all: build

.PHONY:

build: .PHONY
	$(PYTHON) setup.py build_ext -i
	cargo build

clean:
	$(PYTHON) setup.py clean --all
	rm -f sourmash/*.so
	cd doc && make clean

install: all
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

test: all
	pip install -e '.[test]'
	$(PYTHON) -m pytest
	cargo test

doc: .PHONY
	cd doc && make html

include/sourmash.h: crates/sourmash/src/lib.rs \
                    crates/sourmash/src/ffi/minhash.rs \
                    crates/sourmash/src/ffi/signature.rs \
                    crates/sourmash/src/ffi/nodegraph.rs \
                    crates/sourmash/src/errors.rs
	rustup override set nightly
	cd crates/sourmash && \
	RUST_BACKTRACE=1 cbindgen -c cbindgen.toml -o ../../$@
	rustup override set stable

coverage: all
	$(PYTHON) setup.py build_ext -i
	$(PYTHON) -m pytest --cov=. --cov-report term-missing

benchmark:
	asv continuous master `git rev-parse HEAD`
	cargo bench

check:
	cargo build
	cargo test
	cargo bench

last-tag:
	git fetch -p -q; git tag -l | sort -V | tail -1

wasm:
	wasm-pack build

wasi:
	cargo wasi build

FORCE:
