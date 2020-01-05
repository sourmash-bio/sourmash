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

include/sourmash.h: src/core/src/lib.rs \
                    src/core/src/ffi/minhash.rs \
                    src/core/src/ffi/signature.rs \
                    src/core/src/ffi/nodegraph.rs \
                    src/core/src/errors.rs
	rustup override set nightly
	cd src/core && \
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
	wasm-pack build src/core -d ../../pkg

wasi:
	cargo wasi build

FORCE:
