PYTHON ?= python

all: build

.PHONY:

build: .PHONY
	$(PYTHON) -m pip install -e .

clean:
	$(PYTHON) -m pip uninstall -y sourmash
	rm -rf src/sourmash/_lowlevel
	cd doc && make clean

install: build

dist: FORCE
	$(PYTHON) -m build --sdist

test: .PHONY
	tox -e py39
	cargo test

doc: .PHONY
	tox -e docs

include/sourmash.h: src/core/src/lib.rs \
                    src/core/src/ffi/mod.rs \
                    src/core/src/ffi/hyperloglog.rs \
                    src/core/src/ffi/minhash.rs \
                    src/core/src/ffi/signature.rs \
                    src/core/src/ffi/nodegraph.rs \
                    src/core/src/ffi/index/mod.rs \
                    src/core/src/ffi/index/revindex.rs \
                    src/core/src/ffi/storage.rs \
                    src/core/src/errors.rs \
                    src/core/cbindgen.toml
	cd src/core && \
	RUSTC_BOOTSTRAP=1 cbindgen -c cbindgen.toml . -o ../../$@

coverage: all
	tox -e coverage

benchmark:
	tox -e asv
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
