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

offline:
	pip install -e . --no-index --find-links '.' --no-build-isolation

dist: FORCE
	$(PYTHON) -m build --sdist

wheel:
	$(PYTHON) -m maturin build -r

test: .PHONY
	tox -e py311
	cargo nextest run

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
	RUSTC_BOOTSTRAP=1 cbindgen -c cbindgen.toml . -o ../../$@ -v && \
	touch ../../$@

coverage: all
	tox -e coverage

benchmark:
	tox -e asv
	cargo bench

check:
	cargo build
	cargo nextest run
	cargo bench

last-tag:
	git fetch -p -q; git tag -l | sort -V | tail -1

wasm:
	wasm-pack build src/core -d ../../pkg -- --features 'niffler/wasm'

wasm-test:
	wasm-pack test --node src/core -- --features 'niffler/wasm'

wasi:
	cargo wasi build

FORCE:
