PYTHON ?= python

all:
	$(PYTHON) setup.py build_ext -i

clean:
	$(PYTHON) setup.py clean --all

test: all
	$(PYTHON) -m pytest sourmash_lib/*.py

coverage: all
	$(PYTHON) setup.py clean --all
	SOURMASH_COVERAGE=1 $(PYTHON) setup.py build_ext -i
	$(PYTHON) -m pytest --cov=. sourmash_lib/*.py
