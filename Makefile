PYTHON ?= python

all:
	$(PYTHON) setup.py build_ext -i

clean:
	$(PYTHON) setup.py clean --all

test: all
	$(PYTHON) -m pytest sourmash_lib/*.py
