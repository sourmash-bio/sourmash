PYTHON ?= python

all:
	$(PYTHON) setup.py build_ext -i

.PHONY:

clean:
	$(PYTHON) setup.py clean --all
	cd doc && make clean

install:
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

test: all
	$(PYTHON) -m pytest

doc: .PHONY
	cd doc && make html

coverage: all
	$(PYTHON) setup.py clean --all
	SOURMASH_COVERAGE=1 $(PYTHON) setup.py build_ext -i
	$(PYTHON) -m pytest --cov=.

FORCE:
