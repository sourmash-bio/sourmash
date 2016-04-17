all:
	python setup.py build_ext -i

clean:
	python setup.py clean

test:
	py.test sourmash_lib.py sourmash_signature.py
