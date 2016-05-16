all:
	python setup.py build_ext -i

clean:
	python setup.py clean --all

test: all
	py.test sourmash_lib.py sourmash_signature.py test__minhash.py
