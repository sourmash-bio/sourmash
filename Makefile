all:
	python setup.py build_ext -i

clean:
	python setup.py clean
	-rm -f *.so

test: all
	py.test sourmash_lib.py sourmash_signature.py test__minhash.py
