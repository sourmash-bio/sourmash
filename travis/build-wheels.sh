#!/bin/bash
set -e -x

rm -rf /opt/python/cp26-*
rm -rf /opt/python/cp34-*
rm -rf /opt/python/cp35-*

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/requirements-wheel.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
#for whl in wheelhouse/{sourmash,khmer}*.whl; do
for whl in wheelhouse/sourmash*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done
