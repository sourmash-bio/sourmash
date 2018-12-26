#!/usr/bin/env bash

set -e -x

PYTHON_VERSION=${1:-py36}

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    if [[ "${PYBIN}" == *"cp${PYTHON_VERSION:2:2}"* ]]; then
        "${PYBIN}/pip" install Cython
        "${PYBIN}/pip" wheel /io/ -w wheels/
        rm -rf /io/build /io/*.egg-info
    fi
done

# Bundle external shared libraries into the wheels
for whl in wheels/sourmash*.whl; do
    auditwheel show "$whl"
    auditwheel repair "$whl" -w /io/dist/
done
