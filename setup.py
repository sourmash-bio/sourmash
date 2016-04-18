import sys
from setuptools import setup
from setuptools import Extension

VERSION="0.2"

SETUP_METADATA = \
               {
    "name": "sourmash",
    "version": VERSION,
    "description": "tools for comparing DNA sequences with MinHash sketches",
    "url": "https://github.com/dib-lab/sourmash",
    "author": "C. Titus Brown",
    "author_email": "titus@idyll.org",
    "license": "BSD 3-clause",
    "py_modules": ["sourmash_lib","sourmash_signature"],
    "scripts": ["sourmash"],
    "install_requires": ["khmer>=2.0", "PyYAML>=3.11"]
    }

setup(**SETUP_METADATA)

