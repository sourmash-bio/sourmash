from __future__ import print_function
import sys
from setuptools import setup
from setuptools import Extension
import os

VERSION="0.2"

EXTRA_COMPILE_ARGS = ['-std=c++11', '-pedantic']
EXTRA_LINK_ARGS=[]

if os.environ.get('SOURMASH_COVERAGE'):
   print('Turning on coverage analysis.')
   EXTRA_COMPILE_ARGS.extend(['-g', '--coverage', '-lgcov'])
   EXTRA_LINK_ARGS=['--coverage', '-lgcov']
else:
    EXTRA_COMPILE_ARGS.append('-O3')

if sys.platform == 'darwin':
    # force 64bit only builds
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64', '-mmacosx-version-min=10.7',
                               '-stdlib=libc++'])

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
    "ext_modules": [Extension("_minhash",
                              sources=["_minhash.cc",
                                       "third-party/smhasher/MurmurHash3.cc"],
                              depends=["_minhash.hh", "kmer_min_hash.hh"],
                              language="c++",
                              extra_compile_args=EXTRA_COMPILE_ARGS,
                              extra_link_args=EXTRA_LINK_ARGS)],
    "scripts": ["sourmash"],
    "install_requires": ["khmer>=2.0", "PyYAML>=3.11"]
    }

setup(**SETUP_METADATA)

