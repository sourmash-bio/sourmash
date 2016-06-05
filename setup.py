from __future__ import print_function
import sys
from setuptools import setup
from setuptools import Extension
import os

VERSION="0.3"

EXTRA_COMPILE_ARGS = ['-std=c++11', '-pedantic']
EXTRA_LINK_ARGS=[]

if sys.platform == 'darwin':              # Mac OS X?
    # force 64bit only builds
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64', '-mmacosx-version-min=10.7',
                               '-stdlib=libc++'])

else:                                     # ...likely Linux
   if os.environ.get('SOURMASH_COVERAGE'):
      print('Turning on coverage analysis.')
      EXTRA_COMPILE_ARGS.extend(['-g', '--coverage', '-lgcov'])
      EXTRA_LINK_ARGS.extend(['--coverage', '-lgcov'])
   else:
      EXTRA_COMPILE_ARGS.append('-O3')

SETUP_METADATA = \
               {
    "name": "sourmash",
    "version": VERSION,
    "description": "tools for comparing DNA sequences with MinHash sketches",
    "url": "https://github.com/dib-lab/sourmash",
    "author": "C. Titus Brown",
    "author_email": "titus@idyll.org",
    "license": "BSD 3-clause",
    "packages": ["sourmash_lib"],
    "ext_modules": [Extension("sourmash_lib._minhash",
                              sources=["sourmash_lib/_minhash.cc",
                                       "third-party/smhasher/MurmurHash3.cc"],
                              depends=["sourmash_lib/_minhash.hh",
                                       "sourmash_lib/kmer_min_hash.hh"],
                              language="c++",
                              extra_compile_args=EXTRA_COMPILE_ARGS,
                              extra_link_args=EXTRA_LINK_ARGS)],
    "scripts": ["sourmash"],
    "install_requires": ["screed>=0.9", "PyYAML>=3.11"]
    }

setup(**SETUP_METADATA)

