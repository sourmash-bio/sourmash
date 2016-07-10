from __future__ import print_function
import sys
from setuptools import setup
from setuptools import Extension
import os

VERSION="0.9.3"

EXTRA_COMPILE_ARGS = ['-std=c++11', '-pedantic']
EXTRA_LINK_ARGS=[]

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

if "-rc" in VERSION:
    CLASSIFIERS.append("Development Status :: 4 - Beta")
else:
    CLASSIFIERS.append("Development Status :: 5 - Production/Stable")

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
                              include_dirs=["./sourmash_lib",
                                            "./third-party/smhasher/"],
                              language="c++",
                              extra_compile_args=EXTRA_COMPILE_ARGS,
                              extra_link_args=EXTRA_LINK_ARGS)],
    "scripts": ["sourmash"],
    "install_requires": ["screed>=0.9", "PyYAML>=3.11"],
    "extras_require": {
        'test' : ['pytest', 'pytest-cov', 'numpy', 'matplotlib', 'scipy'],
        'fig' : ['numpy', 'matplotlib', 'scipy'],
        'demo' : ['numpy', 'matplotlib', 'scipy', 'jupyter',
                  'jupyter_client', 'ipython'],
        'doc' : ['sphinx'],
        },
    "include_package_data": True,
    "classifiers": CLASSIFIERS
    }

setup(**SETUP_METADATA)

