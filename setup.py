import sys
from setuptools import setup
from setuptools import Extension

VERSION="0.2"

# Don't forget to update lib/Makefile with these flags!
EXTRA_COMPILE_ARGS = ['-g', '-O3', '-std=c++11', '-pedantic']

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
    "packages": ["sourmash_lib"],
    "ext_modules": [Extension("sourmash_lib._minhash",
                              sources=["sourmash_lib/_minhash.cc",
                                       "third-party/smhasher/MurmurHash3.cc"],
                              depends=["sourmash_lib/_minhash.hh",
                                       "sourmash_lib/kmer_min_hash.hh"],
                              language="c++",
                              extra_compile_args=EXTRA_COMPILE_ARGS)],
    "scripts": ["sourmash"],
    "install_requires": ["khmer>=2.0", "PyYAML>=3.11"]
    }

setup(**SETUP_METADATA)

