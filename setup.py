from __future__ import print_function
import sys
from setuptools import setup, find_packages
from setuptools import Extension
import os

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
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

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

with open('README.md', 'r') as readme:
    LONG_DESCRIPTION = readme.read()

SETUP_METADATA = \
               {
    "name": "sourmash",
    "description": "tools for comparing DNA sequences with MinHash sketches",
    "long_description": LONG_DESCRIPTION,
    "long_description_content_type": "text/markdown",
    "url": "https://github.com/dib-lab/sourmash",
    "author": "C. Titus Brown",
    "author_email": "titus@idyll.org",
    "license": "BSD 3-clause",
    "packages": find_packages(exclude=["tests", "benchmarks"]),
    "entry_points": {'console_scripts': [
        'sourmash = sourmash.__main__:main'
        ]
    },
    "ext_modules": [Extension("sourmash._minhash",
                               sources=["sourmash/_minhash.pyx",
                                        "third-party/smhasher/MurmurHash3.cc"],
                               depends=["sourmash/kmer_min_hash.hh"],
                               include_dirs=["./sourmash",
                                             "./third-party/smhasher/"],
                               language="c++",
                               extra_compile_args=EXTRA_COMPILE_ARGS,
                               extra_link_args=EXTRA_LINK_ARGS)],
    "install_requires": ["screed>=0.9", "ijson", "khmer>=2.1"],
    "setup_requires": ['Cython>=0.25.2', "setuptools>=38.6.0",
                       'setuptools_scm', 'setuptools_scm_git_archive'],
    "use_scm_version": {"write_to": "sourmash/version.py"},
    "extras_require": {
        'test' : ['pytest', 'pytest-cov', 'numpy', 'matplotlib', 'scipy','recommonmark'],
        'demo' : ['jupyter', 'jupyter_client', 'ipython'],
        'doc' : ['sphinx'],
        '10x': ['pathos', 'bamnostic>=0.9.2'],
        },
    "include_package_data": True,
    "package_data": {
        "sourmash": ['*.pxd']
    },
    "classifiers": CLASSIFIERS
    }

setup(**SETUP_METADATA)

