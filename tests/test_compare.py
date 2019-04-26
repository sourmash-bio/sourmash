import glob
import os

import pytest

import sourmash
from sourmash.compare import compare_all_pairs, _compare_serial, memmap_siglist
from . import sourmash_tst_utils as utils


@pytest.fixture()
def siglist():
    demo_path = utils.get_test_data("demo")
    filenames = glob.glob(os.path.join(demo_path, "*.sig"))
    sigs = []
    for filename in filenames:
        sigs.extend(sourmash.load_signatures(filename))
    return sigs


@pytest.fixture()
def ignore_abundance(track_abundance):
    return not track_abundance


def test__compare_serial(siglist, ignore_abundance):
    similarities = _compare_serial(siglist, ignore_abundance)
    true_similarities = '''[[1.    0.    0.074 0.    0.    0.078 0.086]
 [0.    1.    0.    0.382 0.364 0.    0.   ]
 [0.074 0.    1.    0.    0.    0.072 0.078]
 [0.    0.382 0.    1.    0.386 0.    0.   ]
 [0.    0.364 0.    0.386 1.    0.    0.   ]
 [0.078 0.    0.072 0.    0.    1.    0.356]
 [0.086 0.    0.078 0.    0.    0.356 1.   ]]'''
    names = ['SRR2255622_1.fastq.gz', 'SRR453566_1.fastq.gz',
             'SRR2241509_1.fastq.gz', 'SRR453569_1.fastq.gz',
             'SRR453570_1.fastq.gz', 'SRR2060939_2.fastq.gz',
             'SRR2060939_1.fastq.gz']
    assert names == [x.name() for x in siglist]
    assert str(similarities) == true_similarities


def test_memmap_siglist(siglist):
    memmapped = memmap_siglist(siglist)

    # Assert that the data didn't change as a result of memory-mapping
    assert memmapped == siglist


def test_compare_all_pairs(siglist, ignore_abundance):
    similarities = compare_all_pairs(siglist, ignore_abundance, n_jobs=2)
    assert str(similarities) == str(_compare_serial(siglist, ignore_abundance))
