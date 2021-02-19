import gzip
from tempfile import NamedTemporaryFile

from screed.fasta import fasta_iter
import pytest

from sourmash.hll import HLL

import sourmash_tst_utils as utils

K = 21  # size of kmer
ERR_RATE = 0.01
N_UNIQUE = 3356
TRANSLATE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}


def test_hll_add_python():
    # test python code to count unique kmers using HyperLogLog.
    # use the lower level add() method, which accepts anything,
    # and compare to an exact count using collections.Counter

    filename = utils.get_test_data('ecoli.genes.fna')
    hll = HLL(ERR_RATE, K)
    counter = set()

    for n, record in enumerate(fasta_iter(open(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            rc = "".join(TRANSLATE[c] for c in kmer[::-1])

            hll.add(kmer)

            if rc in counter:
                kmer = rc
            counter.update([kmer])

    n_unique = len(counter)

    assert n_unique == N_UNIQUE
    assert abs(1 - float(hll.cardinality()) / N_UNIQUE) < ERR_RATE


def test_hll_consume_string():
    # test rust code to count unique kmers using HyperLogLog,
    # using screed to feed each read to the counter.

    filename = utils.get_test_data('ecoli.genes.fna')
    hll = HLL(ERR_RATE, K)
    n_consumed = n = 0
    for n, record in enumerate(fasta_iter(open(filename)), 1):
        hll.add_sequence(record['sequence'])

    assert abs(1 - float(len(hll)) / N_UNIQUE) < ERR_RATE


def test_hll_similarity_containment():
    N_UNIQUE_H1 = 500741
    N_UNIQUE_H2 = 995845
    N_UNIQUE_U = 995845

    SIMILARITY = 0.502783
    CONTAINMENT_H1 = 1.
    CONTAINMENT_H2 = 0.502783

    INTERSECTION = 500838

    hll1 = HLL(ERR_RATE, K)
    hll2 = HLL(ERR_RATE, K)
    hllu = HLL(ERR_RATE, K)

    filename = utils.get_test_data('genome-s10.fa.gz')
    for n, record in enumerate(fasta_iter(gzip.GzipFile(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            hll1.add(kmer)
            hllu.add(kmer)

    filename = utils.get_test_data('genome-s10+s11.fa.gz')
    for n, record in enumerate(fasta_iter(gzip.GzipFile(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            hll2.add(kmer)
            hllu.add(kmer)

    assert abs(1 - float(hll1.cardinality()) / N_UNIQUE_H1) < ERR_RATE
    assert abs(1 - float(hll2.cardinality()) / N_UNIQUE_H2) < ERR_RATE

    assert abs(1 - float(hll1.similarity(hll2)) / SIMILARITY) < ERR_RATE

    assert abs(1 - float(hll1.containment(hll2)) / CONTAINMENT_H1) < ERR_RATE
    assert abs(1 - float(hll2.containment(hll1)) / CONTAINMENT_H2) < ERR_RATE

    assert abs(1 - float(hll1.intersection(hll2)) / INTERSECTION) < ERR_RATE

    """
    hll1.merge(hll2)

    assert abs(1 - float(hllu.similarity(hll1))) < ERR_RATE

    assert abs(1 - float(hllu.containment(hll1))) < ERR_RATE
    assert abs(1 - float(hllu.containment(hll2))) < ERR_RATE

    assert abs(1 - float(hll1.intersection(hllu)) / N_UNIQUE_U) < ERR_RATE
    """

def test_hll_save_load():
    filename = utils.get_test_data('ecoli.genes.fna')
    hll = HLL(ERR_RATE, K)
    n_consumed = n = 0
    for n, record in enumerate(fasta_iter(open(filename)), 1):
        hll.add_sequence(record['sequence'])

    assert abs(1 - float(len(hll)) / N_UNIQUE) < ERR_RATE

    with NamedTemporaryFile() as f:
        hll.save(f.name)

        new_hll = HLL.load(f.name)

    assert len(hll) == len(new_hll)
