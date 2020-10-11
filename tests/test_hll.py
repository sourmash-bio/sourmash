from screed.fasta import fasta_iter
import pytest

from sourmash.hll import HLL

from . import sourmash_tst_utils as utils

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
