from __future__ import unicode_literals

from screed.fasta import fasta_iter
from sourmash_lib._minhash import MinHash
from sourmash_lib.sourmash_tst_utils import get_test_data


def load_sequences(filepath):
    sequences = []
    with open(filepath, 'rb') as f:
        for s in fasta_iter(f):
            sequences.append(s['sequence'])
    return sequences


class TimeMinHashSuite:
    def setup(self):
        self.mh = MinHash(500, 21, track_abundance=False)
        self.sequences = load_sequences(get_test_data('ecoli.genes.fna')) * 10

        self.populated_mh = MinHash(500, 21, track_abundance=False)
        for seq in self.sequences:
            self.populated_mh.add_sequence(seq)

    def time_add_sequence(self):
        mh = self.mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_sequence(seq)

    def time_get_mins(self):
        mh = self.populated_mh
        for i in range(500):
            mh.get_mins()

    def time_add_hash(self):
        mh = self.mh
        for i in range(10000):
            mh.add_hash(i)

    def time_compare(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh.compare(other_mh)

    def time_count_common(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh.count_common(other_mh)

    def time_merge(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh.merge(other_mh)

    def time_copy(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh.__copy__(other_mh)

    def time_concat(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh += other_mh

    # TODO: add_protein


class PeakmemMinHashSuite:
    def setup(self):
        self.mh = MinHash(500, 21, track_abundance=True)
        self.sequences = load_sequences(get_test_data('ecoli.genes.fna'))

    def peakmem_add_sequence(self):
        mh = self.mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_sequence(seq)


####################


class TimeMinAbundanceSuite(TimeMinHashSuite):
    def setup(self):
        TimeMinHashSuite.setup(self)
        self.mh = MinHash(500, 21, track_abundance=True)

        self.populated_mh = MinHash(500, 21, track_abundance=True)
        for seq in self.sequences:
            self.populated_mh.add_sequence(seq)

    def time_get_mins_abundance(self):
        mh = self.populated_mh
        for i in range(500):
            mh.get_mins(with_abundance=True)

    def time_set_abundances(self):
        mh = self.mh
        mins = self.populated_mh.get_mins()
        for i in range(500):
            mh.set_abundances(mins)

class PeakmemMinAbundanceSuite(PeakmemMinHashSuite):
    def setup(self):
        PeakmemMinHashSuite.setup(self)
        self.mh = MinHash(500, 21, track_abundance=True)
