from __future__ import unicode_literals
import random


try:
    from sourmash._minhash import MinHash
except:
    from sourmash.minhash import MinHash


def load_sequences():
    sequences = []
    for i in range(10):
        random_seq = random.sample("A,C,G,T".split(",") * 3000, 300)
        sequences.append("".join(random_seq))
    return sequences


class TimeMinHashSuite:
    def setup(self):
        self.mh = MinHash(500, 21, track_abundance=False)
        self.sequences = load_sequences()

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
        mh = self.populated_mh
        for i in range(500):
            mh.__copy__()

    def time_concat(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh += other_mh

    # TODO: add_protein


class PeakmemMinHashSuite:
    def setup(self):
        self.mh = MinHash(500, 21, track_abundance=True)
        self.sequences = load_sequences()

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
        mins = self.populated_mh.get_mins(with_abundance=True)
        for i in range(500):
            mh.set_abundances(mins)

class PeakmemMinAbundanceSuite(PeakmemMinHashSuite):
    def setup(self):
        PeakmemMinHashSuite.setup(self)
        self.mh = MinHash(500, 21, track_abundance=True)
