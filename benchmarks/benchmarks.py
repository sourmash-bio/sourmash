import os
import random
from pathlib import Path
from tempfile import NamedTemporaryFile


from sourmash.sbt_storage import ZipStorage
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
        self.protein_mh = MinHash(500, 21, is_protein=True, track_abundance=False)
        self.sequences = load_sequences()

        self.populated_mh = MinHash(500, 21, track_abundance=False)
        for seq in self.sequences:
            self.populated_mh.add_sequence(seq)

    def time_add_sequence(self):
        mh = self.mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_sequence(seq)

    def time_add_protein(self):
        mh = self.protein_mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_protein(seq)

    def time_get_mins(self):
        mh = self.populated_mh
        for i in range(500):
            mh.get_mins()

    def time_add_hash(self):
        mh = self.mh
        for i in range(10000):
            mh.add_hash(i)

    def time_add_many(self):
        mh = self.mh
        mh.add_many(list(range(1000)))

    def time_similarity(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(500):
            mh.similarity(other_mh)

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


class PeakmemMinHashSuite:
    def setup(self):
        self.mh = MinHash(500, 21, track_abundance=True)
        self.protein_mh = MinHash(500, 21, is_protein=True, track_abundance=True)
        self.sequences = load_sequences()

    def peakmem_add_sequence(self):
        mh = self.mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_sequence(seq)

    def peakmem_add_protein(self):
        mh = self.protein_mh
        sequences = self.sequences
        for seq in sequences:
            mh.add_protein(seq)

    def peakmem_add_hash(self):
        mh = self.mh
        for i in range(10000):
            mh.add_hash(i)

    def peakmem_add_many(self):
        mh = self.mh
        mh.add_many(list(range(1000)))


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

    def time_set_abundances_noclear(self):
        mh = self.mh
        mins = self.populated_mh.get_mins(with_abundance=True)
        for i in range(500):
            mh.set_abundances(mins, clear=False)

class PeakmemMinAbundanceSuite(PeakmemMinHashSuite):
    def setup(self):
        PeakmemMinHashSuite.setup(self)
        self.mh = MinHash(500, 21, track_abundance=True)

####################

class TimeZipStorageSuite:

    def setup(self):
        import zipfile
        self.zipfile = NamedTemporaryFile()

        with zipfile.ZipFile(self.zipfile, mode='w',
                          compression=zipfile.ZIP_STORED) as storage:
            for i in range(100_000):
                # just so we have lots of entries
                storage.writestr(str(i), b"0")
            # one big-ish entry
            storage.writestr("sig1", b"9" * 1_000_000)

    def time_load_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(20):
                storage.load("sig1")

    def time_load_small_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(20):
                storage.load("99999")

    def teardown(self):
        self.zipfile.close()


class PeakmemZipStorageSuite:
    def setup(self):
        import zipfile
        self.zipfile = NamedTemporaryFile()

        with zipfile.ZipFile(self.zipfile, mode='w',
                          compression=zipfile.ZIP_STORED) as storage:
            for i in range(100_000):
                # just so we have lots of entries
                storage.writestr(str(i), b"0")
            # one big-ish entry
            storage.writestr("sig1", b"9" * 1_000_000)


    def peakmem_load_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(20):
                storage.load("sig1")

    def peakmem_load_small_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(20):
                storage.load("99999")

    def teardown(self):
        self.zipfile.close()
