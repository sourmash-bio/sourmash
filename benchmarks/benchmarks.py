import random
from tempfile import NamedTemporaryFile

from sourmash.sbt_storage import ZipStorage
from sourmash.minhash import MinHash

RANDOM_SEQ_SIZE=3000
RANDOM_SEQ_SAMPLE=300

MINHASH_NUM=500
MINHASH_K=21

GET_MINS_RANGE=500
ADD_HASH_RANGE=10_000
ADD_MANY_RANGE=1000
SIMILARITY_TIMES=500
COUNT_COMMON_TIMES=500
MERGE_TIMES=500
COPY_TIMES=500
CONCAT_TIMES=500
SET_ABUNDANCES_RANGE=500
ZIP_STORAGE_WRITE=100_000
ZIP_STORAGE_LOAD=20


def load_sequences():
    sequences = []
    for i in range(10):
        random_seq = random.sample("A,C,G,T".split(",") * RANDOM_SEQ_SIZE,
                                   RANDOM_SEQ_NUM)
        sequences.append("".join(random_seq))
    return sequences


class TimeMinHashSuite:
    def setup(self):
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=False)
        self.protein_mh = MinHash(MINHASH_NUM, MINHASH_K, is_protein=True,
                                  track_abundance=False)
        self.sequences = load_sequences()

        self.populated_mh = MinHash(MINHASH_NUM, MINHASH_K,
                                    track_abundance=False)
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
        for i in range(GET_MINS_RANGE):
            mh.get_mins()

    def time_add_hash(self):
        mh = self.mh
        for i in range(ADD_HASH_RANGE):
            mh.add_hash(i)

    def time_add_many(self):
        mh = self.mh
        mh.add_many(list(range(ADD_MANY_RANGE)))

    def time_similarity(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(SIMILARITY_TIMES):
            mh.similarity(other_mh)

    def time_count_common(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(COUNT_COMMON_TIMES):
            mh.count_common(other_mh)

    def time_merge(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(MERGE_TIMES):
            mh.merge(other_mh)

    def time_copy(self):
        mh = self.populated_mh
        for i in range(COPY_TIMES):
            mh.__copy__()

    def time_concat(self):
        mh = self.mh
        other_mh = self.populated_mh
        for i in range(CONCAT_TIMES):
            mh += other_mh


class PeakmemMinHashSuite:
    def setup(self):
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)
        self.protein_mh = MinHash(MINHASH_NUM, MINHASH_K,
                                  is_protein=True, track_abundance=True)
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
        for i in range(ADD_HASH_RANGE):
            mh.add_hash(i)

    def peakmem_add_many(self):
        mh = self.mh
        mh.add_many(list(range(ADD_MANY_RANGE)))


####################


class TimeMinAbundanceSuite(TimeMinHashSuite):
    def setup(self):
        TimeMinHashSuite.setup(self)
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)

        self.populated_mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)
        for seq in self.sequences:
            self.populated_mh.add_sequence(seq)

    def time_get_mins_abundance(self):
        mh = self.populated_mh
        for i in range(GET_MINS_RANGE):
            mh.get_mins(with_abundance=True)

    def time_set_abundances(self):
        mh = self.mh
        mins = self.populated_mh.get_mins(with_abundance=True)
        for i in range(SET_ABUNDANCES_RANGE):
            mh.set_abundances(mins)

    def time_set_abundances_noclear(self):
        mh = self.mh
        mins = self.populated_mh.get_mins(with_abundance=True)
        for i in range(SET_ABUNDANCES_RANGE):
            mh.set_abundances(mins, clear=False)

class PeakmemMinAbundanceSuite(PeakmemMinHashSuite):
    def setup(self):
        PeakmemMinHashSuite.setup(self)
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)

####################

class TimeZipStorageSuite:

    def setup(self):
        import zipfile
        self.zipfile = NamedTemporaryFile()

        with zipfile.ZipFile(self.zipfile, mode='w',
                          compression=zipfile.ZIP_STORED) as storage:
            for i in range(ZIP_STORAGE_WRITE):
                # just so we have lots of entries
                storage.writestr(str(i), b"0")
            # one big-ish entry
            storage.writestr("sig1", b"9" * 1_000_000)

    def time_load_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(ZIP_STORAGE_LOAD):
                storage.load("sig1")

    def time_load_small_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(ZIP_STORAGE_LOAD):
                storage.load("99999")

    def teardown(self):
        self.zipfile.close()


class PeakmemZipStorageSuite:
    def setup(self):
        import zipfile
        self.zipfile = NamedTemporaryFile()

        with zipfile.ZipFile(self.zipfile, mode='w',
                          compression=zipfile.ZIP_STORED) as storage:
            for i in range(ZIP_STORAGE_WRITE):
                # just so we have lots of entries
                storage.writestr(str(i), b"0")
            # one big-ish entry
            storage.writestr("sig1", b"9" * 1_000_000)


    def peakmem_load_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(ZIP_STORAGE_LOAD):
                storage.load("sig1")

    def peakmem_load_small_from_zipstorage(self):
        with ZipStorage(self.zipfile.name) as storage:
            for i in range(ZIP_STORAGE_LOAD):
                storage.load("99999")

    def teardown(self):
        self.zipfile.close()
