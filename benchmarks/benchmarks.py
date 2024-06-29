import random
from tempfile import NamedTemporaryFile

from sourmash.sbt_storage import ZipStorage
from sourmash.minhash import MinHash

RANDOM_SEQ_SIZE = 3000
RANDOM_SEQ_NUMBER = 300

MINHASH_NUM = 500
MINHASH_K = 21

GET_MINS_RANGE = 500
ADD_HASH_RANGE = 10_000
ADD_MANY_RANGE = 1000
SIMILARITY_TIMES = 500
COUNT_COMMON_TIMES = 500
MERGE_TIMES = 500
COPY_TIMES = 500
CONCAT_TIMES = 500
SET_ABUNDANCES_RANGE = 500
ZIP_STORAGE_WRITE = 100_000
ZIP_STORAGE_LOAD = 20


def load_sequences():
    sequences = []
    for i in range(10):
        random_seq = random.sample(
            "A,C,G,T".split(",") * RANDOM_SEQ_SIZE, RANDOM_SEQ_NUMBER
        )
        sequences.append("".join(random_seq))
    return sequences


class PeakmemMinHashSuite:
    def setup(self):
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)
        self.protein_mh = MinHash(
            MINHASH_NUM, MINHASH_K, is_protein=True, track_abundance=True
        )
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


class PeakmemMinAbundanceSuite(PeakmemMinHashSuite):
    def setup(self):
        PeakmemMinHashSuite.setup(self)
        self.mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=True)


####################


class PeakmemZipStorageSuite:
    def setup(self):
        import zipfile

        self.zipfile = NamedTemporaryFile()

        with zipfile.ZipFile(
            self.zipfile, mode="w", compression=zipfile.ZIP_STORED
        ) as storage:
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
