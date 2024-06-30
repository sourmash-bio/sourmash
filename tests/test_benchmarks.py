import random
from tempfile import NamedTemporaryFile

import pytest

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
    for _ in range(10):
        random_seq = random.sample(
            "A,C,G,T".split(",") * RANDOM_SEQ_SIZE, RANDOM_SEQ_NUMBER
        )
        sequences.append("".join(random_seq))
    return sequences


@pytest.fixture
def mh():
    return MinHash(MINHASH_NUM, MINHASH_K, track_abundance=False)

@pytest.fixture
def mh_protein():
    return MinHash(
        MINHASH_NUM, MINHASH_K, is_protein=True, track_abundance=False
    )


@pytest.fixture
def sequences():
    return load_sequences()


@pytest.fixture
def populated_mh(sequences):
    populated_mh = MinHash(MINHASH_NUM, MINHASH_K, track_abundance=False)
    for seq in sequences:
        populated_mh.add_sequence(seq)
    return populated_mh


def test_add_sequence(benchmark, mh, sequences):
    @benchmark
    def bench():
        for seq in sequences:
            mh.add_sequence(seq)


def test_add_protein(benchmark, mh_protein, sequences):
    @benchmark
    def bench():
        for seq in sequences:
            mh_protein.add_protein(seq)


def test_get_mins(benchmark, populated_mh):
    benchmark(populated_mh.get_mins)


def test_add_hash(benchmark, mh):
    @benchmark
    def bench():
        for i in range(ADD_HASH_RANGE):
            mh.add_hash(i)


def test_add_many(benchmark, mh):
    benchmark(mh.add_many, list(range(ADD_MANY_RANGE)))


def test_similarity(benchmark, mh, populated_mh):
    benchmark(mh.similarity, populated_mh)

def test_count_common(benchmark, mh, populated_mh):
    benchmark(mh.count_common, populated_mh)


def test_merge(benchmark, mh, populated_mh):
    benchmark(mh.merge, populated_mh)


def test_copy(benchmark, populated_mh):
    benchmark(populated_mh.__copy__)


def test_concat(benchmark, mh, populated_mh):
    benchmark(mh.__iadd__, populated_mh)

####################


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


####################


@pytest.fixture
def zipstore():
    import zipfile

    zf = NamedTemporaryFile()

    with zipfile.ZipFile(
        zf, mode="w", compression=zipfile.ZIP_STORED
    ) as storage:
        for i in range(ZIP_STORAGE_WRITE):
            # just so we have lots of entries
            storage.writestr(str(i), b"0")
        # one big-ish entry
        storage.writestr("sig1", b"9" * 1_000_000)

    yield zf

    zf.close()


def test_load_from_zipstorage(benchmark, zipstore):
    @benchmark
    def bench():
        with ZipStorage(zipstore.name) as storage:
            for _ in range(ZIP_STORAGE_LOAD):
                storage.load("sig1")


def test_load_small_from_zipstorage(benchmark, zipstore):
    @benchmark
    def bench():
        with ZipStorage(zipstore.name) as storage:
            for _ in range(ZIP_STORAGE_LOAD):
                storage.load("99999")
