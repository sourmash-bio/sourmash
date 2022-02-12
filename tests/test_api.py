import pytest
import sourmash

import sourmash_tst_utils as utils


@utils.in_tempdir
def test_sourmash_signature_api(c):
    e = sourmash.MinHash(n=1, ksize=20)
    sig = sourmash.SourmashSignature(e)

    with open(c.output('xxx.sig'), 'wt') as fp:
        sourmash.save_signatures([sig], fp)
    sig_x1 = sourmash.load_one_signature(c.output('xxx.sig'))
    sig_x2 = list(sourmash.load_file_as_signatures(c.output('xxx.sig')))[0]

    assert sig_x1 == sig
    assert sig_x2 == sig


@utils.in_tempdir
def test_load_index_0_no_file(c):
    with pytest.raises(ValueError) as exc:
        idx = sourmash.load_file_as_index(c.output('does-not-exist'))
    assert 'Error while reading signatures from ' in str(exc.value)


def test_load_index_1():
    testfile = utils.get_test_data('prot/protein.sbt.zip')
    idx = sourmash.load_file_as_index(testfile)

    sigs = list(idx.signatures())
    assert len(sigs) == 2


def test_load_index_2():
    testfile = utils.get_test_data('prot/protein.lca.json.gz')
    idx = sourmash.load_file_as_index(testfile)

    sigs = list(idx.signatures())
    assert len(sigs) == 2


def test_load_index_3():
    testfile = utils.get_test_data('prot/protein/')
    idx = sourmash.load_file_as_index(testfile)

    sigs = list(idx.signatures())
    assert len(sigs) == 2


def test_load_index_4():
    testfile = utils.get_test_data('prot/all.zip')
    idx = sourmash.load_file_as_index(testfile)

    sigs = list(idx.signatures())
    assert len(sigs) == 8


def test_load_index_4_b():
    testfile = utils.get_test_data('prot/protein.zip')
    idx = sourmash.load_file_as_index(testfile)

    sigs = list(idx.signatures())
    assert len(sigs) == 2


def test_load_fasta_as_signature():
    # try loading a fasta file - should fail with informative exception
    testfile = utils.get_test_data('short.fa')

    with pytest.raises(ValueError) as exc:
        idx = sourmash.load_file_as_index(testfile)

    print(exc.value)

    assert f"Error while reading signatures from '{testfile}' - got sequences instead! Is this a FASTA/FASTQ file?" in str(exc.value)


def test_load_and_search_sbt_api():
    treefile = utils.get_test_data('prot/protein.sbt.zip')
    queryfile = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')

    tree = sourmash.load_sbt_index(treefile)
    query = sourmash.load_one_signature(queryfile)

    results = list(sourmash.search_sbt_index(tree, query, 0))
    assert len(results) == 2
