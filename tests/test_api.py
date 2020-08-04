import pytest
import sourmash

from . import sourmash_tst_utils as utils


def test_sourmash_signature_api():
    e = sourmash.MinHash(n=1, ksize=20)
    sig = sourmash.SourmashSignature(e)

    s = sourmash.save_signatures([sig])
    sig_x1 = sourmash.load_one_signature(s)
    sig_x2 = list(sourmash.load_signatures(s))[0]

    assert sig_x1 == sig
    assert sig_x2 == sig


@utils.in_tempdir
def test_load_index_0_no_file(c):
    with pytest.raises(OSError) as exc:
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


def test_load_fasta_as_signature():
    # try loading a fasta file - should fail with informative exception
    testfile = utils.get_test_data('short.fa')

    with pytest.raises(OSError) as e:
        idx = sourmash.load_file_as_index(testfile)

    assert "Error while reading signatures from '{}' - got sequences instead! Is this a FASTA/FASTQ file?".format(testfile) in str(e)
