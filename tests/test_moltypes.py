import pytest
from collections import namedtuple

import sourmash
import sourmash_tst_utils as utils
from sourmash.sourmash_args import load_one_signature


MoltypeHolder = namedtuple(
    "MoltypeHolder",
    ["moltype_descr", "genome_sketch", "metag_sketch", "cli_moltype_arg", "molecule"],
)


#@pytest.fixture(scope="session", params=["dna", "protein", "hp", "dayhoff", "skipm1n3"])
@pytest.fixture(scope="session", params=["dna", "protein", "hp", "dayhoff"])
def moltype(request):
    yield request.param


# build and return genome & metagenome sketches of the given moltype
@pytest.fixture(scope="session")
def moltype_sketches(runtmp_session, moltype):
    genome = utils.get_test_data("genome-s10.fa.gz")
    metagenome = utils.get_test_data("genome-s10+s11.fa.gz")

    outfile = runtmp_session.output(f"genome.{moltype}.sig.zip")
    outfile2 = runtmp_session.output(f"metagenome.{moltype}.sig.zip")
        
    # @CTB use match!
    if moltype == "dna":
        runtmp_session.sourmash("sketch", "dna", genome, "-o", outfile)
        runtmp_session.sourmash("sketch", "dna", metagenome, "-o", outfile2)
        mt = MoltypeHolder("dna", outfile, outfile2, "--dna", "DNA")
    elif moltype == "protein":
        runtmp_session.sourmash("sketch", "translate", genome, "-o", outfile)
        runtmp_session.sourmash("sketch", "translate", metagenome, "-o", outfile2)
        mt = MoltypeHolder("dna", outfile, outfile2, "--protein", "protein")
    elif moltype == "hp":
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "hp"
        )
        runtmp_session.sourmash(
            "sketch", "translate", metagenome, "-o", outfile2, "-p", "hp"
        )
        mt = MoltypeHolder("dna", outfile, outfile2, "--hp", "hp")
    elif moltype == "dayhoff":
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "dayhoff"
        )
        runtmp_session.sourmash(
            "sketch", "translate", metagenome, "-o", outfile2, "-p", "dayhoff"
        )
        mt = MoltypeHolder("dna", outfile, outfile2, "--dayhoff", "dayhoff")
    elif moltype == "skipm1n3":
        runtmp_session.sourmash(
            "sketch", "dna", genome, "-o", outfile, "-p", "skipm1n3"
        )
        runtmp_session.sourmash(
            "sketch", "translate", metagenome, "-o", outfile2, "-p", "skipm1n3"
        )
        mt = MoltypeHolder("dna", outfile, "", "--skipm1n3", "skipm1n3")
    else:
        assert 0

    assert mt.genome_sketch
    assert mt.metag_sketch
    assert mt.genome_sketch != mt.metag_sketch

    yield (mt, runtmp_session)


def test_api_load(moltype_sketches):
    # can we load exactly one sketch? yay.
    mt, rts = moltype_sketches

    gsig = load_one_signature(mt.genome_sketch, select_moltype=mt.molecule)
    msig = load_one_signature(mt.metag_sketch, select_moltype=mt.molecule)

    assert gsig.minhash.moltype == mt.molecule
    assert msig.minhash.moltype == mt.molecule


def test_api_overlap(moltype_sketches):
    # test basic overlap calculations
    mt, rts = moltype_sketches

    gsig = load_one_signature(mt.genome_sketch, select_moltype=mt.molecule)
    msig = load_one_signature(mt.metag_sketch, select_moltype=mt.molecule)

    mh1 = gsig.minhash
    mh2 = msig.minhash

    assert mh1.contained_by(mh2) == 1.0
    assert mh2.contained_by(mh1) > 0
    assert mh1.jaccard(mh2) > 0
    assert mh1.jaccard(mh2) < 1


def test_sig_cat(moltype_sketches):
    mt, rts = moltype_sketches
    rts.sourmash(
        "sig",
        "cat",
        mt.cli_moltype_arg,
        mt.genome_sketch,
        "-o",
        rts.output("sig_cat.out.zip"),
    )


def test_sig_describe(moltype_sketches):
    mt, rts = moltype_sketches
    rts.sourmash("sig", "describe", mt.cli_moltype_arg, mt.genome_sketch)
    print(rts.last_result.out)
    assert f"molecule={mt.molecule}" in rts.last_result.out


# def test_sketch_file(runtmp, moltype):
#    genome = utils.get_test_data('genome-s10+s11.fa.gz')
#    runtmp.sourmash('sketch', moltype, genome, '-o', f'xxx.{moltype}.sig.zip')
#
#    print(moltype)
#    assert 0
