import pytest
import sourmash_tst_utils as utils
from collections import namedtuple

MoltypeHolder = namedtuple(
    "MoltypeHolder",
    ["moltype_descr", "genome_sketch", "metag_sketch", "cli_moltype_arg", "molecule"],
)


@pytest.fixture(scope="session", params=["dna", "protein", "hp", "dayhoff", "skipm1n3"])
def moltype(request):
    yield request.param


# build and return genome & metagenome sketches of the given moltype
@pytest.fixture(scope="session")
def moltype_sketches(runtmp_session, moltype):
    genome = utils.get_test_data("genome-s10+s11.fa.gz")

    # @CTB use match!
    if moltype == "dna":
        outfile = "genome.dna.sig.zip"
        runtmp_session.sourmash("sketch", "dna", genome, "-o", outfile)
        mt = MoltypeHolder("dna", outfile, "", "--dna", "DNA")
    elif moltype == "protein":
        outfile = f"genome.{moltype}.sig.zip"
        runtmp_session.sourmash("sketch", "translate", genome, "-o", outfile)
        mt = MoltypeHolder("dna", outfile, "", "--protein", "protein")
    elif moltype == "hp":
        outfile = f"genome.{moltype}.sig.zip"
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "hp"
        )
        mt = MoltypeHolder("dna", outfile, "", "--hp", "hp")
    elif moltype == "dayhoff":
        outfile = f"genome.{moltype}.sig.zip"
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "dayhoff"
        )
        mt = MoltypeHolder("dna", outfile, "", "--dayhoff", "dayhoff")
    elif moltype == "skipm1n3":
        outfile = f"genome.{moltype}.sig.zip"
        runtmp_session.sourmash(
            "sketch", "dna", genome, "-o", outfile, "-p", "skipm1n3"
        )
        mt = MoltypeHolder("dna", outfile, "", "--skipm1n3", "skipm1n3")
    else:
        assert 0

    yield (mt, runtmp_session)


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
