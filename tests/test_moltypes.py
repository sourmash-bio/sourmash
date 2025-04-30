import os
import pytest
from collections import namedtuple
import numpy

import sourmash
import sourmash_tst_utils as utils
from sourmash.sourmash_args import load_one_signature
from sourmash.command_sketch import _signatures_for_sketch_factory, ComputeParameters


MoltypeHolder = namedtuple(
    "MoltypeHolder",
    ["moltype", "genome_sketch", "metag_sketch", "cli_moltype_arg", "molecule"],
)


@pytest.fixture(
    scope="session", params=["dna", "protein", "hp", "dayhoff", "skipm1n3", "skipm2n3"]
)
def moltype(request):
    yield request.param


@pytest.fixture(
    scope="session", params=["dna", "protein", "hp", "dayhoff", "skipm1n3", "skipm2n3"]
)
def moltype2(request):
    yield request.param


# build and return genome & metagenome sketches of the given moltype, using
# MoltypeHolder.
@pytest.fixture(scope="session")
def moltype_sketches(runtmp_session, moltype):
    genome = utils.get_test_data("genome-s10.fa.gz")
    metagenome = utils.get_test_data("genome-s10+s11.fa.gz")

    outfile = runtmp_session.output(f"genome.{moltype}.sig.zip")
    outfile2 = runtmp_session.output(f"metagenome.{moltype}.sig.zip")
    assert not os.path.exists(outfile)
    assert not os.path.exists(outfile2)

    # @CTB use match!
    if moltype == "dna":
        runtmp_session.sourmash("sketch", "dna", genome, "-o", outfile)
        runtmp_session.sourmash("sketch", "dna", metagenome, "-o", outfile2)
        mt = MoltypeHolder("dna", outfile, outfile2, "--dna", "DNA")
    elif moltype == "protein":
        runtmp_session.sourmash("sketch", "translate", genome, "-o", outfile)
        runtmp_session.sourmash("sketch", "translate", metagenome, "-o", outfile2)
        mt = MoltypeHolder("protein", outfile, outfile2, "--protein", "protein")
    elif moltype == "hp":
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "hp"
        )
        runtmp_session.sourmash(
            "sketch", "translate", metagenome, "-o", outfile2, "-p", "hp"
        )
        mt = MoltypeHolder("hp", outfile, outfile2, "--hp", "hp")
    elif moltype == "dayhoff":
        runtmp_session.sourmash(
            "sketch", "translate", genome, "-o", outfile, "-p", "dayhoff"
        )
        runtmp_session.sourmash(
            "sketch", "translate", metagenome, "-o", outfile2, "-p", "dayhoff"
        )
        mt = MoltypeHolder("dayhoff", outfile, outfile2, "--dayhoff", "dayhoff")
    elif moltype == "skipm1n3":
        runtmp_session.sourmash(
            "sketch", "dna", genome, "-o", outfile, "-p", "skipm1n3"
        )
        runtmp_session.sourmash(
            "sketch", "dna", metagenome, "-o", outfile2, "-p", "skipm1n3"
        )
        mt = MoltypeHolder("skipm1n3", outfile, outfile2, "--skipm1n3", "skipm1n3")
    elif moltype == "skipm2n3":
        runtmp_session.sourmash(
            "sketch", "dna", genome, "-o", outfile, "-p", "skipm2n3"
        )
        runtmp_session.sourmash(
            "sketch", "dna", metagenome, "-o", outfile2, "-p", "skipm2n3"
        )
        mt = MoltypeHolder("skipm2n3", outfile, outfile2, "--skipm2n3", "skipm2n3")
    else:
        assert 0, "unknown moltype!?"

    assert mt.genome_sketch
    assert mt.metag_sketch
    assert mt.genome_sketch != mt.metag_sketch

    yield (mt, runtmp_session)


#
# all the actual tests
#


def test_factory_moltype(moltype):
    factory = _signatures_for_sketch_factory([moltype], None)
    params_list = list(factory.get_compute_params())
    assert len(params_list) == 1

    params = params_list[0]
    assert getattr(params, moltype)


def test_factory_moltype_equal(moltype):
    factory1 = _signatures_for_sketch_factory([moltype], None)
    params_list1 = list(factory1.get_compute_params())
    assert len(params_list1) == 1
    params1 = params_list1[0]

    factory2 = _signatures_for_sketch_factory([], moltype)
    params_list2 = list(factory2.get_compute_params())
    assert len(params_list2) == 1
    params2 = params_list2[0]

    assert params1 == params2
    assert repr(params1) == repr(params2)

    descr = repr(params1)
    assert f"{moltype}=True" in descr


def test_factory_moltype_ne(moltype, moltype2):
    # make sure that params are different for different moltypes
    # (tests that __eq__ fails properly ;)
    if moltype == moltype2:
        return

    factory1 = _signatures_for_sketch_factory([moltype], None)
    params_list1 = list(factory1.get_compute_params())
    assert len(params_list1) == 1
    params1 = params_list1[0]

    factory2 = _signatures_for_sketch_factory([moltype2], None)
    params_list2 = list(factory2.get_compute_params())
    assert len(params_list2) == 1
    params2 = params_list2[0]

    assert params1 != params2
    assert repr(params1) != repr(params2)


def test_manifest_row_to_compute_parameters(moltype, moltype2):
    if moltype == moltype2:
        return

    # test ComputeParameters.from_manifest_row with moltype
    if moltype == "dna":
        moltype_str = "DNA"
    else:
        moltype_str = moltype

    row = dict(moltype=moltype_str, ksize=21, num=0, scaled=1000, with_abundance=1)
    p = ComputeParameters.from_manifest_row(row)
    assert getattr(p, moltype)
    assert not getattr(p, moltype2)
    assert p.moltype.lower() == moltype
    assert p.num_hashes == 0
    assert p.scaled == 1000
    assert p.track_abundance
    assert p.seed == 42


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
    # test: cat works
    mt, rts = moltype_sketches
    rts.sourmash(
        "sig",
        "cat",
        mt.cli_moltype_arg,
        mt.genome_sketch,
        "-o",
        rts.output(f"sig_cat.out.{mt.moltype}.zip"),
    )


def test_sig_describe(moltype_sketches):
    # test: describe works
    mt, rts = moltype_sketches
    rts.sourmash("sig", "describe", mt.cli_moltype_arg, mt.genome_sketch)
    print(rts.last_result.out)
    assert f"molecule={mt.molecule}" in rts.last_result.out


def test_search(moltype_sketches):
    # test: search finds matches.
    mt, rts = moltype_sketches

    output = rts.output(f"search.{mt.moltype}.csv")
    assert not os.path.exists(output)

    rts.sourmash(
        "search", mt.cli_moltype_arg, mt.genome_sketch, mt.metag_sketch, "-o", output
    )

    with open(output, newline="") as fp:
        x = fp.readlines()
    print(x)
    assert len(x) == 2


def test_compare(moltype_sketches):
    # test: compare finds matches.
    mt, rts = moltype_sketches

    output = rts.output(f"compare.{mt.moltype}.cmp")
    assert not os.path.exists(output)

    rts.sourmash(
        "compare", mt.cli_moltype_arg, mt.genome_sketch, mt.metag_sketch, "-o", output
    )

    with open(output, "rb") as fp:
        arr = numpy.load(fp)
    print(arr)
    assert arr[0, 1] > 0


def test_index(moltype_sketches, disk_index_type):
    # test: 'index' works
    mt, rts = moltype_sketches

    outfile = rts.output(f"index.{mt.moltype}.{disk_index_type}")
    assert not os.path.exists(outfile)

    rts.sourmash(
        "index", mt.cli_moltype_arg, "-F", disk_index_type, outfile, mt.genome_sketch
    )
    print(rts.last_result.err)
    assert "loaded 1 sigs; saving" in rts.last_result.err


def test_gather(moltype_sketches):
    # test: 'gather' works
    mt, rts = moltype_sketches

    output = rts.output(f"gather.{mt.moltype}.csv")
    assert not os.path.exists(output)

    rts.sourmash(
        "gather", mt.cli_moltype_arg, mt.metag_sketch, mt.genome_sketch, "-o", output
    )

    with open(output, newline="") as fp:
        x = fp.readlines()
    print(x)
    assert len(x) == 2


def test_prefetch(moltype_sketches):
    # test: 'prefetch' works
    mt, rts = moltype_sketches

    output = rts.output(f"prefetch.{mt.moltype}.csv")
    assert not os.path.exists(output)

    rts.sourmash(
        "prefetch", mt.cli_moltype_arg, mt.metag_sketch, mt.genome_sketch, "-o", output
    )

    with open(output, newline="") as fp:
        x = fp.readlines()
    print(x)
    assert len(x) == 2


def test_sig_manifest(moltype_sketches):
    # test: 'sig manifest' works
    mt, rts = moltype_sketches

    output = rts.output(f"manifest.{mt.moltype}.csv")
    assert not os.path.exists(output)

    rts.sourmash("sig", "manifest", mt.genome_sketch, "-o", output)

    with open(output, newline="") as fp:
        x = fp.readlines()
    print(x)
    assert len(x) == 3
    assert f",{mt.molecule}," in x[-1]
