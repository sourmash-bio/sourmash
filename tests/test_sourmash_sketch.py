"""
Tests for sourmash sketch command-line functionality.
"""
import os
import gzip
import shutil
import screed
import glob
import json
import csv
import pytest
import screed

import sourmash_tst_utils as utils
import sourmash
from sourmash import MinHash
from sourmash.sbt import SBT, Node
from sourmash.sbtmh import SigLeaf, load_sbt_index
from sourmash.command_compute import ComputeParameters
from sourmash.cli.compute import subparser
from sourmash.cli import SourmashParser

from sourmash import signature
from sourmash import VERSION

###

from sourmash.command_sketch import _signatures_for_sketch_factory
from sourmash_tst_utils import SourmashCommandFailed


def test_do_sourmash_sketch_check_scaled_bounds_negative(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'translate', '-p', 'scaled=-5', testdata1)
    assert "ERROR: scaled value must be positive" in runtmp.last_result.err


def test_do_sourmash_sketch_check_scaled_bounds_less_than_minimum(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'scaled=50', testdata1)
    assert "WARNING: scaled value should be >= 100. Continuing anyway." in runtmp.last_result.err


def test_do_sourmash_sketch_check_scaled_bounds_more_than_maximum(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'scaled=1000000000', testdata1)
    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in runtmp.last_result.err


def test_do_sourmash_sketch_check_num_bounds_negative(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'translate', '-p', 'num=-5', testdata1)
    assert "ERROR: num value must be positive" in runtmp.last_result.err


def test_do_sourmash_sketch_check_num_bounds_less_than_minimum(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'num=25', testdata1)
    assert "WARNING: num value should be >= 50. Continuing anyway." in runtmp.last_result.err


def test_do_sourmash_sketch_check_num_bounds_more_than_maximum(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'num=100000', testdata1)
    assert "WARNING: num value should be <= 50000. Continuing anyway." in runtmp.last_result.err


def test_dna_defaults():
    factory = _signatures_for_sketch_factory([], 'dna', False)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 1
    params = params_list[0]

    assert params.ksizes == [31]
    assert params.num_hashes == 0
    assert params.scaled == 1000
    assert not params.track_abundance
    assert params.seed == 42
    assert params.dna
    assert not params.dayhoff
    assert not params.hp
    assert not params.protein


def test_dna_override_1():
    factory = _signatures_for_sketch_factory(['k=21,scaled=2000,abund'],
                                             'dna', False)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 1
    params = params_list[0]

    assert params.ksizes == [21]
    assert params.num_hashes == 0
    assert params.scaled == 2000
    assert params.track_abundance
    assert params.seed == 42
    assert params.dna
    assert not params.dayhoff
    assert not params.hp
    assert not params.protein


def test_dna_override_bad_1():
    with pytest.raises(ValueError):
        factory = _signatures_for_sketch_factory(['k=21,scaledFOO=2000,abund'],
                                                 'dna', False)


def test_dna_override_bad_2():
    with pytest.raises(ValueError):
        factory = _signatures_for_sketch_factory(['k=21,protein'],
                                                 'dna', False)

def test_protein_defaults():
    factory = _signatures_for_sketch_factory([], 'protein', True)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 1
    params = params_list[0]

    assert params.ksizes == [30]          # x3 for now
    assert params.num_hashes == 0
    assert params.scaled == 200
    assert not params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert not params.dayhoff
    assert not params.hp
    assert params.protein


def test_protein_override_bad_2():
    with pytest.raises(ValueError):
        factory = _signatures_for_sketch_factory(['k=21,dna'],
                                                 'protein', False)

def test_protein_override_bad_rust_foo():
    # mimic 'sourmash sketch protein -p dna'
    factory = _signatures_for_sketch_factory([], 'protein', False)

    # reach in and avoid error checking to construct a bad params_list.
    factory.params_list = [('dna', {})]

    # now, get sigs...
    siglist = factory()
    assert len(siglist) == 1
    sig = siglist[0]

    # try adding something
    testdata1 = utils.get_test_data('ecoli.faa')
    record = next(iter(screed.open(testdata1)))

    with pytest.raises(ValueError) as exc:
        sig.add_protein(record.sequence)

    assert 'Invalid hash function: "dna"' in str(exc)


def test_dayhoff_defaults():
    factory = _signatures_for_sketch_factory([], 'dayhoff', True)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 1
    params = params_list[0]

    assert params.ksizes == [48]          # x3 for now
    assert params.num_hashes == 0
    assert params.scaled == 200
    assert not params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert params.dayhoff
    assert not params.hp
    assert not params.protein


def test_dayhoff_override_bad_2():
    with pytest.raises(ValueError):
        factory = _signatures_for_sketch_factory(['k=21,dna'],
                                                 'dayhoff', False)

def test_hp_defaults():
    factory = _signatures_for_sketch_factory([], 'hp', True)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 1
    params = params_list[0]

    assert params.ksizes == [126]          # x3 for now
    assert params.num_hashes == 0
    assert params.scaled == 200
    assert not params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert not params.dayhoff
    assert params.hp
    assert not params.protein


def test_hp_override_bad_2():
    with pytest.raises(ValueError):
        factory = _signatures_for_sketch_factory(['k=21,dna'],
                                                 'hp', False)


def test_multiple_moltypes():
    params_foo = ['k=20,num=500,protein',
                  'k=19,num=400,dayhoff,abund',
                  'k=30,scaled=200,hp',
                  'k=30,scaled=200,seed=58']
    factory = _signatures_for_sketch_factory(params_foo, 'protein', True)
    params_list = list(factory.get_compute_params())

    assert len(params_list) == 4

    params = params_list[0]
    assert params.ksizes == [60]          # x3, for now.
    assert params.num_hashes == 500
    assert params.scaled == 0
    assert not params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert not params.dayhoff
    assert not params.hp
    assert params.protein

    params = params_list[1]
    assert params.ksizes == [57]          # x3, for now.
    assert params.num_hashes == 400
    assert params.scaled == 0
    assert params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert params.dayhoff
    assert not params.hp
    assert not params.protein

    params = params_list[2]
    assert params.ksizes == [90]          # x3, for now.
    assert params.num_hashes == 0
    assert params.scaled == 200
    assert not params.track_abundance
    assert params.seed == 42
    assert not params.dna
    assert not params.dayhoff
    assert params.hp
    assert not params.protein

    params = params_list[3]
    assert params.ksizes == [90]          # x3, for now.
    assert params.num_hashes == 0
    assert params.scaled == 200
    assert not params.track_abundance
    assert params.seed == 58
    assert not params.dna
    assert not params.dayhoff
    assert not params.hp
    assert params.protein


### command line tests


@utils.in_thisdir
def test_do_sourmash_sketchdna_empty(c):
    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sketch', 'dna')
    assert 'error: no input filenames provided! nothing to do - exiting.' in c.last_result.err


@utils.in_thisdir
def test_do_sourmash_sketchprotein_empty(c):
    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sketch', 'protein')
    assert 'error: no input filenames provided! nothing to do - exiting.' in c.last_result.err


@utils.in_thisdir
def test_do_sourmash_sketchtranslate_empty(c):
    with pytest.raises(SourmashCommandFailed):
        c.run_sourmash('sketch', 'translate')
    assert 'error: no input filenames provided! nothing to do - exiting.' in c.last_result.err


def test_do_sourmash_sketchdna(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'dna', testdata1)

    sigfile = runtmp.output('short.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('short.fa')


def test_do_sourmash_sketchdna_from_file(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    file_list = runtmp.output("filelist.txt")
    with open(file_list, 'wt') as fp:
        print(testdata1, file=fp)

    runtmp.sourmash('sketch', 'dna', '--from-file', file_list)

    sigfile = runtmp.output('short.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('short.fa')


@utils.in_tempdir
def test_do_sourmash_sketchdna_noinput(c):
    data = ""

    cmd = ['sketch', 'dna', '-', '-o', c.output('xxx.sig')]
    c.run_sourmash(*cmd, stdin_data=data)

    sigfile = c.output('xxx.sig')
    assert not os.path.exists(sigfile)
    assert 'no sequences found' in c.last_result.err


@utils.in_tempdir
def test_do_sourmash_sketchdna_noinput_singleton(c):
    data = ""

    cmd = ['sketch', 'dna', '-', '-o', c.output('xxx.sig'), '--singleton']
    c.run_sourmash(*cmd, stdin_data=data)

    sigfile = c.output('xxx.sig')
    assert not os.path.exists(sigfile)
    assert 'no sequences found' in c.last_result.err


@utils.in_tempdir
def test_do_sourmash_sketchdna_noinput_merge(c):
    data = ""

    cmd = ['sketch', 'dna', '-', '-o', c.output('xxx.sig'), '--merge', 'name']
    c.run_sourmash(*cmd, stdin_data=data)

    sigfile = c.output('xxx.sig')
    assert not os.path.exists(sigfile)
    assert 'no sequences found' in c.last_result.err


@utils.in_tempdir
def test_do_sourmash_sketchdna_outdir(c):
    testdata1 = utils.get_test_data('short.fa')
    status, out, err = utils.runscript('sourmash',
                                       ['sketch', 'dna', testdata1,
                                        '--outdir', c.location])

    sigfile = os.path.join(c.location, 'short.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('short.fa')


def test_do_sourmash_sketchdna_output_valid_file(runtmp):
    """ Trigger bug #123 """
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')
    sigfile = runtmp.output('short.fa.sig')

    runtmp.sourmash('sketch', 'dna', '-o', sigfile, testdata1, testdata2, testdata3)

    assert os.path.exists(sigfile)
    assert not runtmp.last_result.out # stdout should be empty

    # is it valid json?
    with open(sigfile, 'r') as f:
        data = json.load(f)

    filesigs = [sig['filename'] for sig in data]
    assert all(testdata in filesigs
                for testdata in (testdata1, testdata2, testdata3))


def test_do_sourmash_sketchdna_output_zipfile(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    outfile = runtmp.output('shorts.zip')

    runtmp.sourmash('sketch', 'dna', '-o', outfile, testdata1, testdata2, testdata3)

    assert os.path.exists(outfile)
    assert not runtmp.last_result.out # stdout should be empty

    sigs = list(sourmash.load_file_as_signatures(outfile))
    assert len(sigs) == 3


def test_do_sourmash_sketchdna_output_stdout_valid(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch', 'dna', '-o', '-',  testdata1, testdata2, testdata3)

    # is it valid json?
    data = json.loads(runtmp.last_result.out)

    filesigs = [sig['filename'] for sig in data]
    assert all(testdata in filesigs
                for testdata in (testdata1, testdata2, testdata3))


@utils.in_tempdir
def test_do_sourmash_sketchdna_output_and_name_valid_file(c):
    # test --merge of multiple input files
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')
    sigfile = c.output('short.fa.sig')

    c.run_sourmash('sketch', 'dna', '-p', 'num=500', '-o', sigfile, '--merge',
                   '"name"', testdata1, testdata2, testdata3)

    assert os.path.exists(sigfile)
    assert 'calculated 1 signature for 4 sequences taken from 3 files' in c.last_result.err

    # is it valid json?
    with open(sigfile, 'r') as f:
        data = json.load(f)

    assert len(data) == 1

    sigfile_merged = c.output('short.all.fa.sig')
    c.run_sourmash('sketch', 'dna', '-p', 'num=500', '-o', sigfile_merged,
                   '--merge', '"name"', testdata1, testdata2, testdata3)

    with open(sigfile_merged, 'r') as f:
        data_merged = json.load(f)

    assert data[0]['signatures'][0]['mins'] == data_merged[0]['signatures'][0]['mins']


@utils.in_tempdir
def test_do_sourmash_sketchdna_output_and_name_valid_file_outdir(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')
    sigfile = os.path.join(c.location, 'short.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('sketch', 'dna', '-o', sigfile,
                       '--merge', '"name"',
                       testdata1, testdata2, testdata3,
                       '--outdir', c.location)

    errmsg = c.last_result.err
    assert "ERROR: --outdir doesn't make sense with -o/--output" in errmsg


def test_do_sourmash_sketchdna_singleton(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'dna', '--singleton', testdata1)

    sigfile = runtmp.output('short.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('shortName')


def test_do_sourmash_sketchdna_name(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'dna', '--merge', 'foo', testdata1, '-o', 'foo.sig')

    sigfile = runtmp.output('foo.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert sig.name == 'foo'
    
    runtmp.sourmash('sketch', 'dna', '--name', 'foo', testdata1, '-o', 'foo2.sig')

    sigfile2 = runtmp.output('foo2.sig')
    assert os.path.exists(sigfile2)

    sig2 = next(signature.load_signatures(sigfile))
    assert sig2.name == 'foo'
    assert sig.name == sig2.name


def test_do_sourmash_sketchdna_name_fail_no_output(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '--merge', 'foo', testdata1)

    assert runtmp.last_result.status == -1


def test_do_sourmash_sketchdna_fail_no_output(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '--merge', 'foo', testdata1)

    assert runtmp.last_result.status == -1
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '--name', 'foo', testdata1)

    assert runtmp.last_result.status == -1


def test_do_sourmash_sketchdna_name_from_first(runtmp):
    testdata1 = utils.get_test_data('short3.fa')
    runtmp.sourmash('sketch', 'dna', '--name-from-first', testdata1)

    sigfile = runtmp.output('short3.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert sig.name == 'firstname'


def test_do_sourmash_sketchdna_multik(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,k=21', testdata1)

    outfile = runtmp.output('short.fa.sig')
    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2
    ksizes = set([ x.minhash.ksize for x in siglist ])
    assert 21 in ksizes
    assert 31 in ksizes


def test_do_sketch_dna_override_protein_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '-p', 'k=7,num=500,protein', testdata1)

    assert runtmp.last_result.status != 0
    assert 'Error creating signatures: Incompatible sketch type' in runtmp.last_result.err


def test_do_sketch_protein_override_dna_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'protein', '-p', 'k=7,num=500,dna', testdata1)

    assert runtmp.last_result.status != 0
    assert 'Error creating signatures: Incompatible sketch type' in runtmp.last_result.err


def test_do_sketch_translate_multik_with_protein(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500', testdata1)

    outfile = runtmp.output('short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes


def test_do_sketch_translate_multik_with_protein_from_file(runtmp):
    testdata1 = utils.get_test_data('short.fa')

    file_list = runtmp.output("filelist.txt")
    with open(file_list, 'wt') as fp:
        print(testdata1, file=fp)

    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500', '--from-file', file_list)

    outfile = runtmp.output('short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes


def test_do_sketch_translate_multik_with_dayhoff(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500', '--dayhoff', testdata1)

    outfile = runtmp.output('short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes
        assert all(x.minhash.dayhoff for x in siglist)


def test_do_sketch_translate_multik_with_hp(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500', '--hp', testdata1)

    outfile = runtmp.output('short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes
        assert all(x.minhash.hp for x in siglist)


@utils.in_tempdir
def test_do_sourmash_sketch_translate_multik_only_protein(c):
    # check sourmash sketch_translate with only protein, no nucl
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500',
                   testdata1)
    outfile = os.path.join(c.location, 'short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes


def test_do_sourmash_sketch_translate_bad_sequences(runtmp):
    """Proper error handling when Ns in dna sequence"""
    testdata1 = utils.get_test_data('short.bad.fa')
    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,k=10,num=500', testdata1)

    outfile = runtmp.output('short.bad.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes


def test_do_sketch_protein_multik_input(runtmp):
    testdata1 = utils.get_test_data('ecoli.faa')
    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,k=10,num=500', testdata1)

    outfile = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes

        moltype = set([ x.minhash.moltype == 'protein'
                        for x in siglist ])
        assert len(moltype) == 1
        assert True in moltype


def test_do_sketch_protein_multik_input_from_file(runtmp):
    testdata1 = utils.get_test_data('ecoli.faa')

    file_list = runtmp.output("filelist.txt")
    with open(file_list, 'wt') as fp:
        print(testdata1, file=fp)

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,k=10,num=500', '--from-file', file_list)

    outfile = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 7 in ksizes
        assert 10 in ksizes

        moltype = set([ x.minhash.moltype == 'protein'
                        for x in siglist ])
        assert len(moltype) == 1
        assert True in moltype


def test_do_sourmash_sketchdna_multik_outfile(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31', testdata1, '-o', outfile)

    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2
    ksizes = set([ x.minhash.ksize for x in siglist ])
    assert 21 in ksizes
    assert 31 in ksizes


def test_do_sourmash_sketchdna_with_scaled_1(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=1', testdata1, '-o', outfile)

    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2

    scaled_vals = [ x.minhash.scaled for x in siglist ]
    assert len(scaled_vals) == 2
    assert set(scaled_vals) == { 1 }


def test_do_sourmash_sketchdna_with_scaled_2(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=2', testdata1, '-o', outfile)

    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2

    max_hashes = [ x.minhash._max_hash for x in siglist ]
    assert len(max_hashes) == 2
    assert set(max_hashes) == set([ int(2**64 /2.) ])


def test_do_sourmash_sketchdna_with_scaled(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=100', testdata1, '-o', outfile)

    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2

    max_hashes = [ x.minhash._max_hash for x in siglist ]
    assert len(max_hashes) == 2
    assert set(max_hashes) == set([ int(2**64 /100.) ])


def test_do_sourmash_sketchdna_with_bad_scaled(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=-1', testdata1, '-o', outfile)

    assert runtmp.last_result.status != 0
    print(runtmp.last_result.err)
    assert 'ERROR: scaled value must be positive' in runtmp.last_result.err

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=1000.5', testdata1, '-o', outfile)

    assert runtmp.last_result.status != 0
    assert "cannot parse scaled='1000.5' as an integer" in runtmp.last_result.err

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,scaled=1000000000', testdata1, '-o', outfile)

    assert runtmp.last_result.status == 0
    print('XXX')
    print(runtmp.last_result.err)
    assert 'WARNING: scaled value should be <= 1e6. Continuing anyway.' in runtmp.last_result.err


def test_do_sketch_with_seed(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    outfile = runtmp.output('FOO.xxx')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,k=31,seed=43', testdata1, '-o', outfile)

    assert os.path.exists(outfile)

    siglist = list(signature.load_signatures(outfile))
    assert len(siglist) == 2

    seeds = [ x.minhash.seed for x in siglist ]
    assert len(seeds) == 2
    assert set(seeds) == set([ 43 ])


def test_do_sourmash_check_protein_comparisons(runtmp):
    # this test checks 2 x 2 protein comparisons with E. coli genes.
    testdata1 = utils.get_test_data('ecoli.faa')

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,num=500', '--singleton', testdata1)

    sig1 = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(sig1)

    testdata2 = utils.get_test_data('ecoli.genes.fna')
    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,num=500', '--singleton', testdata2)

    sig2 = runtmp.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig2)

    # I'm not sure why load_signatures is randomizing order, but ok.
    x = list(signature.load_signatures(sig1))
    sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name)

    x = list(signature.load_signatures(sig2))
    sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name)

    name1 = sig1_aa.name.split()[0]
    assert name1 == 'NP_414543.1'
    name2 = sig2_aa.name.split()[0]
    assert name2 == 'NP_414544.1'
    name3 = sig1_trans.name.split()[0]
    assert name3 == 'gi|556503834:2801-3733'
    name4 = sig2_trans.name.split()[0]
    assert name4 == 'gi|556503834:337-2799'

    print(name1, name3, round(sig1_aa.similarity(sig1_trans), 3))
    print(name2, name3, round(sig2_aa.similarity(sig1_trans), 3))
    print(name1, name4, round(sig1_aa.similarity(sig2_trans), 3))
    print(name2, name4, round(sig2_aa.similarity(sig2_trans), 3))

    assert round(sig1_aa.similarity(sig1_trans), 3) == 0.0
    assert round(sig2_aa.similarity(sig1_trans), 3) == 0.166
    assert round(sig1_aa.similarity(sig2_trans), 3) == 0.174
    assert round(sig2_aa.similarity(sig2_trans), 3) == 0.0


@utils.in_tempdir
def test_do_sourmash_check_knowngood_dna_comparisons(c):
    # this test checks against a known good signature calculated
    # by utils/compute-dna-mh-another-way.py
    testdata1 = utils.get_test_data('ecoli.genes.fna')
    c.run_sourmash('sketch', 'dna', '-p', 'k=21,num=500',
                   '--singleton', testdata1)
    sig1 = c.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1, sig2 = sorted(x, key=lambda x: x.name)

    print(sig1.name)
    print(sig2.name)

    knowngood = utils.get_test_data('benchmark.dna.sig')
    good = list(signature.load_signatures(knowngood))[0]

    assert sig2.similarity(good) == 1.0


@utils.in_tempdir
def test_do_sourmash_check_knowngood_dna_comparisons_use_rna(c):
    # check the rna ; otherwise identical to previous test.
    testdata1 = utils.get_test_data('ecoli.genes.fna')
    c.run_sourmash('sketch', 'rna', '-p', 'k=21,num=500', '--singleton',
                   testdata1)
    sig1 = c.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1, sig2 = sorted(x, key=lambda x: x.name)

    knowngood = utils.get_test_data('benchmark.dna.sig')
    good = list(signature.load_signatures(knowngood))[0]

    assert sig2.similarity(good) == 1.0


def test_do_sourmash_check_knowngood_input_protein_comparisons(runtmp):
    # this test checks against a known good signature calculated
    # by utils/compute-input-prot-another-way.py
    testdata1 = utils.get_test_data('ecoli.faa')

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,num=500', '--singleton', testdata1)

    sig1 = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name)

    knowngood = utils.get_test_data('benchmark.input_prot.sig')
    good_aa = list(signature.load_signatures(knowngood))[0]

    assert sig1_aa.similarity(good_aa) == 1.0


def test_do_sourmash_check_knowngood_protein_comparisons(runtmp):
    # this test checks against a known good signature calculated
    # by utils/compute-prot-mh-another-way.py
    testdata1 = utils.get_test_data('ecoli.genes.fna')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=7,num=500', '--singleton', testdata1)

    sig1 = runtmp.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name)

    knowngood = utils.get_test_data('benchmark.prot.sig')
    good_trans = list(signature.load_signatures(knowngood))[0]

    assert sig2_trans.similarity(good_trans) == 1.0


def test_protein_with_stop_codons(runtmp):
    # compare protein seq with/without stop codons, via cli and also python
    # apis

    testdata1 = utils.get_test_data('ecoli.faa')
    ecoli_seq = [ record.sequence for record in screed.open(testdata1) ]

    # first, via CLI w/o stop codons
    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1', testdata1)
    sig1 = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(sig1)

    x = signature.load_one_signature(sig1)
    cli_mh1 = x.minhash

    # second, via CLI w/stop codons
    ecoli_stop = runtmp.output('ecoli.stop.faa')
    with open(ecoli_stop, 'wt') as fp:
        for seq in ecoli_seq:
            fp.write(f'>seq\n{seq}*\n')

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1', ecoli_stop)
    sig2 = runtmp.output('ecoli.stop.faa.sig')
    assert os.path.exists(sig2)

    x = signature.load_one_signature(sig2)
    cli_mh2 = x.minhash

    # now calculate sketch with MinHash...
    py_mh1 = MinHash(n=0, ksize=7, is_protein=True, scaled=1)
    for seq in ecoli_seq:
        py_mh1.add_protein(seq)

    # now calculate sketch with MinHash and stop codons...
    py_mh2 = MinHash(n=0, ksize=7, is_protein=True, scaled=1)
    for seq in ecoli_seq:
        py_mh2.add_protein(seq + '*')

    # and, last, calculate hashes separately with seq_to_hashes
    h_mh1 = MinHash(n=0, ksize=7, is_protein=True, scaled=1)
    h_mh2 = MinHash(n=0, ksize=7, is_protein=True, scaled=1)

    for seq in ecoli_seq:
        h = h_mh1.seq_to_hashes(seq, is_protein=1)
        h_mh1.add_many(h)

        h = h_mh2.seq_to_hashes(seq + '*', is_protein=1)
        h_mh2.add_many(h)

    # check!
    assert cli_mh1 == py_mh1
    assert cli_mh2 == py_mh2

    assert cli_mh1 == h_mh1
    assert cli_mh2 == h_mh2

    assert cli_mh1.contained_by(cli_mh2) == 1.0
    assert py_mh1.contained_by(cli_mh2) == 1.0
    assert h_mh1.contained_by(h_mh2) == 1.0

    assert cli_mh2.contained_by(cli_mh1) < 1
    assert py_mh2.contained_by(cli_mh1) < 1
    assert h_mh2.contained_by(h_mh1) < 1


def test_hp_with_stop_codons(runtmp):
    # compare hp seq with/without stop codons, via cli and also python
    # apis

    testdata1 = utils.get_test_data('ecoli.faa')
    ecoli_seq = [ record.sequence for record in screed.open(testdata1) ]

    # first, via CLI w/o stop codons
    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1,hp', testdata1)
    sig1 = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(sig1)

    x = signature.load_one_signature(sig1)
    cli_mh1 = x.minhash

    # second, via CLI w/stop codons
    ecoli_stop = runtmp.output('ecoli.stop.faa')
    with open(ecoli_stop, 'wt') as fp:
        for seq in ecoli_seq:
            fp.write(f'>seq\n{seq}*\n')

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1,hp', ecoli_stop)
    sig2 = runtmp.output('ecoli.stop.faa.sig')
    assert os.path.exists(sig2)

    x = signature.load_one_signature(sig2)
    cli_mh2 = x.minhash

    # now calculate sketch with MinHash...
    py_mh1 = MinHash(n=0, ksize=7, hp=True, scaled=1)
    for seq in ecoli_seq:
        py_mh1.add_protein(seq)

    # now calculate sketch with MinHash and stop codons...
    py_mh2 = MinHash(n=0, ksize=7, hp=True, scaled=1)
    for seq in ecoli_seq:
        py_mh2.add_protein(seq + '*')

    # and, last, calculate hashes separately with seq_to_hashes
    h_mh1 = MinHash(n=0, ksize=7, hp=True, scaled=1)
    h_mh2 = MinHash(n=0, ksize=7, hp=True, scaled=1)

    for seq in ecoli_seq:
        h = h_mh1.seq_to_hashes(seq, is_protein=1)
        h_mh1.add_many(h)

        h = h_mh2.seq_to_hashes(seq + '*', is_protein=1)
        h_mh2.add_many(h)

    # check!
    assert cli_mh1 == py_mh1
    assert cli_mh2 == py_mh2

    assert cli_mh1 == h_mh1
    assert cli_mh2 == h_mh2

    assert cli_mh1.contained_by(cli_mh2) == 1.0
    assert py_mh1.contained_by(cli_mh2) == 1.0
    assert h_mh1.contained_by(h_mh2) == 1.0

    assert cli_mh2.contained_by(cli_mh1) < 1
    assert py_mh2.contained_by(cli_mh1) < 1
    assert h_mh2.contained_by(h_mh1) < 1


def test_dayhoff_with_stop_codons(runtmp):
    # compare dayhoff seq with/without stop codons, via cli and also python
    # apis

    testdata1 = utils.get_test_data('ecoli.faa')
    ecoli_seq = [ record.sequence for record in screed.open(testdata1) ]

    # first, via CLI w/o stop codons
    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1,dayhoff', testdata1)
    sig1 = runtmp.output('ecoli.faa.sig')
    assert os.path.exists(sig1)

    x = signature.load_one_signature(sig1)
    cli_mh1 = x.minhash

    # second, via CLI w/stop codons
    ecoli_stop = runtmp.output('ecoli.stop.faa')
    with open(ecoli_stop, 'wt') as fp:
        for seq in ecoli_seq:
            fp.write(f'>seq\n{seq}*\n')

    runtmp.sourmash('sketch', 'protein', '-p', 'k=7,scaled=1,dayhoff', ecoli_stop)
    sig2 = runtmp.output('ecoli.stop.faa.sig')
    assert os.path.exists(sig2)

    x = signature.load_one_signature(sig2)
    cli_mh2 = x.minhash

    # now calculate sketch with MinHash...
    py_mh1 = MinHash(n=0, ksize=7, dayhoff=True, scaled=1)
    for seq in ecoli_seq:
        py_mh1.add_protein(seq)

    # now calculate sketch with MinHash and stop codons...
    py_mh2 = MinHash(n=0, ksize=7, dayhoff=True, scaled=1)
    for seq in ecoli_seq:
        py_mh2.add_protein(seq + '*')

    # and, last, calculate hashes separately with seq_to_hashes
    h_mh1 = MinHash(n=0, ksize=7, dayhoff=True, scaled=1)
    h_mh2 = MinHash(n=0, ksize=7, dayhoff=True, scaled=1)

    for seq in ecoli_seq:
        h = h_mh1.seq_to_hashes(seq, is_protein=1)
        h_mh1.add_many(h)

        h = h_mh2.seq_to_hashes(seq + '*', is_protein=1)
        h_mh2.add_many(h)

    # check!
    assert cli_mh1 == py_mh1
    assert cli_mh2 == py_mh2

    assert cli_mh1 == h_mh1
    assert cli_mh2 == h_mh2

    assert cli_mh1.contained_by(cli_mh2) == 1.0
    assert py_mh1.contained_by(cli_mh2) == 1.0
    assert h_mh1.contained_by(h_mh2) == 1.0

    assert cli_mh2.contained_by(cli_mh1) < 1
    assert py_mh2.contained_by(cli_mh1) < 1
    assert h_mh2.contained_by(h_mh1) < 1
