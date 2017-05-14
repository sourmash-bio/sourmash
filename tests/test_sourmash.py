"""
Tests for the 'sourmash' command line.
"""
from __future__ import print_function, unicode_literals
import os
import gzip
import shutil
import time
import screed
import glob

from . import sourmash_tst_utils as utils
from sourmash_lib import MinHash
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass

from sourmash_lib import signature
from sourmash_lib import VERSION

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_badcmd():
    status, out, err = utils.runscript('sourmash', ['foobarbaz'], fail_ok=True)
    assert status != 0                    # bad arg!
    assert "Unrecognized command" in err

def test_sourmash_info():
    status, out, err = utils.runscript('sourmash', ['info'], fail_ok=False)

    # no output to stdout
    assert not out
    assert "sourmash version" in err
    assert "loaded from path" in err
    assert VERSION in err


def test_sourmash_info_verbose():
    status, out, err = utils.runscript('sourmash', ['info', '-v'])

    # no output to stdout
    assert not out
    assert "khmer version" in err
    assert "screed version" in err
    assert "loaded from path" in err


def test_do_sourmash_compute():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name().endswith('short.fa')


def test_do_sourmash_compute_output_valid_file():
    """ Trigger bug #123 """
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        sigfile = os.path.join(location, 'short.fa.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-o', sigfile,
                                            testdata1,
                                            testdata2, testdata3],
                                           in_directory=location)

        assert os.path.exists(sigfile)

        # is it valid json?
        import json
        with open(sigfile, 'r') as f:
            data = json.load(f)

        filesigs = [sig['filename'] for sig in data]
        assert all(testdata in filesigs
                   for testdata in (testdata1, testdata2, testdata3))


def test_do_sourmash_compute_output_and_name_valid_file():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        sigfile = os.path.join(location, 'short.fa.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-o', sigfile,
                                            '--merge', '"name"',
                                            testdata1,
                                            testdata2, testdata3],
                                           in_directory=location)

        assert os.path.exists(sigfile)

        # is it valid json?
        import json
        with open(sigfile, 'r') as f:
            data = json.load(f)

        assert len(data) == 1

        all_testdata = " ".join([testdata1, testdata2, testdata3])
        sigfile_merged = os.path.join(location, 'short.all.fa.sig')
        #cmd = "cat {} | {}/sourmash compute -o {} -".format(
        cmd = "cat {} | {}/sourmash compute -o {} -".format(
                all_testdata, utils.scriptpath(), sigfile_merged)
        status, out, err = utils.run_shell_cmd(cmd, in_directory=location)

        with open(sigfile_merged, 'r') as f:
            data_merged = json.load(f)

        assert data[0]['signatures'][0]['mins'] == data_merged[0]['signatures'][0]['mins']


def test_do_sourmash_compute_singleton():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--singleton',
                                            testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name().endswith('shortName')


def test_do_sourmash_compute_name():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--merge', 'foo',
                                            testdata1, '-o', 'foo.sig'],
                                           in_directory=location)

        sigfile = os.path.join(location, 'foo.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name() == 'foo'

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--name', 'foo',
                                            testdata1, '-o', 'foo2.sig'],
                                           in_directory=location)

        sigfile2 = os.path.join(location, 'foo2.sig')
        assert os.path.exists(sigfile)

        sig2 = next(signature.load_signatures(sigfile))
        assert sig2.name() == 'foo'
        assert sig.name() == sig2.name()


def test_do_sourmash_compute_name_fail_no_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--merge', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1


def test_do_sourmash_compute_merge_fail_no_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--merge', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--name', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1


def test_do_sourmash_compute_name_from_first():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--name-from-first',
                                            testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short3.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name() == 'firstname'


def test_do_sourmash_compute_multik():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 31 in ksizes


def test_do_sourmash_compute_multik_with_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 4
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes


def test_do_sourmash_compute_multik_with_nothing():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--no-protein', '--no-dna',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        outfile = os.path.join(location, 'short.fa.sig')
        assert not os.path.exists(outfile)


def test_do_sourmash_compute_multik_protein_bad_ksize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '20,32',
                                            '--protein', '--no-dna',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        outfile = os.path.join(location, 'short.fa.sig')
        assert not os.path.exists(outfile)
        assert 'protein ksizes must be divisible by 3' in err


def test_do_sourmash_compute_multik_only_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--protein', '--no-dna',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes


def test_do_sourmash_compute_multik_input_is_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--input-is-protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes

            moltype = set([ x.minhash.is_molecule_type('protein')
                            for x in siglist ])
            assert len(moltype) == 1
            assert True in moltype


def test_do_sourmash_compute_multik_outfile():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            testdata1, '-o', outfile],
                                           in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 31 in ksizes


def test_do_sourmash_compute_with_scaled_1():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '1',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == { 0 }


def test_do_sourmash_compute_with_scaled_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '2',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == set([ int(2**64 /2.) ])


def test_do_sourmash_compute_with_scaled():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '100',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == set([ int(2**64 /100.) ])


def test_do_sourmash_compute_with_seed():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--seed', '43',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        seeds = [ x.minhash.seed for x in siglist ]
        assert len(seeds) == 2
        assert set(seeds) == set([ 43 ])


def test_do_sourmash_check_protein_comparisons():
    # this test checks 2 x 2 protein comparisons with E. coli genes.
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--input-is-protein',
                                            '--singleton',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(sig1)

        testdata2 = utils.get_test_data('ecoli.genes.fna')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--protein', '--no-dna',
                                            '--singleton',
                                            testdata2],
                                           in_directory=location)
        sig2 = os.path.join(location, 'ecoli.genes.fna.sig')
        assert os.path.exists(sig2)

        # I'm not sure why load_signatures is randomizing order, but ok.
        x = list(signature.load_signatures(sig1))
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name())

        x = list(signature.load_signatures(sig2))
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name())

        name1 = sig1_aa.name().split()[0]
        assert name1 == 'NP_414543.1'
        name2 = sig2_aa.name().split()[0]
        assert name2 == 'NP_414544.1'
        name3 = sig1_trans.name().split()[0]
        assert name3 == 'gi|556503834:2801-3733'
        name4 = sig2_trans.name().split()[0]
        assert name4 == 'gi|556503834:337-2799'

        print(name1, name3, round(sig1_aa.similarity(sig1_trans), 3))
        print(name2, name3, round(sig2_aa.similarity(sig1_trans), 3))
        print(name1, name4, round(sig1_aa.similarity(sig2_trans), 3))
        print(name2, name4, round(sig2_aa.similarity(sig2_trans), 3))

        assert round(sig1_aa.similarity(sig1_trans), 3) == 0.0
        assert round(sig2_aa.similarity(sig1_trans), 3) == 0.166
        assert round(sig1_aa.similarity(sig2_trans), 3) == 0.174
        assert round(sig2_aa.similarity(sig2_trans), 3) == 0.0


def test_do_sourmash_check_knowngood_dna_comparisons():
    # this test checks against a known good signature calculated
    # by utils/compute-dna-mh-another-way.py
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.genes.fna')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--singleton', '--dna',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.genes.fna.sig')
        assert os.path.exists(sig1)

        x = list(signature.load_signatures(sig1))
        sig1, sig2 = sorted(x, key=lambda x: x.name())

        knowngood = utils.get_test_data('benchmark.dna.sig')
        good = list(signature.load_signatures(knowngood))[0]

        assert sig2.similarity(good) == 1.0


def test_do_sourmash_check_knowngood_input_protein_comparisons():
    # this test checks against a known good signature calculated
    # by utils/compute-input-prot-another-way.py
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--input-is-protein',
                                            '--singleton',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(sig1)

        x = list(signature.load_signatures(sig1))
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name())

        knowngood = utils.get_test_data('benchmark.input_prot.sig')
        good_aa = list(signature.load_signatures(knowngood))[0]

        assert sig1_aa.similarity(good_aa) == 1.0


def test_do_sourmash_check_knowngood_protein_comparisons():
    # this test checks against a known good signature calculated
    # by utils/compute-prot-mh-another-way.py
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.genes.fna')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--singleton', '--protein',
                                            '--no-dna',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.genes.fna.sig')
        assert os.path.exists(sig1)

        x = list(signature.load_signatures(sig1))
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name())

        knowngood = utils.get_test_data('benchmark.prot.sig')
        good_trans = list(signature.load_signatures(knowngood))[0]

        assert sig2_trans.similarity(good_trans) == 1.0


def test_do_plot_comparison():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash', ['plot', 'cmp'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


def test_do_plot_comparison_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--pdf'],
                                           in_directory=location)
        assert os.path.exists(os.path.join(location, "cmp.dendro.pdf"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.pdf"))


def test_do_plot_comparison_3():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--labels'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


def test_search():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in err
        assert '0.930' in out


def test_search_gzip():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        data = open(os.path.join(location, 'short.fa.sig'), 'rb').read()
        with gzip.open(os.path.join(location, 'zzz.gz'), 'wb') as fp:
            fp.write(data)

        data = open(os.path.join(location, 'short2.fa.sig'), 'rb').read()
        with gzip.open(os.path.join(location, 'yyy.gz'), 'wb') as fp:
            fp.write(data)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'zzz.gz',
                                            'yyy.gz'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in err
        assert '0.930' in out


def test_search_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            testdata3],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig', 'short3.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '2 matches' in err
        assert '0.930' in out
        assert '0.896' in out


def test_search_3():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            testdata3],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', '-n', '1',
                                            'short.fa.sig',
                                            'short2.fa.sig', 'short3.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '2 matches; showing first 1' in err


def test_mash_csv_to_sig():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa.msh.dump')
        testdata2 = utils.get_test_data('short.fa')

        status, out, err = utils.runscript('sourmash', ['import_csv',
                                                        testdata1,
                                                        '-o', 'xxx.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            '-n', '970', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '31',
                                            'short.fa.sig', 'xxx.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches:' in err
        assert 'short.fa \t 1.000 \t xxx.sig' in out


def test_do_sourmash_sbt_index_bad_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig',
                                            '--dna', '--protein'],
                                           in_directory=location, fail_ok=True)

        assert "cannot specify both --dna and --protein!" in err
        assert status != 0


def test_do_sourmash_sbt_search():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)
        print(out)

        assert testdata1 in out
        assert testdata2 in out


def test_do_sourmash_sbt_index_single():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)
        print(out)

        assert testdata1 in out


def test_do_sourmash_sbt_search_selectprot():
    # sbt_index should fail when run on signatures with multiple types
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        args = ['compute', testdata1, testdata2,
                '--protein', '--dna', '-k', '30']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['sbt_index', 'zzz', 'short.fa.sig', 'short2.fa.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)

        print(out)
        print(err)
        assert status != 0


def test_do_sourmash_sbt_search_dnaprotquery():
    # sbt_search should fail if non-single query sig given
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        args = ['compute', testdata1, testdata2,
                '--protein', '--dna', '-k', '30']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['sbt_index', 'zzz', '--protein', 'short.fa.sig',
                'short2.fa.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        args = ['sbt_search', 'zzz', 'short.fa.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)
        assert status != 0
        print(out)
        print(err)
        assert 'need exactly one' in err


def test_do_sourmash_sbt_index_traverse():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            '--traverse-dir', '.'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))
        assert 'loaded 2 sigs; saving SBT under' in err

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)
        print(out)

        assert testdata1 in out
        assert testdata2 in out


def test_do_sourmash_sbt_combine():
    with utils.TempDirectory() as location:
        files = [utils.get_test_data(f) for f in utils.SIG_FILES]

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz'] + files,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_combine', 'joined',
                                            'zzz.sbt.json', 'zzz.sbt.json'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'joined.sbt.json'))

        filename = os.path.splitext(os.path.basename(utils.SIG_FILES[0]))[0]

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', 'zzz'] + [files[0]],
                                           in_directory=location)
        print(out)

        assert out.count(filename) == 1

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', 'joined'] + [files[0]],
                                           in_directory=location)
        print(out)

        # TODO: signature is loaded more than once,
        # so checking if we get two results.
        # If we ever start reporting only one match (even if appears repeated),
        # change this test too!
        assert out.count(filename) == 2


def test_do_sourmash_sbt_search_otherdir():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'xxx/zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'xxx', 'zzz.sbt.json'))

        sbt_name = os.path.join(location,'xxx', 'zzz',)
        sig_loc = os.path.join(location, 'short.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', sbt_name, sig_loc])
        print(out)

        assert testdata1 in out
        assert testdata2 in out


def test_do_sourmash_sbt_search_bestonly():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_search', '--best-only', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)
        print(out)

        assert testdata1 in out


def test_compare_with_abundance_1():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature('', E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature('', E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '1.000' in out


def test_compare_with_abundance_2():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature('', E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature('', E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '1.0' in out


def test_compare_with_abundance_3():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')
        E1.add_sequence('GGACA')

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature('', E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature('', E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '0.705' in out


def test_sbt_gather():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_gather', 'zzz',
                                            'query.fa.sig', '--csv',
                                            'foo.csv'],
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found: 1.00 1.00 ' in err


def test_sbt_gather_file_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_gather', 'zzz',
                                            'query.fa.sig',
                                            '-o', 'foo.out'],
                                           in_directory=location)

        with open(os.path.join(location, 'foo.out')) as f:
            output = f.read()
            print(output)
            assert '1.00 1.00 ' in output


def test_sbt_gather_metagenome():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['sbt_index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_gather', 'gcf_all',
                                            query_sig, '-k', '21'],
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 11 matches total' in err
        assert 'the recovered matches hit 100.0% of the query' in err


def test_sbt_gather_error_no_cardinality_query():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata3],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_gather', 'zzz',
                                            'short3.fa.sig'],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1
        assert "query signature needs to be created with --scaled" in err


def test_sbt_categorize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['categorize', 'zzz', '--traverse-directory', '.',
                '--ksize', '21', '--dna', '--csv', 'out.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)

        # mash dist genome-s10.fa.gz genome-s10+s11.fa.gz
        # yields 521/1000 ==> ~0.5
        assert 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err

        out_csv = open(os.path.join(location, 'out.csv')).read()
        assert './4.sig,genome-s10.fa.gz,0.50' in out_csv


def test_sbt_categorize_already_done():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'in.csv'), 'wt') as fp:
            fp.write('./4.sig,genome-s10.fa.gz,0.50')

        args = ['categorize', 'zzz', './2.sig', './4.sig',
                '--ksize', '21', '--dna', '--load-csv', 'in.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'for genome-s11.fa.gz, no match found'
        assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err


def test_sbt_categorize_already_done_traverse():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'in.csv'), 'wt') as fp:
            fp.write('./4.sig,genome-s10.fa.gz,0.50')

        args = ['categorize', 'zzz', '--traverse-directory', '.',
                '--ksize', '21', '--dna', '--load-csv', 'in.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'for genome-s11.fa.gz, no match found'
        assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err


def test_sbt_categorize_multiple_ksizes_moltypes():
    # 'categorize' should fail when there are multiple ksizes or moltypes
    # present
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))

        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['categorize', 'zzz', '--traverse-directory', '.']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)

        assert status != 0
        assert 'multiple k-mer sizes/molecule types present' in err


def test_watch():
    with utils.TempDirectory() as location:
        testdata0 = utils.get_test_data('genome-s10.fa.gz')
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))

        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        cmd = """

             gunzip -c {} | {}/sourmash watch --ksize 21 --dna zzz

        """.format(testdata0, utils.scriptpath())
        status, out, err = utils.run_shell_cmd(cmd, in_directory=location)

        print(out)
        print(err)
        assert 'FOUND: genome-s10.fa.gz, at 1.000' in err


def test_watch_coverage():
    with utils.TempDirectory() as location:
        testdata0 = utils.get_test_data('genome-s10.fa.gz')
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))

        args = ['sbt_index', '--dna', '-k', '21', 'zzz', '1.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'query.fa'), 'wt') as fp:
            record = list(screed.open(testdata0))[0]
            for start in range(0, len(record), 100):
                fp.write('>{}\n{}\n'.format(start,
                                            record.sequence[start:start+500]))

        args = ['watch', '--ksize', '21', '--dna', 'zzz', 'query.fa']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'FOUND: genome-s10.fa.gz, at 1.000' in err


def test_mash_yaml_to_json():
    with utils.TempDirectory() as location:
        orig_sig = utils.get_test_data('genome-s10.fa.gz.sig')
        shutil.copy(orig_sig, location)
        test_sig = os.path.join(location, os.path.basename(orig_sig))

        # create directory
        os.mkdir(os.path.join(location, "foo"))
        shutil.copy(orig_sig, os.path.join(location, "foo"))

        assert not os.path.exists(test_sig + ".json")
        status, out, err = utils.runscript('sourmash', ['convert',
                                                        test_sig,
                                                        os.path.join(location, "foo")],
                                           in_directory=location)
        # check success
        assert status == 0
        # check existence of JSON files
        assert os.path.exists(test_sig + ".json")
        assert os.path.exists(os.path.join(location, "foo", os.path.basename(orig_sig)) + ".json")

        # check that the files can be read (as JSON)
        with open(test_sig + ".json") as fh:
            sig = signature.signature_json.load_signatures_json(fh)
        with open(os.path.join(location, "foo", os.path.basename(orig_sig)) + ".json") as fh:
            sig = signature.signature_json.load_signatures_json(fh)

        # try again: will fail because .json already found
        status, out, err = utils.runscript('sourmash', ['convert',
                                                        test_sig],
                                           in_directory=location,
                                           fail_ok=True)

        assert status == 1

        timestamp = os.path.getmtime(test_sig + ".json")

        # briefly sleep to make sure the clock ticks at least once
        # between the two timestamps
        time.sleep(1)

        # try again: will not fail when .json already found because of --force
        time.sleep(1)
        status, out, err = utils.runscript('sourmash', ['convert',
                                                        '--force',
                                                        test_sig],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == 0
        # check that --force overwrote the file
        assert int(timestamp) != int(os.path.getmtime(test_sig + ".json"))
        # check that the file can be read (as JSON)
        with open(test_sig + ".json") as fh:
            sig = signature.signature_json.load_signatures_json(fh)
