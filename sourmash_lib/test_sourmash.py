from __future__ import print_function
import os
import glob
import gzip

from . import sourmash_tst_utils as utils
from . import Estimators
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass

from sourmash_lib import signature

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


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
                                           ['compute', '--name', 'foo',
                                            testdata1, '-o', 'foo.sig'],
                                           in_directory=location)

        sigfile = os.path.join(location, 'foo.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name() == 'foo'


def test_do_sourmash_compute_name_fail_no_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
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
        ksizes = set([ x.estimator.ksize for x in siglist ])
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
            ksizes = set([ x.estimator.ksize for x in siglist ])
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
            ksizes = set([ x.estimator.ksize for x in siglist ])
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
            ksizes = set([ x.estimator.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes

            moltype = set([ x.estimator.is_molecule_type('protein')
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
        ksizes = set([ x.estimator.ksize for x in siglist ])
        assert 21 in ksizes
        assert 31 in ksizes


def test_do_sourmash_compute_with_cardinality():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--with-cardinality',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        cards = [ x.estimator.hll.estimate_cardinality() for x in siglist ]
        assert len(cards) == 2
        assert set(cards) == set([ 966, 986 ])


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

        cards = [ x.estimator.max_hash for x in siglist ]
        assert len(cards) == 2
        assert set(cards) == set([ 2**64/100 ])


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
        assert round(sig2_aa.similarity(sig1_trans), 3) == 0.273
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


def test_sourmash_search():
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
        assert '1 matches' in out
        assert '0.958' in out


def test_sourmash_search_gzip():
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
        assert '1 matches' in out
        assert '0.958' in out


def test_sourmash_search_2():
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
        assert '2 matches' in out
        assert '0.958' in out


def test_mash_csv_to_sig():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa.msh.dump')
        testdata2 = utils.get_test_data('short.fa')

        status, out, err = utils.runscript('sourmash', ['import_csv',
                                                        testdata1,
                                                        '-o', 'xxx.sig'],
                                           in_directory=location)
        
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '31',
                                            'short.fa.sig', 'xxx.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches; showing 3:' in out


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


def test_sourmash_compare_with_abundance_1():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = Estimators(ksize=5, n=5, protein=False,
                                     track_abundance=True)
        E2 = Estimators(ksize=5, n=5, protein=False,
                                     track_abundance=True)

        E1.mh.add_sequence('ATGGA')
        E2.mh.add_sequence('ATGGA')

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


def test_sourmash_compare_with_abundance_2():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = Estimators(ksize=5, n=5, protein=False,
                        track_abundance=True)
        E2 = Estimators(ksize=5, n=5, protein=False,
                        track_abundance=True)

        E1.mh.add_sequence('ATGGA')

        E1.mh.add_sequence('ATGGA')
        E2.mh.add_sequence('ATGGA')

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


def test_sourmash_compare_with_abundance_3():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = Estimators(ksize=5, n=5, protein=False,
                        track_abundance=True)
        E2 = Estimators(ksize=5, n=5, protein=False,
                        track_abundance=True)

        E1.mh.add_sequence('ATGGA')
        E1.mh.add_sequence('GGACA')

        E1.mh.add_sequence('ATGGA')
        E2.mh.add_sequence('ATGGA')

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
