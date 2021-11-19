"""
Tests for the 'sourmash' command line.
"""
import argparse
import os
import gzip
import shutil
import screed
import glob
import json
import csv
import pytest
import sys
import zipfile

import sourmash_tst_utils as utils

import sourmash
from sourmash import MinHash
from sourmash.sbt import SBT, Node
from sourmash.sbtmh import SigLeaf, load_sbt_index
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass

from sourmash import signature
from sourmash import VERSION
from sourmash.sourmash_args import load_pathlist_from_file
from sourmash_tst_utils import SourmashCommandFailed


def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_badcmd():
    status, out, err = utils.runscript('sourmash', ['foobarbaz'], fail_ok=True)
    assert status != 0                    # bad arg!
    assert "cmd: invalid choice" in err


def test_run_sourmash_subcmd_help():
    status, out, err = utils.runscript('sourmash', ['sbt'], fail_ok=True)
    print(out)
    print(err)
    assert status != 0               # should fail

    assert "invalid choice:" in err
    assert "'sbt' (choose from" in err

    # should not have printed a Traceback
    assert any("Traceback" not in o for o in (err, out))


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


def test_load_pathlist_from_file_does_not_exist():
    from sourmash.sourmash_args import load_pathlist_from_file
    with pytest.raises(ValueError) as e:
        load_pathlist_from_file("")
    assert "file '' does not exist" in str(e.value)


@utils.in_tempdir
def test_load_pathlist_from_file_empty(c):
    file_list = c.output("file_list")
    with open(file_list, "w") as fp:
        fp.write("")
    with pytest.raises(ValueError) as e:
        load_pathlist_from_file(file_list)
    assert "pathlist is empty" in str(e.value)


@utils.in_tempdir
def test_load_pathlist_from_file_badly_formatted(c):
    file_list = c.output("file_list")
    with open(file_list, "w") as fp:
        fp.write("{'a':1}")
    with pytest.raises(ValueError) as e:
        load_pathlist_from_file(file_list)
    assert "file '{'a':1}' inside the pathlist does not exist" in str(e.value)


@utils.in_tempdir
def test_load_pathlist_from_file_badly_formatted_2(c):
    file_list = c.output("file_list")
    sig1 = utils.get_test_data('compare/genome-s10.fa.gz.sig')
    with open(file_list, "w") as fp:
        fp.write(sig1 + "\n")
        fp.write("{'a':1}")
    with pytest.raises(ValueError) as e:
        load_pathlist_from_file(file_list)
    assert "file '{'a':1}' inside the pathlist does not exist" in str(e.value)


@utils.in_tempdir
def test_load_pathlist_from_file_duplicate(c):
    file_list = c.output("file_list")
    sig1 = utils.get_test_data('compare/genome-s10.fa.gz.sig')
    with open(file_list, "w") as fp:
        fp.write(sig1 + "\n")
        fp.write(sig1 + "\n")
    check = load_pathlist_from_file(file_list)
    print (check)
    assert len(check) == 1


@utils.in_tempdir
def test_do_serial_compare(c):
    # try doing a compare serial
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--dna', *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile)

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

        sigs = []
        for fn in testsigs:
            sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                    select_moltype='dna'))
    assert (cmp_out == cmp_calc).all()


@utils.in_tempdir
def test_do_compare_parallel(c):
    # try doing a compare parallel
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--dna',
                   "--processes", "2", *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile)

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

        sigs = []
        for fn in testsigs:
            sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                    select_moltype='dna'))
    assert (cmp_out == cmp_calc).all()


@utils.in_tempdir
def test_do_serial_compare_with_from_file(c):
    # try doing a compare serial
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    file_list = c.output('file.list')
    with open(file_list, 'wt') as fp:
        print("\n".join(testsigs), file=fp)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--dna',
                   '--from-file', file_list)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile)

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

        sigs = []
        for fn in testsigs:
            sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                    select_moltype='dna'))

    assert numpy.array_equal(numpy.sort(cmp_out.flat), numpy.sort(cmp_calc.flat))


@utils.in_tempdir
def test_do_basic_compare_using_rna_arg(c):
    # try doing a basic compare using --rna instead of --dna
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--rna', *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile)

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

    assert (cmp_out == cmp_calc).all()


def test_do_basic_compare_using_nucleotide_arg(runtmp):
    # try doing a basic compare using --nucleotide instead of --dna/--rna
    c=runtmp
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--nucleotide', *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile)

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

    assert (cmp_out == cmp_calc).all()


@utils.in_tempdir
def test_do_compare_quiet(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('compare', 'short.fa.sig',
                   'short2.fa.sig', '--csv', 'xxx', '-q')
    assert not c.last_result.out
    assert not c.last_result.err

@utils.in_tempdir
def test_do_traverse_directory_compare(c):
    c.run_sourmash('compare', '-k 21',
                   '--dna', utils.get_test_data('compare'))
    print(c.last_result.out)
    assert 'genome-s10.fa.gz' in c.last_result.out
    assert 'genome-s11.fa.gz' in c.last_result.out


@utils.in_tempdir
def test_do_traverse_directory_compare_force(c):
    sig1 = utils.get_test_data('compare/genome-s10.fa.gz.sig')
    sig2 = utils.get_test_data('compare/genome-s11.fa.gz.sig')
    newdir = c.output('newdir')
    os.mkdir(newdir)

    shutil.copyfile(sig1, os.path.join(newdir, 'sig1'))
    shutil.copyfile(sig2, os.path.join(newdir, 'sig2'))

    c.run_sourmash('compare', '-k 21',
                   '--dna', newdir, '-f')
    print(c.last_result.out)
    assert 'genome-s10.fa.gz' in c.last_result.out
    assert 'genome-s11.fa.gz' in c.last_result.out


@utils.in_tempdir
def test_do_compare_output_csv(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1, testdata2)
    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '--csv', 'xxx')

    with open(c.output('xxx')) as fp:
        r = iter(csv.reader(fp))
        row = next(r)
        print(row)
        row = next(r)
        print(row)
        assert float(row[0]) == 1.0
        assert float(row[1]) == 0.93
        row = next(r)
        assert float(row[0]) == 0.93
        assert float(row[1]) == 1.0

        # exactly three lines
        with pytest.raises(StopIteration) as e:
            next(r)


@utils.in_tempdir
def test_do_compare_downsample(c):
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=31,scaled=200', testdata1)

    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=31,scaled=100', testdata2)

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '--csv', 'xxx')

    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert 'downsampling to scaled value of 200' in c.last_result.err
    with open(c.output('xxx')) as fp:
        lines = fp.readlines()
        assert len(lines) == 3
        assert lines[1].startswith('1.0,0.6666')
        assert lines[2].startswith('0.6666')


@utils.in_tempdir
def test_do_compare_output_multiple_k(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', testdata1)
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '--csv', 'xxx',
                       fail_ok=True)

    print(c.last_result.status, c.last_result.out, c.last_result.err)

    assert c.last_result.status == -1
    assert 'multiple k-mer sizes loaded; please specify one' in c.last_result.err
    assert '(saw k-mer sizes 21, 31)' in c.last_result.err


@utils.in_tempdir
def test_do_compare_output_multiple_moltype(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=21,num=500', testdata1)
    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '--csv', 'xxx',
                       fail_ok=True)

    assert c.last_result.status == -1
    print(c.last_result.err)
    assert 'multiple molecule types loaded;' in c.last_result.err


@utils.in_tempdir
def test_do_compare_dayhoff(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', '--dayhoff', testdata1)
    assert c.last_result.status == 0

    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', '--dayhoff', testdata2)
    assert c.last_result.status == 0

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig',
                   '--dayhoff', '--csv', 'xxx')
    true_out = '''[1.   0.94]
[0.94 1.  ]
min similarity in matrix: 0.940'''.splitlines()
    for line in c.last_result.out:
        cleaned_line = line.split('...')[-1].strip()
        cleaned_line in true_out
    assert c.last_result.status == 0


@utils.in_tempdir
def test_do_compare_hp(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', '--hp', testdata1)
    assert c.last_result.status == 0

    c.run_sourmash('sketch', 'translate', '-p', 'k=21,num=500', '--hp', testdata2)
    assert c.last_result.status == 0

    c.run_sourmash('compare', 'short.fa.sig',
                   'short2.fa.sig', '--hp',  '--csv', 'xxx')
    true_out = '''[1.   0.94]
[0.94 1.  ]
min similarity in matrix: 0.940'''.splitlines()
    for line in c.last_result.out:
        cleaned_line = line.split('...')[-1].strip()
        cleaned_line in true_out
    assert c.last_result.status == 0


@utils.in_tempdir
def test_compare_containment(c):
    import numpy

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    c.run_sourmash('compare', '--containment', '-k', '31',
                   '--csv', 'output.csv', *testdata_sigs)

    # load the matrix output of compare --containment
    with open(c.output('output.csv'), 'rt') as fp:
        r = iter(csv.reader(fp))
        headers = next(r)

        mat = numpy.zeros((len(headers), len(headers)))
        for i, row in enumerate(r):
            for j, val in enumerate(row):
                mat[i][j] = float(val)

        print(mat)

    # load in all the input signatures
    idx_to_sig = dict()
    for idx, filename in enumerate(testdata_sigs):
        ss = sourmash.load_one_signature(filename, ksize=31)
        idx_to_sig[idx] = ss

    # check explicit containment against output of compare
    for i in range(len(idx_to_sig)):
        ss_i = idx_to_sig[i]
        for j in range(len(idx_to_sig)):
            ss_j = idx_to_sig[j]
            containment = ss_j.contained_by(ss_i)
            containment = round(containment, 3)
            mat_val = round(mat[i][j], 3)

            assert containment == mat_val, (i, j)


@utils.in_tempdir
def test_compare_max_containment(c):
    import numpy

    testdata_glob = utils.get_test_data('scaled/*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    c.run_sourmash('compare', '--max-containment', '-k', '31',
                   '--csv', 'output.csv', *testdata_sigs)

    # load the matrix output of compare --containment
    with open(c.output('output.csv'), 'rt') as fp:
        r = iter(csv.reader(fp))
        headers = next(r)

        mat = numpy.zeros((len(headers), len(headers)))
        for i, row in enumerate(r):
            for j, val in enumerate(row):
                mat[i][j] = float(val)

        print(mat)

    # load in all the input signatures
    idx_to_sig = dict()
    for idx, filename in enumerate(testdata_sigs):
        ss = sourmash.load_one_signature(filename, ksize=31)
        idx_to_sig[idx] = ss

    # check explicit containment against output of compare
    for i in range(len(idx_to_sig)):
        ss_i = idx_to_sig[i]
        for j in range(len(idx_to_sig)):
            ss_j = idx_to_sig[j]
            containment = ss_j.max_containment(ss_i)
            containment = round(containment, 3)
            mat_val = round(mat[i][j], 3)

            assert containment == mat_val, (i, j)


@utils.in_tempdir
def test_compare_max_containment_and_containment(c):
    testdata_glob = utils.get_test_data('scaled/*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('compare', '--max-containment', '-k', '31',
                       '--containment',
                       '--csv', 'output.csv', *testdata_sigs)

    print(c.last_result.err)
    assert "ERROR: cannot specify both --containment and --max-containment!" in c.last_result.err


@utils.in_tempdir
def test_compare_containment_abund_flatten(c):
    s47 = utils.get_test_data('track_abund/47.fa.sig')
    s63 = utils.get_test_data('track_abund/63.fa.sig')

    c.run_sourmash('compare', '--containment', '-k', '31', s47, s63)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'NOTE: --containment and --max-containment ignore signature abundances.' in \
        c.last_result.err


@utils.in_tempdir
def test_compare_containment_require_scaled(c):
    s47 = utils.get_test_data('num/47.fa.sig')
    s63 = utils.get_test_data('num/63.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('compare', '--containment', '-k', '31', s47, s63,
                       fail_ok=True)

    assert 'must use scaled signatures with --containment and --max-containment' in \
        c.last_result.err
    assert c.last_result.status != 0


@utils.in_tempdir
def test_do_plot_comparison(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '-o', 'cmp')

    c.run_sourmash('plot', 'cmp')

    assert os.path.exists(c.output("cmp.dendro.png"))
    assert os.path.exists(c.output("cmp.matrix.png"))


@utils.in_tempdir
def test_do_plot_comparison_2(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '-o', 'cmp')

    c.run_sourmash('plot', 'cmp', '--pdf')
    assert os.path.exists(c.output("cmp.dendro.pdf"))
    assert os.path.exists(c.output("cmp.matrix.pdf"))


@utils.in_tempdir
def test_do_plot_comparison_3(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '-o', 'cmp')

    c.run_sourmash('plot', 'cmp', '--labels')

    assert os.path.exists(c.output("cmp.dendro.png"))
    assert os.path.exists(c.output("cmp.matrix.png"))


@utils.in_tempdir
def test_do_plot_comparison_4_output_dir(c):
    output_dir = c.output('xyz_test')

    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('compare', 'short.fa.sig', 'short2.fa.sig', '-o', 'cmp')

    c.run_sourmash('plot', 'cmp', '--labels', '--output-dir', output_dir)

    assert os.path.exists(os.path.join(output_dir, "cmp.dendro.png"))
    assert os.path.exists(os.path.join(output_dir, "cmp.matrix.png"))


@utils.in_tempdir
def test_do_plot_comparison_5_force(c):
    import numpy
    D = numpy.zeros([2, 2])
    D[0, 0] = 5
    with open(c.output('cmp'), 'wb') as fp:
        numpy.save(fp, D)

    with open(c.output('cmp.labels.txt'), 'wt') as fp:
        fp.write("a\nb\n")

    c.run_sourmash('plot', 'cmp', '--labels', '-f')
    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert c.last_result.status == 0


@utils.in_tempdir
def test_do_plot_comparison_4_fail_not_distance(c):
    import numpy
    D = numpy.zeros([2, 2])
    D[0, 0] = 5
    with open(c.output('cmp'), 'wb') as fp:
        numpy.save(fp, D)

    with open(c.output('cmp.labels.txt'), 'wt') as fp:
        fp.write("a\nb\n")

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('plot', 'cmp', '--labels', fail_ok=True)

    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert c.last_result.status != 0


def test_plot_override_labeltext(runtmp):
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        runtmp.run_sourmash('compare', testdata1, testdata2, testdata3, testdata4, '-o', 'cmp', '-k', '21', '--dna')

        with open(runtmp.output('new.labels.txt'), 'wt') as fp:
            fp.write('a\nb\nc\nd\n')

        runtmp.sourmash('plot', 'cmp', '--labeltext', 'new.labels.txt')

        print(runtmp.last_result.out)

        assert 'loading labels from new.labels.txt' in runtmp.last_result.err

        expected = """\
0\ta
1\tb
2\tc
3\td"""
        assert expected in runtmp.last_result.out


def test_plot_override_labeltext_fail(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    runtmp.sourmash('compare', testdata1, testdata2, testdata3, testdata4, '-o', 'cmp', '-k', '21', '--dna')

    with open(runtmp.output('new.labels.txt'), 'wt') as fp:
        fp.write('a\nb\nc\n')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('plot', 'cmp', '--labeltext', 'new.labels.txt')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0
    assert 'loading labels from new.labels.txt' in runtmp.last_result.err
    assert '3 labels != matrix size, exiting' in runtmp.last_result.err


@utils.in_tempdir
def test_plot_reordered_labels_csv(c):
    files = utils.get_test_data('demo/*.sig')
    files = glob.glob(files)
    files.sort()
    assert len(files) == 7

    c.run_sourmash('compare', '-o', 'cmp', *files)
    c.run_sourmash('plot', 'cmp', '--csv', 'neworder.csv')

    with open(c.output('neworder.csv'), 'rt') as fp:
        out_mat = fp.readlines()

    # turns out to be hard to guarantee output order, so... just make sure
    # matrix labels are in different order than inputs!

    header = out_mat[0].strip().split(',')

    files = [ os.path.basename(x)[:-4] + '.fastq.gz' for x in files ]

    print(files)
    print(header)

    assert set(files) == set(header) # same file names...
    assert files != header           # ...different order.


def test_plot_subsample_1(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    runtmp.sourmash('compare', testdata1, testdata2, testdata3, testdata4, '-o', 'cmp', '-k', '21', '--dna')

    runtmp.sourmash('plot', 'cmp', '--subsample', '3')

    print(runtmp.last_result.out)

    expected = """\
0\tgenome-s10+s11
1\tgenome-s12
2\tgenome-s10"""
    assert expected in runtmp.last_result.out


def test_plot_subsample_2(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    runtmp.sourmash('compare', testdata1, testdata2, testdata3, testdata4, '-o', 'cmp', '-k', '21', '--dna')

    runtmp.sourmash('plot', 'cmp', '--subsample', '3', '--subsample-seed=2')

    print(runtmp.last_result.out)
    expected = """\
0\tgenome-s12
1\tgenome-s10+s11
2\tgenome-s11"""
    assert expected in runtmp.last_result.out


@utils.in_tempdir
def test_search_query_sig_does_not_exist(c):
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('search', 'short2.fa.sig', 'short.fa.sig', fail_ok=True)

    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert c.last_result.status == -1
    assert "Cannot open query file 'short2.fa.sig'" in c.last_result.err
    assert len(c.last_result.err.split('\n\r')) < 5


@utils.in_tempdir
def test_search_subject_sig_does_not_exist(c):
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('search', 'short.fa.sig', 'short2.fa.sig', fail_ok=True)

    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert c.last_result.status == -1
    assert "Error while reading signatures from 'short2.fa.sig'" in c.last_result.err


@utils.in_tempdir
def test_search_second_subject_sig_does_not_exist(c):
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('search', 'short.fa.sig', 'short.fa.sig',
                       'short2.fa.sig', fail_ok=True)

    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert c.last_result.status == -1
    assert "Error while reading signatures from 'short2.fa.sig'." in c.last_result.err


@utils.in_tempdir
def test_search(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('search', 'short.fa.sig', 'short2.fa.sig')
    print(c.last_result.status, c.last_result.out, c.last_result.err)
    assert '1 matches' in c.last_result.out
    assert '93.0%' in c.last_result.out


def test_search_ignore_abundance(runtmp):
    # note: uses num signatures.
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    runtmp.sourmash('sketch', 'dna', '-p','k=31,num=500,abund', testdata1, testdata2)

    # Make sure there's different percent matches when using or
    # not using abundance
    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig')
    out1 = runtmp.last_result.out
    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '81.5%' in runtmp.last_result.out

    runtmp.sourmash('search', '--ignore-abundance', 'short.fa.sig', 'short2.fa.sig')
    out2 = runtmp.last_result.out
    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '93.0%' in runtmp.last_result.out

    # Make sure results are different!
    assert out1 != out2


@utils.in_tempdir
def test_search_csv(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1, testdata2)

    c.run_sourmash('search', 'short.fa.sig', 'short2.fa.sig', '-o', 'xxx.csv')
    print(c.last_result.status, c.last_result.out, c.last_result.err)

    csv_file = c.output('xxx.csv')

    with open(csv_file) as fp:
        reader = csv.DictReader(fp)
        row = next(reader)
        print(row)
        assert float(row['similarity']) == 0.93
        assert row['filename'].endswith('short2.fa.sig')
        assert row['md5'] == '914591cd1130aa915fe0c0c63db8f19d'
        assert row['query_filename'].endswith('short.fa')
        assert row['query_name'] == ''
        assert row['query_md5'] == 'e26a306d'


@utils.in_tempdir
def test_search_lca_db(c):
    # can we do a 'sourmash search' on an LCA database?
    query = utils.get_test_data('47.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    c.run_sourmash('search', query, lca_db)
    print(c)
    assert 'NC_009665.1 Shewanella baltica OS185, complete genome' in str(c)


def test_search_query_db_md5(runtmp):
    # pull a search query out of a database with an md5sum
    db = utils.get_test_data('prot/protein.sbt.zip')
    runtmp.run_sourmash('search', db, db, '--md5', '16869d2c8a1')

    assert '100.0%       GCA_001593925' in str(runtmp)


def test_gather_query_db_md5(runtmp, linear_gather, prefetch_gather):
    # pull a search query out of a database with an md5sum
    db = utils.get_test_data('prot/protein.sbt.zip')
    runtmp.run_sourmash('gather', db, db, '--md5', '16869d2c8a1',
                        linear_gather, prefetch_gather)

    assert '340.9 kbp    100.0%  100.0%    GCA_001593925' in str(runtmp)


def test_gather_query_db_md5_ambiguous(runtmp, linear_gather, prefetch_gather):
    c = runtmp
    # what if we give an ambiguous md5 prefix?
    db = utils.get_test_data('prot/protein.sbt.zip')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('gather', db, db, '--md5', '1', linear_gather,
                       prefetch_gather)

    err = c.last_result.err
    assert "Error! Multiple signatures start with md5 '1'" in err


def test_gather_lca_db(runtmp, linear_gather, prefetch_gather):
    # can we do a 'sourmash gather' on an LCA database?
    query = utils.get_test_data('47+63.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    runtmp.sourmash('gather', query, lca_db, linear_gather, prefetch_gather)
    print(runtmp)
    out = runtmp.last_result.out

    assert 'NC_009665.1 Shewanella baltica OS185' in out
    assert 'WARNING: final scaled was 10000, vs query scaled of 1000' in out


def test_gather_csv_output_filename_bug(runtmp, linear_gather, prefetch_gather):
    c = runtmp

    # check a bug where the database filename in the output CSV was incorrect
    query = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db_1 = utils.get_test_data('lca/delmont-1.lca.json')
    lca_db_2 = utils.get_test_data('lca/delmont-2.lca.json')

    c.run_sourmash('gather', query, lca_db_1, lca_db_2, '-o', 'out.csv',
                   linear_gather, prefetch_gather)
    with open(c.output('out.csv'), 'rt') as fp:
        r = csv.DictReader(fp)
        row = next(r)
        assert row['filename'] == lca_db_1


@utils.in_tempdir
def test_compare_no_such_file(c):
    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('compare', 'nosuchfile.sig')

    assert "Error while reading signatures from 'nosuchfile.sig'." in c.last_result.err


@utils.in_tempdir
def test_compare_no_such_file_force(c):
    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('compare', 'nosuchfile.sig', '-f')

    print(c.last_result.err)
    assert "Error while reading signatures from 'nosuchfile.sig'."


@utils.in_tempdir
def test_compare_no_matching_sigs(c):
    query = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.last_result.status, c.last_result.out, c.last_result.err = \
            c.run_sourmash('compare', '-k', '100', query, fail_ok=True)

    print(c.last_result.out)
    print(c.last_result.err)
    assert c.last_result.status
    assert 'warning: no signatures loaded at given ksize/molecule type' in c.last_result.err
    assert 'no signatures found! exiting.' in c.last_result.err


def test_compare_deduce_molecule(runtmp):
    # deduce DNA vs protein from query, if it is unique
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=10,num=500', testdata1,testdata2)

    runtmp.sourmash('compare', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert 'min similarity in matrix: 0.91' in runtmp.last_result.out


def test_compare_choose_molecule_dna(runtmp):
    # choose molecule type
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('compute', '-k', '30', '--dna', '--protein', testdata1, testdata2)

    runtmp.sourmash('compare', '--dna', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert 'min similarity in matrix: 0.938' in runtmp.last_result.out


def test_compare_choose_molecule_protein(runtmp):
    # choose molecule type
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('compute', '-k', '30', '--dna', '--protein', testdata1, testdata2)

    runtmp.sourmash('compare', '--protein', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert 'min similarity in matrix: 0.91' in runtmp.last_result.out


def test_compare_no_choose_molecule_fail(runtmp):
    # choose molecule type
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=30,num=500',testdata1)

    runtmp.sourmash('sketch', 'protein', '-p', 'k=30,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('compare', 'short.fa.sig', 'short2.fa.sig')

    assert 'multiple molecule types loaded; please specify' in runtmp.last_result.err
    assert runtmp.last_result.status != 0


def test_compare_deduce_ksize(runtmp):
    # deduce ksize, if it is unique
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=29,num=500', testdata1, testdata2)

    runtmp.sourmash('compare', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert 'min similarity in matrix: 0.938' in runtmp.last_result.out


def test_search_deduce_molecule(runtmp):
    # deduce DNA vs protein from query, if it is unique
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=10,num=500', testdata1, testdata2)

    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '(k=10, protein)' in runtmp.last_result.err


def test_search_deduce_ksize(runtmp):
    # deduce ksize from query, if it is unique
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=23,num=500', testdata1, testdata2)

    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert 'k=23' in runtmp.last_result.err


def test_do_sourmash_index_multik_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch', 'translate', '-p', 'k=32,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig',  'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1


def test_do_sourmash_index_multimol_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', testdata1)

    runtmp.sourmash('sketch', 'translate', '-p', 'k=30,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1


def test_do_sourmash_index_multinum_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch', 'translate', '-p', 'k=31,num=1000', testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1
    assert 'trying to build an SBT with incompatible signatures.' in runtmp.last_result.err


def test_do_sourmash_index_multiscaled_fail(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=10', testdata1)

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=1', testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1
    assert 'trying to build an SBT with incompatible signatures.' in runtmp.last_result.err


@utils.in_tempdir
def test_do_sourmash_index_multiscaled_rescale(c):
    # test sourmash index --scaled
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('sketch', 'dna', '-p', 'scaled=10', testdata1)
    c.run_sourmash('sketch', 'dna', '-p', 'scaled=1', testdata2)

    c.run_sourmash('index', 'zzz',
                   'short.fa.sig',
                   'short2.fa.sig',
                   '-k', '31',
                   '--scaled', '10')

    print(c)
    assert c.last_result.status == 0


@utils.in_tempdir
def test_do_sourmash_index_multiscaled_rescale_fail(c):
    # test sourmash index --scaled with invalid rescaling (10 -> 5)
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('sketch', 'dna', '-p', 'scaled=10', testdata1)
    c.run_sourmash('sketch', 'dna', '-p', 'scaled=1', testdata2)
    # this should fail: cannot go from a scaled value of 10 to 5

    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('index', 'zzz',
                       'short.fa.sig',
                       'short2.fa.sig',
                       '-k', '31',
                       '--scaled', '5')

    print(e.value)
    assert c.last_result.status == -1
    assert 'new scaled 5 is lower than current sample scaled 10' in c.last_result.err


def test_do_sourmash_sbt_search_output(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1,testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz', '-o', 'foo')

    outfile = open(runtmp.output('foo'))
    output = outfile.read()
    print(output)
    assert 'e26a306d26512' in output
    assert '914591cd1130aa915' in output


# check against a bug in sbt search triggered by incorrect max Jaccard
# calculation.
def test_do_sourmash_sbt_search_check_bug(runtmp):
    # mins: 431
    testdata1 = utils.get_test_data('sbt-search-bug/nano.sig')

    # mins: 6264
    testdata2 = utils.get_test_data('sbt-search-bug/bacteroides.sig')

    runtmp.sourmash('index', 'zzz', testdata1, testdata2, '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', testdata1, 'zzz')

    assert '1 matches:' in runtmp.last_result.out

    tree = load_sbt_index(runtmp.output('zzz.sbt.zip'))
    assert tree._nodes[0].metadata['min_n_below'] == 431


def test_do_sourmash_sbt_search_empty_sig(runtmp):
    # mins: 431
    testdata1 = utils.get_test_data('sbt-search-bug/nano.sig')

    # mins: 0
    testdata2 = utils.get_test_data('sbt-search-bug/empty.sig')

    runtmp.sourmash('index', 'zzz', testdata1, testdata2, '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', testdata1, 'zzz')

    assert '1 matches:' in runtmp.last_result.out

    tree = load_sbt_index(runtmp.output('zzz.sbt.zip'))
    assert tree._nodes[0].metadata['min_n_below'] == 1


def test_do_sourmash_sbt_move_and_search_output(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1,testdata2)

    runtmp.sourmash('index', 'zzz.sbt.json', 'short.fa.sig', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.json'))

    print(runtmp.last_result.out)

    with open(runtmp.output('zzz.sbt.json')) as fp:
        d = json.load(fp)
        assert d['storage']['args']['path'] == '.sbt.zzz'

    newpath = runtmp.output('subdir')
    os.mkdir(newpath)

    # move both JSON file and subdirectory.
    shutil.move(runtmp.output('zzz.sbt.json'), newpath)
    shutil.move(runtmp.output('.sbt.zzz'), newpath)

    status, out, err = utils.runscript('sourmash',
                                        ['search', '../short.fa.sig',
                                        'zzz.sbt.json', '-o', 'foo'],
                                        in_directory=newpath)

    outfile = open(os.path.join(newpath, 'foo'))
    output = outfile.read()
    print(output)
    assert '914591cd1130aa91' in output
    assert 'e26a306d2651' in output


def test_search_deduce_ksize_and_select_appropriate(runtmp):
    # deduce ksize from query and select correct signature from DB
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=24,num=500', testdata1)

    # The DB contains signatres for multiple ksizes
    runtmp.sourmash('sketch', 'translate', '-p', 'k=23,num=500', '-p', 'k=24,num=500', testdata2)

    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert 'k=24' in runtmp.last_result.err


def test_search_deduce_ksize_not_unique(runtmp):
    # deduce ksize from query, fail because it is not unique
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=23,num=500', '-p', 'k=25,num=500', testdata1, testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1
    assert '2 signatures matching ksize' in runtmp.last_result.err


@utils.in_tempdir
def test_search_deduce_ksize_no_match(c):
    # no matching sigs in search sig list
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('sketch', 'translate', '-p', 'k=23,num=500', testdata1)
    c.run_sourmash('sketch', 'translate', '-p', 'k=25,num=500', testdata2)

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('search', 'short.fa.sig', 'short2.fa.sig')
    assert "no compatible signatures found in 'short2.fa.sig'" in str(exc.value)


def test_search_deduce_ksize_vs_user_specified(runtmp):
    # user specified ksize is not available
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=23,num=500', testdata1, testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', '-k', '24', 'short.fa.sig', 'short2.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert runtmp.last_result.status == -1
    assert '0 signatures matching ksize' in runtmp.last_result.err


def test_search_containment(runtmp):
    # search with --containment in signatures
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=1', testdata1, testdata2)

    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig', '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '95.6%' in runtmp.last_result.out


def test_search_containment_sbt(runtmp):
    # search with --containment in an SBT
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=1', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz', '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '95.6%' in runtmp.last_result.out


def test_search_containment_s10(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/genome-s10-small.fa.gz.sig')

    runtmp.sourmash('search', q1, q2, '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '16.7%' in runtmp.last_result.out


def test_search_containment_s10_no_max(run):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/genome-s10-small.fa.gz.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        run.run_sourmash('search', q1, q2, '--containment',
                       '--max-containment')

    print(run.last_result.out)
    print(run.last_result.err)
    assert "ERROR: cannot specify both --containment and --max-containment!" in run.last_result.err


def test_search_max_containment_s10_pairwise(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/genome-s10-small.fa.gz.sig')

    runtmp.sourmash('search', q1, q2,'--max-containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '100.0%' in runtmp.last_result.out


def test_search_containment_s10_siglist(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/*.sig')
    q2 = glob.glob(q2)

    runtmp.sourmash('search', q1, *q2, '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert ' 16.7%       ../genome-s10-small.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10+s11.fa.gz' in runtmp.last_result.out


def test_search_max_containment_s10_siglist(runtmp):
    # check --max-containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/*.sig')
    q2 = glob.glob(q2)

    runtmp.sourmash('search', q1, *q2, '--max-containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert '100.0%       ../genome-s10-small.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10+s11.fa.gz' in runtmp.last_result.out


def test_search_containment_s10_sbt(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert '100.0%       ../genome-s10+s11.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10.fa.gz' in runtmp.last_result.out
    assert ' 16.7%       ../genome-s10-small.fa.gz' in runtmp.last_result.out


def test_search_containment_s10_sbt_best_only(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--containment', '--best-only')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '100.0%       ' in runtmp.last_result.out # there are at least two perfect matches!

    assert runtmp.last_result.status == 0


def test_search_containment_s10_sbt_empty(runtmp):
    # check --containment for s10/s10-small at absurd scaled/empty mh
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--scaled', '1e7', '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '0 matches' in runtmp.last_result.out


def test_search_max_containment_s10_sbt(runtmp):
    # check --max-containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--max-containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert '100.0%       ../genome-s10-small.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10.fa.gz' in runtmp.last_result.out
    assert '100.0%       ../genome-s10+s11.fa.gz' in runtmp.last_result.out


def test_search_max_containment_s10_sbt_best_only(runtmp):
    # check --max-containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--max-containment', '--best-only')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert runtmp.last_result.status == 0


def test_search_max_containment_s10_sbt_empty(runtmp):
    # check --max-containment for s10/s10-small at absurd scaled/empty mh.
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.sbt.zip')

    runtmp.sourmash('search', q1, q2, '--scaled', '1e7', '--max-containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '0 matches' in runtmp.last_result.out


def test_search_containment_s10_lca(runtmp):
    # check --containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.lca.json')

    runtmp.sourmash('search', q1, q2, '--containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert '100.0%       455c2f95' in runtmp.last_result.out
    assert '100.0%       684aa226' in runtmp.last_result.out
    assert ' 16.7%       7f7835d2' in runtmp.last_result.out


def test_search_max_containment_s10_lca(runtmp):
    # check --max-containment for s10/s10-small
    q1 = utils.get_test_data('scaled/genome-s10.fa.gz.sig')
    q2 = utils.get_test_data('scaled/all.lca.json')

    runtmp.sourmash('search', q1, q2, '--max-containment')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '3 matches' in runtmp.last_result.out
    assert '100.0%       455c2f95' in runtmp.last_result.out
    assert '100.0%       684aa226' in runtmp.last_result.out
    assert '100.0%       7f7835d2' in runtmp.last_result.out


def test_search_gzip(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    data = open(runtmp.output('short.fa.sig'), 'rb').read()
    with gzip.open(runtmp.output('zzz.gz'), 'wb') as fp:
        fp.write(data)

    data = open(runtmp.output('short2.fa.sig'), 'rb').read()
    with gzip.open(runtmp.output('yyy.gz'), 'wb') as fp:
        fp.write(data)

    runtmp.sourmash('search', 'zzz.gz', 'yyy.gz')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches' in runtmp.last_result.out
    assert '93.0%' in runtmp.last_result.out


def test_search_2(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2, testdata3)

    runtmp.sourmash('search', 'short.fa.sig', 'short2.fa.sig', 'short3.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '2 matches' in runtmp.last_result.out
    assert '93.0%' in runtmp.last_result.out
    assert '89.6%' in runtmp.last_result.out


def test_search_3(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2, testdata3)

    runtmp.sourmash('search', '-n', '1', 'short.fa.sig', 'short2.fa.sig', 'short3.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '2 matches; showing first 1' in runtmp.last_result.out


def test_search_4(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2, testdata3)

    runtmp.sourmash('search', '-n', '0', 'short.fa.sig', 'short2.fa.sig', 'short3.fa.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '2 matches:' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out
    assert 'short3.fa' in runtmp.last_result.out


def test_index_check_scaled_bounds_negative(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31', '--scaled', '-5', '--dna')

    print(runtmp.last_result.err)

    assert "ERROR: scaled value must be positive" in runtmp.last_result.err


def test_index_check_scaled_bounds_less_than_minimum(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31', '--scaled', '50', '--dna')

    assert "WARNING: scaled value should be >= 100. Continuing anyway." in runtmp.last_result.err


def test_index_check_scaled_bounds_more_than_maximum(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31', '--scaled', '1e9', '--dna')

    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in runtmp.last_result.err


@utils.in_tempdir
def test_index_metagenome_fromfile(c):
    # test index --from-file
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    # construct a file list
    with open(c.output('sig.list'), 'wt') as fp:
        fp.write("\n".join(testdata_sigs))

    cmd = ['index', 'gcf_all', testdata_sigs[0], '-k', '21',
           '--from-file', c.output('sig.list')]
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    cmd = 'search {} gcf_all -k 21'.format(query_sig)
    cmd = cmd.split()
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    print(c.last_result.err)

    assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in out
    assert '12 matches; showing first 3:' in out

@utils.in_tempdir
def test_index_metagenome_fromfile_no_cmdline_sig(c):
    # test index --from-file
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    # construct a file list
    with open(c.output('sig.list'), 'wt') as fp:
        fp.write("\n".join(testdata_sigs))

    cmd = ['index', 'gcf_all', '-k', '21',
           '--from-file', c.output('sig.list')]
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    cmd = 'search {} gcf_all -k 21'.format(query_sig)
    cmd = cmd.split()
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    print(c.last_result.err)

    assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in out
    assert '12 matches; showing first 3:' in out


def test_search_metagenome(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in runtmp.last_result.out
    assert '12 matches; showing first 3:' in runtmp.last_result.out


def test_search_metagenome_traverse(runtmp):
    testdata_dir = utils.get_test_data('gather')

    query_sig = utils.get_test_data('gather/combined.sig')

    runtmp.sourmash('search', query_sig, testdata_dir, '-k', '21')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in runtmp.last_result.out
    assert '13 matches; showing first 3:' in runtmp.last_result.out


def test_search_metagenome_traverse_check_csv(runtmp):
    # this test confirms that the CSV 'filename' output for signatures loaded
    # via directory traversal properly contains the actual path to the
    # signature file from which the signature was loaded.
    testdata_dir = utils.get_test_data('gather')

    query_sig = utils.get_test_data('gather/combined.sig')
    out_csv = runtmp.output('out.csv')

    runtmp.sourmash('search', query_sig, testdata_dir, '-k', '21', '-o', out_csv)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(out_csv, 'rt') as fp:
        prefix_len = len(testdata_dir)
        r = csv.DictReader(fp)
        for row in r:
            filename = row['filename']
            assert filename.startswith(testdata_dir), filename
            # should have full path to file sig was loaded from
            assert len(filename) > prefix_len

    assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in runtmp.last_result.out
    assert '13 matches; showing first 3:' in runtmp.last_result.out


@utils.in_thisdir
def test_search_incompatible(c):
    num_sig = utils.get_test_data('num/47.fa.sig')
    scaled_sig = utils.get_test_data('47.fa.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash("search", scaled_sig, num_sig, fail_ok=True)
    assert c.last_result.status != 0
    print(c.last_result.out)
    print(c.last_result.err)

    assert "no compatible signatures found in " in c.last_result.err


@utils.in_tempdir
def test_search_traverse_incompatible(c):
    # build a directory with some signatures in it, search for compatible
    # signatures.
    searchdir = c.output('searchme')
    os.mkdir(searchdir)

    num_sig = utils.get_test_data('num/47.fa.sig')
    scaled_sig = utils.get_test_data('47.fa.sig')
    shutil.copyfile(num_sig, c.output('searchme/num.sig'))
    shutil.copyfile(scaled_sig, c.output('searchme/scaled.sig'))

    c.run_sourmash("search", scaled_sig, c.output('searchme'))
    assert '100.0%       NC_009665.1 Shewanella baltica OS185, complete genome' in c.last_result.out


def test_search_check_scaled_bounds_negative(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21', '--scaled', '-5')

    assert "ERROR: scaled value must be positive" in runtmp.last_result.err


def test_search_check_scaled_bounds_less_than_minimum(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21', '--scaled', '50')

    assert "WARNING: scaled value should be >= 100. Continuing anyway." in runtmp.last_result.err


def test_search_check_scaled_bounds_more_than_maximum(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21', '--scaled', '1e9')

    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in runtmp.last_result.err


# explanation: you cannot downsample a scaled SBT to match a scaled
# signature, so make sure that when you try such a search, it fails!
# (you *can* downsample a signature to match an SBT.)
def test_search_metagenome_downsample(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    cmd = 'search {} gcf_all -k 21 --scaled 100000'.format(query_sig)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21', '--scaled', '100000')

    assert runtmp.last_result.status == -1

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "ERROR: cannot use 'gcf_all' for this query." in runtmp.last_result.err
    assert "search scaled value 100000 is less than database scaled value of 10000" in runtmp.last_result.err


def test_search_metagenome_downsample_containment(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('search', query_sig, 'gcf_all', '-k', '21', '--scaled', '100000', '--containment')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert ' 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in runtmp.last_result.out
    assert '12 matches; showing first 3:' in runtmp.last_result.out


@utils.in_tempdir
def test_search_metagenome_downsample_index(c):
    # does same search as search_metagenome_downsample_containment but
    # rescales during indexing

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    # downscale during indexing, rather than during search.
    c.run_sourmash('index', 'gcf_all', *testdata_sigs, '-k', '21',
                   '--scaled', '100000')

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    c.run_sourmash('search', query_sig, 'gcf_all', '-k', '21',
                   '--containment')
    print(c)

    assert ' 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in str(
        c)
    assert ' 29.7%       NC_003197.2 Salmonella enterica subsp. enterica serovar T...' in str(
        c)
    assert '12 matches; showing first 3:' in str(c)


def test_search_with_picklist(runtmp):
    # test 'sourmash search' with picklists
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('search', metag_sig, *gcf_sigs, '--containment',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 3 matches to 9 distinct values" in err
    # these are the different ksizes
    assert "WARNING: 6 missing picklist values." in err

    out = runtmp.last_result.out
    print(out)
    assert "3 matches:" in out
    assert "13.1%       NC_000853.1 Thermotoga" in out
    assert "13.0%       NC_009486.1 Thermotoga" in out
    assert "12.8%       NC_011978.1 Thermotoga" in out


def test_search_with_picklist_exclude(runtmp):
    # test 'sourmash search' with picklists
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('search', metag_sig, *gcf_sigs, '--containment',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5:exclude")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 9 matches by excluding 9 distinct values" in err
    # these are the different ksizes

    out = runtmp.last_result.out
    print(out)
    assert "9 matches; showing first 3:" in out
    assert "33.2%       NC_003198.1 Salmonella" in out
    assert "33.1%       NC_003197.2 Salmonella" in out
    assert "32.2%       NC_006905.1 Salmonella" in out


def test_mash_csv_to_sig(runtmp):
    testdata1 = utils.get_test_data('short.fa.msh.dump')
    testdata2 = utils.get_test_data('short.fa')

    runtmp.sourmash('import_csv', testdata1, '-o', 'xxx.sig')

    runtmp.sourmash('sketch', 'dna', '-p','k=31,num=970',testdata2)

    runtmp.sourmash('search', '-k', '31', 'short.fa.sig', 'xxx.sig')

    print(runtmp.last_result.status, runtmp.last_result.out, runtmp.last_result.err)
    assert '1 matches:' in runtmp.last_result.out
    assert '100.0%       short.fa' in runtmp.last_result.out


def test_do_sourmash_index_bad_args(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31', '--dna', '--protein')

    print(runtmp.last_result.out, runtmp.last_result.err)
    assert 'cannot specify more than one of --dna/--rna/--nucleotide/--protein/--hp/--dayhoff' in runtmp.last_result.err
    assert runtmp.last_result.status != 0


def test_do_sourmash_sbt_search(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_wrong_ksize(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=31,num=500', '-p', 'k=51,num=500', testdata1, testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', '-k', '51', 'short.fa.sig', 'zzz')

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "ERROR: cannot use 'zzz' for this query." in runtmp.last_result.err
    assert "search ksize 51 is different from database ksize 31" in runtmp.last_result.err


def test_do_sourmash_sbt_search_multiple(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('index', 'zzz2', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz2.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz', 'zzz2')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_and_sigs(runtmp):
    # search an SBT and a signature at same time.
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz', 'short2.fa.sig')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_downsample(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,scaled=10', testdata1, testdata2)

    testdata1 = utils.get_test_data('short.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,scaled=5', '-o', 'query.sig', testdata1)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'query.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_downsample_2(runtmp):
    testdata1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
    testdata2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')

    sbtname = 'foo'

    runtmp.sourmash('index', '-k', '31', sbtname, testdata2)

    assert runtmp.last_result.status == 0

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', testdata1, sbtname, '--scaled=100000', '--threshold=0.01')

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert "ERROR: cannot use 'foo' for this query." in runtmp.last_result.err
    assert "search scaled value 100000 is less than database scaled value of 2000" in runtmp.last_result.err


@utils.in_tempdir
def test_do_sourmash_index_abund(c):
    # 'sourmash index' should flatten signatures w/track_abund.
    testdata2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')

    with open(testdata2, 'rt') as fp:
        ss = sourmash.load_one_signature(testdata2, ksize=31)
        assert ss.minhash.track_abundance == True

    sbtname = 'foo'

    c.run_sourmash('index', '-k', '31', sbtname, testdata2)

    for kk in sourmash.load_file_as_signatures(c.output(sbtname)):
        assert kk.minhash.track_abundance == False


def test_do_sourmash_index_single(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', 'short.fa.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_selectprot(runtmp):
    # index should fail when run on signatures with multiple types
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    args = ['sketch', 'dna', '-p', 'k=30,num=500',testdata1, testdata2]

    runtmp.sourmash(*args)

    args = ['index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig']

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0


def test_do_sourmash_search_multimoltype_query(runtmp):
    # 'search' should fail if multiple sigs are given as query, due to
    # multiple molecule types.
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    # first, calculate signatures with multiple molecule types
    args = ['sketch', 'translate', testdata1, testdata2,
            '-p', 'protein', '-p', 'dayhoff']
    runtmp.sourmash(*args)

    # now, index one of 'em
    args = ['index', 'zzz', 'short.fa.sig', 'short2.fa.sig', '--protein']
    runtmp.sourmash(*args)

    # output exists, yes?
    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    # now, try searching. Should raise error.
    args = ['search', 'short.fa.sig', 'zzz']
    with pytest.raises(SourmashCommandFailed) as exc:
        runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert 'need exactly one' in runtmp.last_result.err


def test_do_sourmash_index_traverse(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', '.')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))
    assert 'loaded 2 sigs; saving SBT under' in runtmp.last_result.err

    runtmp.sourmash('search', 'short.fa.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


@utils.in_tempdir
def test_do_sourmash_index_traverse_force(c):
    # test loading of files that don't end with .sig with -f
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    outdir = c.output('sigs')
    os.mkdir(outdir)
    out1 = os.path.join(outdir, 'short1')
    out2 = os.path.join(outdir, 'short2')

    c.run_sourmash('sketch','dna','-p','k=31,scaled=5', '-o', out1, testdata1)
    c.run_sourmash('sketch','dna','-p','k=31,scaled=5', '-o', out2, testdata2)

    c.run_sourmash('index', '-k', '31', 'zzz', '.', '-f')

    err = c.last_result.err
    assert os.path.exists(c.output('zzz.sbt.zip'))
    assert 'loaded 2 sigs; saving SBT under' in err

    c.run_sourmash('search', out1, 'zzz')

    out = c.last_result.out
    print(out)

    assert 'short.fa' in out
    assert 'short2.fa' in out


def test_do_sourmash_index_sparseness(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz.sbt.json', '.', '--sparseness', '1.0')

    assert os.path.exists(runtmp.output('zzz.sbt.json'))
    assert 'loaded 2 sigs; saving SBT under' in runtmp.last_result.err

    runtmp.sourmash('search', 'short.fa.sig', 'zzz.sbt.json')

    print(runtmp.last_result.out)

    assert len(glob.glob(runtmp.output('.sbt.zzz/*'))) == 3
    assert not glob.glob(runtmp.output('.sbt.zzz/*internal*'))

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_combine(runtmp):
    files = [utils.get_test_data(f) for f in utils.SIG_FILES]

    runtmp.sourmash('index', '-k', '31', 'zzz', *files)

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('sbt_combine', 'joined', 'zzz.sbt.zip', 'zzz.sbt.zip')

    assert os.path.exists(runtmp.output('joined.sbt.zip'))

    filename = os.path.splitext(os.path.basename(utils.SIG_FILES[0]))[0]

    runtmp.sourmash('search', files[0], 'zzz')

    print(runtmp.last_result.out)

    # we get notification of signature loading, too - so notify + result.
    assert runtmp.last_result.out.count(filename) == 1

    runtmp.sourmash('search', files[0], 'joined')

    print(runtmp.last_result.out)

    assert runtmp.last_result.out.count(filename) == 1


def test_do_sourmash_index_append(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1, testdata2, testdata3)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    sbt_name = runtmp.output('zzz',)
    sig_loc = runtmp.output('short3.fa.sig')

    runtmp.sourmash('search', sig_loc, sbt_name)

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out
    assert 'short3.fa' not in runtmp.last_result.out

    runtmp.sourmash('index', '-k', '31', '--append', 'zzz', 'short3.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    sbt_name = runtmp.output('zzz',)
    sig_loc = runtmp.output('short3.fa.sig')

    runtmp.sourmash('search', '--threshold', '0.95', sig_loc, sbt_name)

    print(runtmp.last_result.out)

    assert 'short.fa' not in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out
    assert 'short3.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_otherdir(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'xxx/zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('xxx/zzz.sbt.zip'))

    sbt_name = runtmp.output('xxx/zzz',)
    sig_loc = runtmp.output('short.fa.sig')

    runtmp.sourmash('search', sig_loc, sbt_name)

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out
    assert 'short2.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_scaled_vs_num_1(runtmp):
    # should not work: scaled query against num tree
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch','dna', '-p', 'scaled=1000', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    sbt_name = runtmp.output('zzz',)
    sig_loc = runtmp.output('short2.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', sig_loc, sbt_name)

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert "ERROR: cannot use '" in runtmp.last_result.err
    assert "this database was created with 'num' MinHash sketches, not 'scaled'" in runtmp.last_result.err


def test_do_sourmash_sbt_search_scaled_vs_num_2(runtmp):
    # should not work: num query against scaled tree
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch','dna', '-p', 'scaled=1000', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    sbt_name = runtmp.output('zzz',)
    sig_loc = runtmp.output('short.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', sig_loc, sbt_name)

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert "ERROR: cannot use '" in runtmp.last_result.err
    assert "this database was created with 'scaled' MinHash sketches, not 'num'" in runtmp.last_result.err


def test_do_sourmash_sbt_search_scaled_vs_num_3(runtmp):
    # should not work: scaled query against num signature
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch','dna', '-p', 'scaled=1000', testdata2)

    sig_loc = runtmp.output('short.fa.sig')
    sig_loc2 = runtmp.output('short2.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', sig_loc, sig_loc2)

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert "no compatible signatures found in " in runtmp.last_result.err


def test_do_sourmash_sbt_search_scaled_vs_num_4(runtmp):
    # should not work: num query against scaled signature
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1)

    runtmp.sourmash('sketch','dna', '-p', 'scaled=1000', testdata2)

    sig_loc = runtmp.output('short.fa.sig')
    sig_loc2 = runtmp.output('short2.fa.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('search', sig_loc2, sig_loc)

    assert runtmp.last_result.status == -1
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert "no compatible signatures found in " in runtmp.last_result.err


def test_do_sourmash_check_search_vs_actual_similarity(runtmp):
    files = [utils.get_test_data(f) for f in utils.SIG_FILES]

    runtmp.sourmash('index', '-k', '31', 'zzz', *files)

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    filename = os.path.splitext(os.path.basename(utils.SIG_FILES[0]))[0]

    runtmp.sourmash('search', files[0], 'zzz')

    assert runtmp.last_result.status == 0


def test_do_sourmash_check_sbt_filenames(runtmp):
    files = [utils.get_test_data(f) for f in utils.SIG_FILES]

    runtmp.sourmash('index', '-k', '31', 'zzz.sbt.json', *files)

    assert os.path.exists(runtmp.output('zzz.sbt.json'))

    sig_names = set()
    sig_md5s = set()
    for f in files:
        sig = signature.load_one_signature(f)
        sig_names.add(sig.name)
        sig_md5s.add(sig.md5sum())

    sbt_files = glob.glob(runtmp.output('.sbt.zzz/*'))
    assert len(sbt_files) == 14

    for f in sbt_files:
        if 'internal' in f or f.endswith('zzz.manifest.csv'):
            continue
        f = os.path.basename(f)
        assert f not in sig_names
        assert f in sig_md5s


def test_do_sourmash_sbt_search_bestonly(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'k=31,num=500', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', '--best-only', 'short.fa.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out


def test_do_sourmash_sbt_search_bestonly_scaled(runtmp):
    # as currently implemented, the query signature will be automatically
    # downsampled to match the tree.
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'scaled=1', testdata1, testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig', '--scaled', '10')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('search', '--best-only', 'short.fa.sig', 'zzz')

    print(runtmp.last_result.out)

    assert 'short.fa' in runtmp.last_result.out


def test_sbt_search_order_dependence(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz')
    testdata2 = utils.get_test_data('genome-s11.fa.gz')
    testdata3 = utils.get_test_data('genome-s12.fa.gz')
    testdata4 = utils.get_test_data('genome-s10+s11.fa.gz')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,scaled=10000', '-p', 'k=31,scaled=10000', testdata1, testdata2, testdata3, testdata4)

    runtmp.sourmash('index', '-k', '21', '134', 'genome-s10+s11.fa.gz.sig', 'genome-s11.fa.gz.sig', 'genome-s12.fa.gz.sig')

    runtmp.sourmash('search', '-k', '21', 'genome-s11.fa.gz.sig', '134', '--best-only', '-k', '21', '--dna')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert '100.0%' in runtmp.last_result.out


def test_sbt_search_order_dependence_2(runtmp):
    # *should* return the same result as test_sbt_search_order_dependence,
    # but does not due to a bug.
    testdata1 = utils.get_test_data('genome-s10.fa.gz')
    testdata2 = utils.get_test_data('genome-s11.fa.gz')
    testdata3 = utils.get_test_data('genome-s12.fa.gz')
    testdata4 = utils.get_test_data('genome-s10+s11.fa.gz')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=21,scaled=10000', '-p', 'k=31,scaled=10000', testdata1, testdata2, testdata3, testdata4)

    runtmp.sourmash('index', '-k', '21', '314', 'genome-s11.fa.gz.sig', 'genome-s10+s11.fa.gz.sig', 'genome-s12.fa.gz.sig')

    runtmp.sourmash('search', '-k', '21', 'genome-s11.fa.gz.sig', '314', '--best-only', '--dna')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert '100.0%' in runtmp.last_result.out


def test_compare_with_abundance_1(runtmp):
    # create two signatures
    E1 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)
    E2 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)

    E1.add_sequence('ATGGA')
    E2.add_sequence('ATGGA')

    s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
    s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

    signature.save_signatures([s1],
                                open(runtmp.output('e1.sig'), 'w'))
    signature.save_signatures([s2],
                                open(runtmp.output('e2.sig'), 'w'))

    runtmp.sourmash('search', 'e1.sig', 'e2.sig', '-k', '5')

    assert '100.0%' in runtmp.last_result.out


def test_compare_with_abundance_2(runtmp):
    # create two signatures
    E1 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)
    E2 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)

    E1.add_sequence('ATGGA')

    E1.add_sequence('ATGGA')
    E2.add_sequence('ATGGA')

    s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
    s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

    signature.save_signatures([s1],
                                open(runtmp.output('e1.sig'), 'w'))
    signature.save_signatures([s2],
                                open(runtmp.output('e2.sig'), 'w'))

    runtmp.sourmash('search', 'e1.sig', 'e2.sig', '-k', '5')

    assert '100.0%' in runtmp.last_result.out


def test_compare_with_abundance_3(runtmp):
    # create two signatures
    E1 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)
    E2 = MinHash(ksize=5, n=5, is_protein=False,
                    track_abundance=True)

    E1.add_sequence('ATGGA')
    E1.add_sequence('GGACA')

    E1.add_sequence('ATGGA')
    E2.add_sequence('ATGGA')

    s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
    s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

    signature.save_signatures([s1],
                                open(runtmp.output('e1.sig'), 'w'))
    signature.save_signatures([s2],
                                open(runtmp.output('e2.sig'), 'w'))

    runtmp.sourmash('search', 'e1.sig', 'e2.sig', '-k', '5')

    assert '70.5%' in runtmp.last_result.out


def test_compare_with_picklist(runtmp):
    # test 'sourmash compare' with picklists
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('compare', *gcf_sigs,
                    '-k', '21', '--picklist', f"{picklist}:md5:md5")

    err = runtmp.last_result.err
    out = runtmp.last_result.out
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "for given picklist, found 3 matches to 9 distinct values" in err
    assert "WARNING: 6 missing picklist values." in err

    assert "NC_009486.1 The..." in out
    assert "NC_000853.1 The..." in out
    assert "NC_011978.1 The..." in out


def test_compare_with_picklist_exclude(runtmp):
    # test 'sourmash compare' with picklists - exclude
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('compare', *gcf_sigs,
                    '-k', '21', '--picklist', f"{picklist}:md5:md5:exclude")

    err = runtmp.last_result.err
    out = runtmp.last_result.out
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "for given picklist, found 9 matches by excluding 9 distinct values" in err

    assert "NC_004631.1 Sal..." in out
    assert "NC_006905.1 Sal..." in out
    assert "NC_003198.1 Sal..." in out
    assert "NC_002163.1 Cam..." in out
    assert "NC_011294.1 Sal..." in out


def test_gather(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', '-o', 'foo.csv', '--threshold-bp=1', linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out


def test_gather_csv(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna','-p','scaled=10', '--name-from-first', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', '--name-from-first', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', '-o', 'foo.csv', '--threshold-bp=1', linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    csv_file = runtmp.output('foo.csv')

    with open(csv_file) as fp:
        reader = csv.DictReader(fp)
        row = next(reader)
        print(row)
        assert float(row['intersect_bp']) == 910
        assert float(row['unique_intersect_bp']) == 910
        assert float(row['remaining_bp']) == 0
        assert float(row['f_orig_query']) == 1.0
        assert float(row['f_unique_to_query']) == 1.0
        assert float(row['f_match']) == 1.0
        assert row['filename'] == 'zzz'
        assert row['name'] == 'tr1 4'
        assert row['md5'] == 'c9d5a795eeaaf58e286fb299133e1938'
        assert row['gather_result_rank'] == '0'
        assert row['query_filename'].endswith('short2.fa')
        assert row['query_name'] == 'tr1 4'
        assert row['query_md5'] == 'c9d5a795'
        assert row['query_bp'] == '910'


def test_gather_abund_x_abund(runtmp, prefetch_gather, linear_gather):
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    runtmp.sourmash('gather', sig47, sig63, linear_gather, prefetch_gather)

    assert '2.5 Mbp       49.2%   48.3%       1.0    NC_011663.1' in runtmp.last_result.out


def test_gather_multiple_sbts(runtmp, prefetch_gather, linear_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('index', 'zzz2', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', 'zzz2', '-o', 'foo.csv', '--threshold-bp=1', linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out


def test_gather_multiple_sbts_save_prefetch(runtmp, linear_gather):
    # test --save-prefetch with multiple databases
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('index', 'zzz2', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', 'zzz2', '-o', 'foo.csv', '--save-prefetch', 'out.zip', '--threshold-bp=1', linear_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out
    assert os.path.exists(runtmp.output('out.zip'))


def test_gather_multiple_sbts_save_prefetch_csv(runtmp, linear_gather):
    # test --save-prefetch-csv with multiple databases
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('index', 'zzz2', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', 'zzz2', '-o', 'foo.csv', '--save-prefetch-csv', 'prefetch.csv', '--threshold-bp=1', linear_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out
    assert os.path.exists(runtmp.output('prefetch.csv'))
    with open(runtmp.output('prefetch.csv')) as f:
        output = f.read()
        print((output,))
        assert '870,0.925531914893617,0.9666666666666667' in output


def test_gather_multiple_sbts_save_prefetch_and_prefetch_csv(runtmp, linear_gather):
    # test --save-prefetch-csv with multiple databases
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch','dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('index', 'zzz2', 'short2.fa.sig', '-k', '31')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', 'zzz2', '-o', 'foo.csv', '--save-prefetch', 'out.zip', '--save-prefetch-csv', 'prefetch.csv', '--threshold-bp=1', linear_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out
    assert os.path.exists(runtmp.output('prefetch.csv'))
    with open(runtmp.output('prefetch.csv')) as f:
        output = f.read()
        print((output,))
        assert '870,0.925531914893617,0.9666666666666667' in output
    assert os.path.exists(runtmp.output('out.zip'))


def test_gather_sbt_and_sigs(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', 'short2.fa.sig', '-o', 'foo.csv', linear_gather, prefetch_gather, '--threshold-bp=1')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out


def test_gather_file_output(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', '--threshold-bp=500', linear_gather, prefetch_gather, '-o', 'foo.out')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out
    with open(runtmp.output('foo.out')) as f:
        output = f.read()
        print((output,))
        assert '910,1.0,1.0' in output


def test_gather_f_match_orig(runtmp, linear_gather, prefetch_gather):
    import copy

    testdata_combined = utils.get_test_data('gather/combined.sig')
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    runtmp.sourmash('gather', testdata_combined, '-o', 'out.csv',
                    *testdata_sigs, linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    combined_sig = sourmash.load_one_signature(testdata_combined, ksize=21)
    remaining_mh = combined_sig.minhash.to_mutable()

    def approx_equal(a, b, n=5):
        return round(a, n) == round(b, n)

    with open(runtmp.output('out.csv'), 'rt') as fp:
        r = csv.DictReader(fp)
        for n, row in enumerate(r):
            print(n, row['f_match'], row['f_match_orig'])

            # each match is completely in the original query
            assert row['f_match_orig'] == "1.0"

            # double check -- should match 'search --containment'.
            # (this is kind of useless for a 1.0 contained_by, I guess)
            filename = row['filename']
            match = sourmash.load_one_signature(filename, ksize=21)
            assert match.contained_by(combined_sig) == 1.0

            # check other fields, too.
            f_orig_query = float(row['f_orig_query'])
            f_match_orig = float(row['f_match_orig'])
            f_match = float(row['f_match'])
            f_unique_to_query = float(row['f_unique_to_query'])

            # f_orig_query is the containment of the query by the match.
            # (note, this only works because containment is 100% in combined).
            assert approx_equal(combined_sig.contained_by(match), f_orig_query)

            # just redoing above, for completeness; this is always 1.0 for
            # this data set.
            assert approx_equal(match.contained_by(combined_sig), f_match_orig)

            # f_match is how much of the match is in the unallocated hashes
            assert approx_equal(match.minhash.contained_by(remaining_mh),
                                f_match)

            # f_unique_to_query is how much of the match is unique wrt
            # the original query.
            a = set(remaining_mh.hashes.keys())
            b = set(match.minhash.hashes.keys())
            n_intersect = len(a.intersection(b))
            f_intersect = n_intersect / float(len(combined_sig.minhash))
            assert approx_equal(f_unique_to_query, f_intersect)

            # now, subtract current match from remaining... and iterate!
            remaining_mh.remove_many(match.minhash.hashes.keys())


def test_gather_nomatch(runtmp):
    testdata_query = utils.get_test_data(
        'gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig')
    testdata_match = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

    runtmp.sourmash('gather', testdata_query, testdata_match)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 0 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 0.0% of the query' in runtmp.last_result.out


def test_gather_metagenome(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('gather', query_sig, 'gcf_all', '-k', '21', '--threshold-bp=0')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))
    assert all(('4.7 Mbp        0.5%    1.5%' in runtmp.last_result.out,
                'NC_011294.1 Salmonella enterica subsp...' in runtmp.last_result.out))


@utils.in_tempdir
def test_gather_metagenome_num_results(c):
    # set a threshold on the number of results to be reported by gather
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    cmd = 'gather {} gcf_all -k 21 --num-results 10'.format(query_sig)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    print(c.last_result.out)
    print(c.last_result.err)

    out = c.last_result.out

    assert 'found 10 matches total' in out
    assert '(truncated gather because --num-results=10)' in out
    assert 'the recovered matches hit 99.4% of the query' in out
    assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
    assert '4.3 Mbp        2.1%    7.3%    NC_006511.1 Salmonella enterica subsp' in out


def test_gather_metagenome_threshold_bp(runtmp):
    # set a threshold on the gather output
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('gather', query_sig, 'gcf_all', '-k',  '21', '--threshold-bp', '2e6')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 1 matches total' in runtmp.last_result.out
    assert 'found less than 2.0 Mbp in common. => exiting' in runtmp.last_result.err
    assert 'the recovered matches hit 33.2% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))


def test_multigather_metagenome(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('multigather', '--query', query_sig,  '--db', 'gcf_all', '-k', '21', '--threshold-bp=0')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))
    assert all(('4.7 Mbp        0.5%    1.5%' in runtmp.last_result.out,
                'NC_011294.1 Salmonella enterica subsp...' in runtmp.last_result.out))


@utils.in_tempdir
def test_multigather_check_scaled_bounds_negative(c):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    cmd = 'multigather --query {} --db gcf_all -k 21 --scaled -5 --threshold-bp=0'.format(query_sig)
    cmd = cmd.split(' ')
    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash(*cmd)

    assert "ERROR: scaled value must be positive" in str(exc.value)


@utils.in_tempdir
def test_multigather_check_scaled_bounds_less_than_minimum(c):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    cmd = 'multigather --query {} --db gcf_all -k 21 --scaled 50 --threshold-bp=0'.format(query_sig)
    cmd = cmd.split(' ')
    # Note: this is the value error that is emitted, but we want the Warning from below to be generated instead. (ValueError: new scaled 50.0 is lower than current sample scaled 10000)
    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash(*cmd)

    assert "WARNING: scaled value should be >= 100. Continuing anyway." in str(exc.value)


@utils.in_tempdir
def test_multigather_check_scaled_bounds_more_than_maximum(c):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    cmd = 'multigather --query {} --db gcf_all -k 21 --scaled 1e9 --threshold-bp=0'.format(query_sig)
    cmd = cmd.split(' ')

    c.run_sourmash(*cmd)

    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in c.last_result.err


@utils.in_tempdir
def test_multigather_metagenome_query_from_file(c):
    # test multigather --query-from-file
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    # make list w/query sig
    query_list = c.output('query.list')
    with open(query_list, 'wt') as fp:
        print(query_sig, file=fp)

    cmd = 'multigather --query-from-file {} --db gcf_all -k 21 --threshold-bp=0'.format(query_list)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'found 12 matches total' in out
    assert 'the recovered matches hit 100.0% of the query' in out
    assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
    assert all(('4.7 Mbp        0.5%    1.5%' in out,
                'NC_011294.1 Salmonella enterica subsp...' in out))


@utils.in_tempdir
def test_multigather_metagenome_query_with_sbt(c):

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    cmd = 'multigather --query gcf_all.sbt.zip --db gcf_all.sbt.zip -k 21 --threshold-bp=0'
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 12 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
    assert all(('4.7 Mbp      100.0%  100.0%'  in out,
                'NC_011080.1 Salmonella enterica subsp...' in out))
    assert all(('4.5 Mbp      100.0%  100.0%' in out,
                'NC_004631.1 Salmonella enterica subsp...' in out))
    assert all (('1.6 Mbp      100.0%  100.0%' in out,
                 'NC_002163.1 Campylobacter jejuni subs...' in out))
    assert all(('1.9 Mbp      100.0%  100.0%' in out,
                'NC_000853.1 Thermotoga maritima MSB8 ...' in out))


@utils.in_tempdir
def test_multigather_metagenome_query_with_lca(c):

    testdata_glob = utils.get_test_data('47*.fa.sig')
    testdata_sigs = glob.glob(testdata_glob)

    lca_db = utils.get_test_data('lca/47+63.lca.json')

    cmd = ['index', '47+63.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '31'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('47+63.sbt.zip'))

    cmd = 'multigather --query {} --db 47+63.sbt.zip -k 31 --threshold-bp=0'.format(lca_db)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 2 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
#    assert '5.1 Mbp      100.0%   64.9%    491c0a81'  in out
    assert '5.5 Mbp      100.0%   69.4%    491c0a81'  in out


@utils.in_tempdir
def test_multigather_metagenome_query_on_lca_db(c):

    testdata_sig1 = utils.get_test_data('47.fa.sig')
    testdata_sig2 = utils.get_test_data('63.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    cmd = 'multigather --query {} {} --db {} -k 31 --threshold-bp=0'.format(testdata_sig1, testdata_sig2, lca_db)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 2 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
    assert all(('5.1 Mbp      100.0%  100.0%'  in out,
                'NC_009665.1 Shewanella baltica OS185,...' in out))
    assert all(('5.5 Mbp      100.0%  100.0%' in out,
                'NC_011663.1 Shewanella baltica OS223,...' in out))


@utils.in_tempdir
def test_multigather_metagenome_query_with_sbt_addl_query(c):

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    another_query = utils.get_test_data('gather/GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig')

    cmd = 'multigather --query {} gcf_all.sbt.zip --db gcf_all.sbt.zip -k 21 --threshold-bp=0'.format(another_query)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 13 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
    #check for matches to some of the sbt signatures
    assert all(('4.7 Mbp      100.0%  100.0%'  in out,
                'NC_011080.1 Salmonella enterica subsp...' in out))
    assert all(('4.5 Mbp      100.0%  100.0%' in out,
                'NC_004631.1 Salmonella enterica subsp...' in out))
    assert all (('1.6 Mbp      100.0%  100.0%' in out,
                 'NC_002163.1 Campylobacter jejuni subs...' in out))
    assert all(('1.9 Mbp      100.0%  100.0%' in out,
                'NC_000853.1 Thermotoga maritima MSB8 ...' in out))

    #check additional query sig
    assert all(('4.9 Mbp      100.0%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))


@utils.in_tempdir
def test_multigather_metagenome_sbt_query_from_file_with_addl_query(c):

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    # make list w/query sbt
    query_list = c.output('query.list')
    with open(query_list, 'wt') as fp:
        print('gcf_all.sbt.zip', file=fp)

    another_query = utils.get_test_data('gather/GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig')

    cmd = 'multigather --query {} --query-from-file {} --db gcf_all.sbt.zip -k 21 --threshold-bp=0'.format(another_query, query_list)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 13 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
    #check for matches to some of the sbt signatures
    assert all(('4.7 Mbp      100.0%  100.0%'  in out,
                'NC_011080.1 Salmonella enterica subsp...' in out))
    assert all(('4.5 Mbp      100.0%  100.0%' in out,
                'NC_004631.1 Salmonella enterica subsp...' in out))
    assert all (('1.6 Mbp      100.0%  100.0%' in out,
                 'NC_002163.1 Campylobacter jejuni subs...' in out))
    assert all(('1.9 Mbp      100.0%  100.0%' in out,
                'NC_000853.1 Thermotoga maritima MSB8 ...' in out))

    #check additional query sig
    assert all(('4.9 Mbp      100.0%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))


@utils.in_tempdir
def test_multigather_metagenome_sbt_query_from_file_incorrect(c):

    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    # incorrectly query with sbt using `--query-from-file`
    cmd = 'multigather --query-from-file gcf_all.sbt.zip --db gcf_all.sbt.zip -k 21 --threshold-bp=0'
    cmd = cmd.split(' ')

    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash(*cmd)

    print(c.last_result.out)
    print(c.last_result.err)


@utils.in_tempdir
def test_multigather_metagenome_lca_query_from_file(c):
    testdata_glob = utils.get_test_data('47*.fa.sig')
    testdata_sigs = glob.glob(testdata_glob)

    lca_db = utils.get_test_data('lca/47+63.lca.json')

    cmd = ['index', '47+63.sbt.zip']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '31'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('47+63.sbt.zip'))

    # make list w/query sig
    query_list = c.output('query.list')
    with open(query_list, 'wt') as fp:
        print(lca_db, file=fp)

    cmd = 'multigather --query-from-file {} --db 47+63.sbt.zip -k 31 --threshold-bp=0'.format(query_list)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'conducted gather searches on 2 signatures' in err
    assert 'the recovered matches hit 100.0% of the query' in out
#    assert '5.1 Mbp      100.0%   64.9%    491c0a81'  in out
    assert '5.5 Mbp      100.0%   69.4%    491c0a81'  in out


@utils.in_tempdir
def test_multigather_metagenome_query_from_file_with_addl_query(c):
    # test multigather --query-from-file and --query too
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])
    c.run_sourmash(*cmd)

    assert os.path.exists(c.output('gcf_all.sbt.zip'))

    # make list w/query sig
    query_list = c.output('query.list')
    with open(query_list, 'wt') as fp:
        print(query_sig, file=fp)

    another_query = utils.get_test_data('gather/GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig')

    cmd = 'multigather --query-from-file {} --query {} --db gcf_all -k 21 --threshold-bp=0'.format(query_list, another_query)
    cmd = cmd.split(' ')
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    # first gather query
    assert 'found 12 matches total' in out
    assert 'the recovered matches hit 100.0% of the query' in out
    assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
    assert all(('4.7 Mbp        0.5%    1.5%' in out,
                'NC_011294.1 Salmonella enterica subsp...' in out))

    # second gather query
    assert '4.9 Mbp      100.0%  100.0%    NC_003198.1 Salmonella enterica subsp...' in out
    assert 'found 1 matches total;' in out
    assert 'the recovered matches hit 100.0% of the query' in out


def test_gather_metagenome_traverse(runtmp, linear_gather, prefetch_gather):
    # set up a directory $location/gather that contains
    # everything in the 'tests/test-data/gather' directory
    # *except* the query sequence, which is 'combined.sig'.
    testdata_dir = utils.get_test_data('gather')
    copy_testdata = runtmp.output('somesigs')
    shutil.copytree(testdata_dir, copy_testdata)
    os.unlink(os.path.join(copy_testdata, 'combined.sig'))

    query_sig = utils.get_test_data('gather/combined.sig')

    # now, feed in the new directory --
    runtmp.sourmash('gather', query_sig, copy_testdata, '-k', '21', '--threshold-bp=0', linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))
    assert all(('4.7 Mbp        0.5%    1.5%' in runtmp.last_result.out,
                'NC_011294.1 Salmonella enterica subsp...' in runtmp.last_result.out))


def test_gather_metagenome_traverse_check_csv(runtmp, linear_gather, prefetch_gather):
    # this test confirms that the CSV 'filename' output for signatures loaded
    # via directory traversal properly contains the actual path to the
    # signature file from which the signature was loaded.
    # set up a directory $location/gather that contains
    # everything in the 'tests/test-data/gather' directory
    # *except* the query sequence, which is 'combined.sig'.
    testdata_dir = utils.get_test_data('gather')
    copy_testdata = runtmp.output('somesigs')
    shutil.copytree(testdata_dir, copy_testdata)
    os.unlink(os.path.join(copy_testdata, 'combined.sig'))

    query_sig = utils.get_test_data('gather/combined.sig')
    out_csv = runtmp.output('out.csv')

    # now, feed in the new directory --
    runtmp.sourmash('gather', query_sig, copy_testdata, '-k', '21', '--threshold-bp=0', '-o', out_csv, linear_gather, prefetch_gather)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    with open(out_csv, 'rt') as fp:
        prefix_len = len(copy_testdata)
        r = csv.DictReader(fp)
        for row in r:
            filename = row['filename']
            assert filename.startswith(copy_testdata)
            # should have full path to file sig was loaded from
            assert len(filename) > prefix_len

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))
    assert all(('4.7 Mbp        0.5%    1.5%' in runtmp.last_result.out,
                'NC_011294.1 Salmonella enterica subsp...' in runtmp.last_result.out))


@utils.in_tempdir
def test_gather_traverse_incompatible(c):
    searchdir = c.output('searchme')
    os.mkdir(searchdir)

    num_sig = utils.get_test_data('num/47.fa.sig')
    scaled_sig = utils.get_test_data('47.fa.sig')
    shutil.copyfile(num_sig, c.output('searchme/num.sig'))
    shutil.copyfile(scaled_sig, c.output('searchme/scaled.sig'))

    c.run_sourmash("gather", scaled_sig, c.output('searchme'))
    print(c.last_result.out)
    print(c.last_result.err)
    assert "5.2 Mbp      100.0%  100.0%    NC_009665.1 Shewanella baltica OS185,..." in c.last_result.out


def test_gather_metagenome_output_unassigned(runtmp):
    testdata_glob = utils.get_test_data('gather/GCF_000195995*g')
    testdata_sigs = glob.glob(testdata_glob)[0]

    query_sig = utils.get_test_data('gather/combined.sig')

    runtmp.sourmash('gather', query_sig, testdata_sigs, '-k', '21', '--output-unassigned=unassigned.sig')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 1 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 33.2% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))

    # now examine unassigned
    testdata2_glob = utils.get_test_data('gather/GCF_000009505.1*.sig')
    testdata2_sigs = glob.glob(testdata2_glob)[0]

    runtmp.sourmash('gather', 'unassigned.sig', testdata_sigs, testdata2_sigs, '-k', '21')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert all(('1.3 Mbp       13.6%   28.2%' in runtmp.last_result.out,
                'NC_011294.1' in runtmp.last_result.out))


def test_gather_metagenome_output_unassigned_none(runtmp):
    # test what happens when there's nothing unassigned to output
    testdata_glob = utils.get_test_data('gather/GCF_*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    runtmp.sourmash('gather', query_sig, *testdata_sigs, '-k', '21', '--output-unassigned=unassigned.sig', '--threshold=0')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('4.9 Mbp       33.2%  100.0%' in runtmp.last_result.out,
                'NC_003198.1 Salmonella enterica subsp...' in runtmp.last_result.out))
    assert all(('4.5 Mbp        0.1%    0.4%' in runtmp.last_result.out,
                'NC_004631.1 Salmonella enterica subsp...' in runtmp.last_result.out))

    # now examine unassigned
    assert not os.path.exists(runtmp.output('unassigned.sig'))
    assert 'no unassigned hashes to save with --output-unassigned!' in runtmp.last_result.err


def test_gather_metagenome_output_unassigned_nomatches(runtmp, prefetch_gather, linear_gather):
    c = runtmp

    # test --output-unassigned when there are no matches
    query_sig = utils.get_test_data('2.fa.sig')
    against_sig = utils.get_test_data('47.fa.sig')

    c.run_sourmash('gather', query_sig, against_sig,
                   '--output-unassigned', 'foo.sig', linear_gather,
                   prefetch_gather)

    print(c.last_result.out)
    assert 'found 0 matches total;' in c.last_result.out

    x = sourmash.load_one_signature(query_sig, ksize=31)
    y = sourmash.load_one_signature(c.output('foo.sig'))

    assert x.minhash == y.minhash


def test_gather_metagenome_output_unassigned_nomatches_protein(runtmp, linear_gather, prefetch_gather):
    c = runtmp

    # test --output-unassigned with protein signatures
    query_sig = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    against_sig = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    c.run_sourmash('gather', query_sig, against_sig,
                   '--output-unassigned', 'foo.sig', linear_gather,
                   prefetch_gather)

    print(c.last_result.out)
    assert 'found 0 matches total;' in c.last_result.out

    c.run_sourmash('sig', 'describe', c.output('foo.sig'))
    print(c.last_result.out)

    x = sourmash.load_one_signature(query_sig, ksize=57)
    y = sourmash.load_one_signature(c.output('foo.sig'))

    assert x.minhash == y.minhash
    assert y.minhash.moltype == "protein"


def test_gather_check_scaled_bounds_negative(runtmp, prefetch_gather, linear_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('gather', query_sig, prefetch_gather, linear_gather, 'gcf_all', '-k', '21', '--scaled', '-5', '--threshold-bp', '50000')

    assert "ERROR: scaled value must be positive" in runtmp.last_result.err


def test_gather_check_scaled_bounds_less_than_minimum(runtmp, prefetch_gather, linear_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('gather', query_sig, prefetch_gather, linear_gather, 'gcf_all', '-k', '21', '--scaled', '50', '--threshold-bp', '50000')

    assert "WARNING: scaled value should be >= 100. Continuing anyway." in runtmp.last_result.err


def test_gather_check_scaled_bounds_more_than_maximum(runtmp, prefetch_gather, linear_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('gather', query_sig, prefetch_gather, linear_gather, '-k', '21', '--scaled', '1e9', '--threshold-bp', '50000')

    assert "WARNING: scaled value should be <= 1e6. Continuing anyway." in runtmp.last_result.err


def test_gather_metagenome_downsample(runtmp, prefetch_gather, linear_gather):
    # downsample w/scaled of 100,000
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('gather', query_sig, 'gcf_all', '-k', '21', '--scaled', '100000', prefetch_gather, linear_gather, '--threshold-bp', '50000')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 11 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert all(('5.2 Mbp       32.9%  100.0%' in runtmp.last_result.out,
                'NC_003198.1' in runtmp.last_result.out))
    assert all(('4.1 Mbp        0.6%    2.4%' in runtmp.last_result.out,
                '4.1 Mbp        4.4%   17.1%' in runtmp.last_result.out))


def test_gather_query_downsample(runtmp, linear_gather, prefetch_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)
    print(testdata_sigs)

    query_sig = utils.get_test_data('GCF_000006945.2-s500.sig')

    runtmp.sourmash('gather', '-k', '31', linear_gather, prefetch_gather, query_sig, *testdata_sigs)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 12 signatures' in runtmp.last_result.err
    assert all(('4.9 Mbp      100.0%  100.0%' in runtmp.last_result.out,
                'NC_003197.2' in runtmp.last_result.out))

    assert 'WARNING: final scaled was 10000, vs query scaled of 500' in runtmp.last_result.out


def test_gather_query_downsample_explicit(runtmp, linear_gather, prefetch_gather):
    # do an explicit downsampling to fix `test_gather_query_downsample`
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('GCF_000006945.2-s500.sig')

    runtmp.sourmash('gather', '-k', '31', '--scaled', '10000', linear_gather, prefetch_gather, query_sig, *testdata_sigs)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loaded 12 signatures' in runtmp.last_result.err
    assert all(('4.9 Mbp      100.0%  100.0%' in runtmp.last_result.out,
                'NC_003197.2' in runtmp.last_result.out))


def test_gather_with_picklist(runtmp, linear_gather, prefetch_gather):
    # test 'sourmash gather' with picklists
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('gather', metag_sig, *gcf_sigs, '--threshold-bp=0',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5",
                    linear_gather, prefetch_gather)

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 3 matches to 9 distinct values" in err
    # these are the different ksizes
    assert "WARNING: 6 missing picklist values." in err

    out = runtmp.last_result.out
    print(out)
    assert "found 3 matches total;" in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 Thermotoga" in out
    assert "1.9 Mbp       11.5%   89.9%    NC_011978.1 Thermotoga" in out
    assert "1.9 Mbp        6.3%   48.4%    NC_009486.1 Thermotoga" in out


def test_gather_with_picklist_exclude(runtmp, linear_gather, prefetch_gather):
    # test 'sourmash gather' with picklists - exclude
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    runtmp.sourmash('gather', metag_sig, *gcf_sigs, '--threshold-bp=0',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5:exclude",
                    linear_gather, prefetch_gather)

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 9 matches by excluding 9 distinct values" in err
    # these are the different ksizes

    out = runtmp.last_result.out
    print(out)
    assert "found 9 matches total;" in out
    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 Salmonella enterica subsp..." in out
    assert "1.6 Mbp       10.7%  100.0%    NC_002163.1 Campylobacter jejuni subs..." in out
    assert "4.8 Mbp       10.4%   31.3%    NC_003197.2 Salmonella enterica subsp..." in out
    assert "4.7 Mbp        5.2%   16.1%    NC_006905.1 Salmonella enterica subsp..." in out
    assert "4.7 Mbp        4.0%   12.6%    NC_011080.1 Salmonella enterica subsp..." in out
    assert "4.6 Mbp        2.9%    9.2%    NC_011274.1 Salmonella enterica subsp..." in out
    assert "4.3 Mbp        2.1%    7.3%    NC_006511.1 Salmonella enterica subsp..." in out
    assert "4.7 Mbp        0.5%    1.5%    NC_011294.1 Salmonella enterica subsp..." in out
    assert "4.5 Mbp        0.1%    0.4%    NC_004631.1 Salmonella enterica subsp..." in out


def test_gather_save_matches(runtmp, linear_gather, prefetch_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('gather', query_sig, 'gcf_all', '-k', '21', '--save-matches', 'save.sigs', linear_gather, prefetch_gather, '--threshold-bp', '0')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out
    assert os.path.exists(runtmp.output('save.sigs'))


def test_gather_save_matches_and_save_prefetch(runtmp, linear_gather):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    cmd = ['index', 'gcf_all']
    cmd.extend(testdata_sigs)
    cmd.extend(['-k', '21'])

    runtmp.sourmash(*cmd)

    assert os.path.exists(runtmp.output('gcf_all.sbt.zip'))

    runtmp.sourmash('gather', query_sig, 'gcf_all', '-k', '21', '--save-matches', 'save.sigs', '--save-prefetch', 'save2.sigs', linear_gather, '--threshold-bp', '0')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'found 12 matches total' in runtmp.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in runtmp.last_result.out

    matches_save = runtmp.output('save.sigs')
    prefetch_save = runtmp.output('save2.sigs')
    assert os.path.exists(matches_save)
    assert os.path.exists(prefetch_save)

    matches = list(sourmash.load_file_as_signatures(matches_save))
    prefetch = list(sourmash.load_file_as_signatures(prefetch_save))

    assert set(matches) == set(prefetch)


@utils.in_tempdir
def test_gather_error_no_sigs_traverse(c):
    # test gather applied to a directory
    query = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')

    emptydir = c.output('')

    with pytest.raises(SourmashCommandFailed) as e:
        c.run_sourmash('gather', query, emptydir)

    err = c.last_result.err
    print(err)
    assert f"Error while reading signatures from '{emptydir}'" in err
    assert not 'found 0 matches total;' in err


def test_gather_error_no_cardinality_query(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=31,num=500', testdata1, testdata2)

    testdata3 = utils.get_test_data('short3.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=31,num=500', testdata3)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('gather', 'short3.fa.sig', 'zzz', linear_gather, prefetch_gather)

    assert runtmp.last_result.status == -1
    assert "query signature needs to be created with --scaled" in runtmp.last_result.err


def test_gather_deduce_ksize(runtmp, prefetch_gather, linear_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'dna', '-p', 'k=23,scaled=10', testdata1, testdata2)

    runtmp.sourmash('sketch','dna','-p','k=23,scaled=10', '-o', 'query.fa.sig', testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', prefetch_gather, linear_gather, '--threshold-bp=1')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '0.9 kbp      100.0%  100.0%' in runtmp.last_result.out


def test_gather_deduce_moltype(runtmp, linear_gather, prefetch_gather):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    runtmp.sourmash('sketch', 'translate', '-p', 'k=10,scaled=10', testdata1,testdata2)

    runtmp.sourmash('sketch', 'translate', '-p', 'k=10,scaled=10', '-o', 'query.fa.sig',testdata2)

    runtmp.sourmash('index', 'zzz', 'short.fa.sig', 'short2.fa.sig')

    assert os.path.exists(runtmp.output('zzz.sbt.zip'))

    runtmp.sourmash('gather', 'query.fa.sig', 'zzz', linear_gather, prefetch_gather, '--threshold-bp=1')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert '1.9 kbp      100.0%  100.0%' in runtmp.last_result.out


def test_gather_abund_1_1(runtmp, linear_gather, prefetch_gather):
    c = runtmp
    #
    # make r1.fa with 2x coverage of genome s10
    # make r2.fa with 20x coverage of genome s10.
    # make r3.fa with 2x coverage of genome s11.
    #
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa

    #
    # make signature s10-s11 with r1 and r3, i.e. 1:1 abundance
    # make signature s10x10-s11 with r2 and r3, i.e. 10:1 abundance
    #
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    query = utils.get_test_data('gather-abund/reads-s10-s11.sig')
    against_list = ['genome-s10', 'genome-s11', 'genome-s12']
    against_list = ['gather-abund/' + i + '.fa.gz.sig'
                    for i in against_list]
    against_list = [utils.get_test_data(i) for i in against_list]

    status, out, err = c.run_sourmash('gather', query, *against_list,
                                      linear_gather, prefetch_gather)

    print(out)
    print(err)

    # when we project s10-s11 (r1+r3), 1:1 abundance,
    # onto s10 and s11 genomes with gather, we get:
    # * approximately 50% of each query matching (first column, p_query)
    # * approximately 80% of subject genomes contents being matched
    #   (this is due to the low coverage of 2 used to build queries)
    # * approximately 2.0 abundance (third column, avg_abund)

    assert '49.6%   78.5%       1.8    tests/test-data/genome-s10.fa.gz' in out
    assert '50.4%   80.0%       1.9    tests/test-data/genome-s11.fa.gz' in out
    assert 'genome-s12.fa.gz' not in out


def test_gather_abund_10_1(runtmp, prefetch_gather, linear_gather):
    c = runtmp
    # see comments in test_gather_abund_1_1, above.
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against_list = ['genome-s10', 'genome-s11', 'genome-s12']
    against_list = ['gather-abund/' + i + '.fa.gz.sig'
                    for i in against_list]
    against_list = [utils.get_test_data(i) for i in against_list]

    status, out, err = c.run_sourmash('gather', query, '-o', 'xxx.csv',
                                      *against_list, linear_gather,
                                      prefetch_gather)

    print(out)
    print(err)

    # when we project s10x10-s11 (r2+r3), 10:1 abundance,
    # onto s10 and s11 genomes with gather, we get:
    # * approximately 91% of s10 matching
    # * approximately 9% of s11 matching
    # * approximately 100% of the high coverage genome being matched,
    #       with only 80% of the low coverage genome
    # * approximately 2.0 abundance (third column, avg_abund) for s11,
    #   and (very) approximately 20x abundance for genome s10.

    assert '91.0%  100.0%      14.5    tests/test-data/genome-s10.fa.gz' in out
    assert '9.0%   80.0%       1.9    tests/test-data/genome-s11.fa.gz' in out
    assert 'genome-s12.fa.gz' not in out

    # check the calculations behind the above output by looking into
    # the CSV.
    with open(c.output('xxx.csv'), 'rt') as fp:
        r = csv.DictReader(fp)

        overlaps = []
        unique_overlaps = []
        f_weighted_list = []
        average_abunds = []
        remaining_bps = []

        for n, row in enumerate(r):
            assert int(row['gather_result_rank']) == n
            overlap = float(row['intersect_bp'])
            remaining_bp = float(row['remaining_bp'])
            unique_overlap = float(row['unique_intersect_bp'])
            f_weighted = float(row['f_unique_weighted'])
            average_abund = float(row['average_abund'])

            overlaps.append(overlap)
            unique_overlaps.append(unique_overlap)
            f_weighted_list.append(f_weighted)
            average_abunds.append(average_abund)
            remaining_bps.append(remaining_bp)

    weighted_calc = []
    for (overlap, average_abund) in zip(overlaps, average_abunds):
        prod = overlap*average_abund
        weighted_calc.append(prod)

    total_weighted = sum(weighted_calc)
    for prod, f_weighted in zip(weighted_calc, f_weighted_list):
        assert prod / total_weighted == f_weighted, (prod, f_weighted)

    query_sig = sourmash.load_one_signature(query)
    query_mh = query_sig.minhash

    total_bp_analyzed = sum(unique_overlaps) + remaining_bps[-1]
    total_query_bp = len(query_mh) * query_mh.scaled
    assert total_bp_analyzed == total_query_bp


def test_gather_abund_10_1_ignore_abundance(runtmp, linear_gather, prefetch_gather):
    c = runtmp
    # see comments in test_gather_abund_1_1, above.
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against_list = ['genome-s10', 'genome-s11', 'genome-s12']
    against_list = ['gather-abund/' + i + '.fa.gz.sig'
                    for i in against_list]
    against_list = [utils.get_test_data(i) for i in against_list]

    status, out, err = c.run_sourmash('gather', query,
                                      '--ignore-abundance',
                                      *against_list,
                                      linear_gather, prefetch_gather,
                                      '-o', c.output('results.csv'))


    print(out)
    print(err)

    # when we project s10x10-s11 (r2+r3), 10:1 abundance,
    # onto s10 and s11 genomes with gather --ignore-abundance, we get:
    # * approximately 50% of s10 and s11 matching (first column)
    # * approximately 100% of the high coverage genome being matched,
    #       with only 80% of the low coverage genome

    assert all(('57.2%  100.0%', 'tests/test-data/genome-s10.fa.gz' in out))
    assert all(('42.8%   80.0%', 'tests/test-data/genome-s11.fa.gz' in out))
    assert 'genome-s12.fa.gz' not in out

    with open(c.output('results.csv'), 'rt') as fp:
        r = csv.DictReader(fp)
        some_results = False
        for row in r:
            some_results = True
            assert row['average_abund'] == ''
            assert row['median_abund'] == ''
            assert row['std_abund'] == ''

        assert some_results


def test_gather_output_unassigned_with_abundance(runtmp, prefetch_gather, linear_gather):
    c = runtmp
    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against = utils.get_test_data('gather-abund/genome-s10.fa.gz.sig')

    c.run_sourmash('gather', query, against, '--output-unassigned',
                   c.output('unassigned.sig'), linear_gather, prefetch_gather)

    assert os.path.exists(c.output('unassigned.sig'))

    nomatch = sourmash.load_one_signature(c.output('unassigned.sig'))
    assert nomatch.minhash.track_abundance

    query_ss = sourmash.load_one_signature(query)
    against_ss = sourmash.load_one_signature(against)

    # unassigned should have nothing that is in the database
    nomatch_mh = nomatch.minhash
    for hashval in against_ss.minhash.hashes:
        assert hashval not in nomatch_mh.hashes

    # unassigned should have abundances from original query, if not in database
    for hashval, abund in query_ss.minhash.hashes.items():
        if hashval not in against_ss.minhash.hashes:
            assert nomatch_mh.hashes[hashval] == abund


def test_sbt_categorize(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    # all four in the current directory for categorize .
    shutil.copyfile(testdata1, runtmp.output('1.sig'))
    shutil.copyfile(testdata2, runtmp.output('2.sig'))
    shutil.copyfile(testdata3, runtmp.output('3.sig'))
    shutil.copyfile(testdata4, runtmp.output('4.sig'))

    # omit 3
    args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
    runtmp.sourmash(*args)

    # categorize all of the ones that were copied to 'location'
    args = ['categorize', 'zzz', '.',
            '--ksize', '21', '--dna', '--csv', 'out.csv']
    runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    # mash dist genome-s10.fa.gz genome-s10+s11.fa.gz
    # yields 521/1000 ==> ~0.5
    assert 'for genome-s10+s11, found: 0.50 genome-s10' in runtmp.last_result.err

    out_csv = open(runtmp.output('out.csv')).read()
    print(out_csv)
    assert './4.sig,genome-s10+s11,genome-s10,0.504' in out_csv


def test_sbt_categorize_ignore_abundance_1(runtmp):
    # --- Categorize without ignoring abundance ---
    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against_list = ['reads-s10-s11']
    against_list = ['gather-abund/' + i + '.sig'
                    for i in against_list]
    against_list = [utils.get_test_data(i) for i in against_list]

    # omit 3
    args = ['index', '--dna', '-k', '21', 'thebestdatabase'] + against_list
    runtmp.sourmash(*args)

    args = ['categorize', 'thebestdatabase',
            '--ksize', '21', '--dna', '--csv', 'out3.csv', query]
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(*args)

    assert runtmp.last_result.status != 0

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "ERROR: this search cannot be done on signatures calculated with abundance." in runtmp.last_result.err
    assert "ERROR: please specify --ignore-abundance." in runtmp.last_result.err


def test_sbt_categorize_ignore_abundance_3(runtmp):
    # --- Now categorize with ignored abundance ---
    query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
    against_list = ['reads-s10-s11']
    against_list = ['gather-abund/' + i + '.sig'
                    for i in against_list]
    against_list = [utils.get_test_data(i) for i in against_list]

    # omit 3
    args = ['index', '--dna', '-k', '21', 'thebestdatabase'] + against_list
    runtmp.sourmash(*args)

    args = ['categorize', '--ignore-abundance',
            '--ksize', '21', '--dna', '--csv', 'out4.csv',
            'thebestdatabase', query]
    runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'for 1-1, found: 0.88 1-1' in runtmp.last_result.err

    out_csv4 = open(runtmp.output('out4.csv')).read()
    assert 'reads-s10x10-s11.sig,1-1,1-1,0.87699' in out_csv4


def test_sbt_categorize_already_done(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    shutil.copyfile(testdata1, runtmp.output('1.sig'))
    shutil.copyfile(testdata2, runtmp.output('2.sig'))
    shutil.copyfile(testdata3, runtmp.output('3.sig'))
    shutil.copyfile(testdata4, runtmp.output('4.sig'))

    # omit 3
    args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
    runtmp.sourmash(*args)

    with open(runtmp.output('in.csv'), 'wt') as fp:
        fp.write('./4.sig,genome-s10.fa.gz,0.50')

    args = ['categorize', 'zzz', './2.sig', './4.sig',
            '--ksize', '21', '--dna', '--load-csv', 'in.csv']
    runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert 'for genome-s11.fa.gz, no match found'
    assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in runtmp.last_result.err


def test_sbt_categorize_already_done_traverse(runtmp):
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
    testdata4 = utils.get_test_data('genome-s10+s11.sig')

    shutil.copyfile(testdata1, runtmp.output('1.sig'))
    shutil.copyfile(testdata2, runtmp.output('2.sig'))
    shutil.copyfile(testdata3, runtmp.output('3.sig'))
    shutil.copyfile(testdata4, runtmp.output('4.sig'))

    # omit 3
    args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
    runtmp.sourmash(*args)

    with open(runtmp.output('in.csv'), 'wt') as fp:
        fp.write('./4.sig,genome-s10.fa.gz,0.50')

    args = ['categorize', 'zzz', '.',
            '--ksize', '21', '--dna', '--load-csv', 'in.csv']
    runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert 'for genome-s11.fa.gz, no match found'
    assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in runtmp.last_result.err


def test_sbt_categorize_multiple_ksizes_moltypes(runtmp):
    # 'categorize' works fine with multiple moltypes/ksizes
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
    testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')

    shutil.copyfile(testdata1, runtmp.output('1.sig'))
    shutil.copyfile(testdata2, runtmp.output('2.sig'))
    shutil.copyfile(testdata3, runtmp.output('3.sig'))

    args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
    runtmp.sourmash(*args)

    args = ['categorize', 'zzz', '.']
    runtmp.sourmash(*args)


def test_watch_check_num_bounds_negative(runtmp):
    c=runtmp
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    shutil.copyfile(testdata1, c.output('1.sig'))

    c.run_sourmash('index', '--dna', '-k', '21', 'zzz', '1.sig')

    with pytest.raises(SourmashCommandFailed) as exc:
        c.run_sourmash('watch', '--ksize', '21', '-n', '-5', '--dna', 'zzz', testdata0)

    assert "ERROR: num value must be positive" in c.last_result.err


def test_watch_check_num_bounds_less_than_minimum(runtmp):
    c=runtmp
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    shutil.copyfile(testdata1, c.output('1.sig'))

    c.run_sourmash('index', '--dna', '-k', '21', 'zzz', '1.sig')

    c.run_sourmash('watch', '--ksize', '21', '-n', '25', '--dna', 'zzz', testdata0)

    assert "WARNING: num value should be >= 50. Continuing anyway." in c.last_result.err


def test_watch_check_num_bounds_more_than_maximum(runtmp):
    c=runtmp
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    shutil.copyfile(testdata1, c.output('1.sig'))

    c.run_sourmash('index', '--dna', '-k', '21', 'zzz', '1.sig')

    c.run_sourmash('watch', '--ksize', '21', '-n', '100000', '--dna', 'zzz', testdata0)

    assert "WARNING: num value should be <= 50000. Continuing anyway." in c.last_result.err


@utils.in_tempdir
def test_watch(c):
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    shutil.copyfile(testdata1, c.output('1.sig'))

    c.run_sourmash('index', '--dna', '-k', '21', 'zzz', '1.sig')

    c.run_sourmash('watch', '--ksize', '21', '--dna', 'zzz', testdata0)

    print(c.last_result.out)
    print(c.last_result.err)
    assert 'FOUND: genome-s10, at 1.000' in c.last_result.out


@utils.in_tempdir
def test_watch_deduce_ksize(c):
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    c.run_sourmash('sketch','dna','-p','k=29,num=500', '-o', '1.sig', testdata0)

    c.run_sourmash('index', '--dna', '-k', '29', 'zzz', '1.sig')

    c.run_sourmash('watch', '--dna', 'zzz', testdata0)

    print(c.last_result.out)
    print(c.last_result.err)
    assert 'Computing signature for k=29' in c.last_result.err
    assert 'genome-s10.fa.gz, at 1.000' in c.last_result.out


def test_watch_coverage(runtmp):
    testdata0 = utils.get_test_data('genome-s10.fa.gz')
    testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
    shutil.copyfile(testdata1, runtmp.output('1.sig'))

    args = ['index', '--dna', '-k', '21', 'zzz', '1.sig']
    runtmp.sourmash(*args)

    with open(runtmp.output('query.fa'), 'wt') as fp:
        record = list(screed.open(testdata0))[0]
        for start in range(0, len(record), 100):
            fp.write('>{}\n{}\n'.format(start,
                                        record.sequence[start:start+500]))

    args = ['watch', '--ksize', '21', '--dna', 'zzz', 'query.fa']
    runtmp.sourmash(*args)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert 'FOUND: genome-s10, at 1.000' in runtmp.last_result.out


def test_storage_convert():
    import pytest

    with utils.TempDirectory() as location:
        testdata = utils.get_test_data('v2.sbt.json')
        shutil.copyfile(testdata, os.path.join(location, 'v2.sbt.json'))
        shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v2'),
                        os.path.join(location, '.sbt.v2'))
        testsbt = os.path.join(location, 'v2.sbt.json')

        original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        args = ['storage', 'convert', '-b', 'ipfs', testsbt]
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)
        if not status and "ipfs.exceptions.ConnectionError" in err:
            raise pytest.xfail('ipfs probably not running')

        ipfs = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        assert len(original) == len(ipfs)
        assert all(n1[1].name == n2[1].name
                   for (n1, n2) in zip(sorted(original), sorted(ipfs)))

        args = ['storage', 'convert',
                '-b', """'ZipStorage("{}")'""".format(
                    os.path.join(location, 'v2.sbt.zip')),
                testsbt]
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)
        tar = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        assert len(original) == len(tar)
        assert all(n1[1].name == n2[1].name
                   for (n1, n2) in zip(sorted(original), sorted(tar)))


def test_storage_convert_identity(runtmp):
    testdata = utils.get_test_data('v2.sbt.json')
    shutil.copyfile(testdata, runtmp.output('v2.sbt.json'))
    shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v2'),
                    runtmp.output('.sbt.v2'))
    testsbt = runtmp.output('v2.sbt.json')

    original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    args = ['storage', 'convert', '-b', 'fsstorage', testsbt]
    runtmp.sourmash(*args)

    identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    assert len(original) == len(identity)
    assert all(n1[1].name == n2[1].name
                for (n1, n2) in zip(sorted(original), sorted(identity)))


def test_storage_convert_fsstorage_newpath(runtmp):
    testdata = utils.get_test_data('v2.sbt.json')
    shutil.copyfile(testdata, runtmp.output('v2.sbt.json'))
    shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v2'),
                    runtmp.output('.sbt.v2'))
    testsbt = runtmp.output('v2.sbt.json')

    original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    args = ['storage', 'convert',
                        '-b', 'fsstorage({})'.format(runtmp.output('v3')),
                        testsbt]
    runtmp.sourmash(*args)

    identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    assert len(original) == len(identity)
    assert all(n1[1].name == n2[1].name
                for (n1, n2) in zip(sorted(original), sorted(identity)))


def test_migrate(runtmp):
    testdata = utils.get_test_data('v3.sbt.json')
    shutil.copyfile(testdata, runtmp.output('v3.sbt.json'))
    shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v3'),
                        runtmp.output('.sbt.v3'))
    testsbt = runtmp.output('v3.sbt.json')

    original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    runtmp.sourmash('migrate', testsbt)

    identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

    assert len(original) == len(identity)
    assert all(n1[1].name == n2[1].name
                for (n1, n2) in zip(sorted(original),
                                    sorted(identity)))

    assert "this is an old index version" not in runtmp.last_result.err
    assert all('min_n_below' in node.metadata
                for node in identity
                if isinstance(node, Node))


def test_license_cc0(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    runtmp.sourmash('sketch','translate', '-p', 'k=31', testdata1)

    sigfile = runtmp.output('short.fa.sig')
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('short.fa')

    assert sig.license == 'CC0'


def test_license_non_cc0(runtmp):
    testdata1 = utils.get_test_data('short.fa')
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('sketch', 'translate', '-p','k=31', '--license', 'GPL', testdata1)

    assert runtmp.last_result.status != 0
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert 'sourmash only supports CC0' in runtmp.last_result.err


def test_license_load_non_cc0():
    sigfile = utils.get_test_data('bad-license.sig')

    try:
        sig = next(signature.load_signatures(sigfile, do_raise=True))
    except Exception as e:
        assert "sourmash only supports CC0-licensed signatures" in str(e)


@utils.in_tempdir
def test_do_sourmash_index_zipfile(c):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    c.run_sourmash('index', '-k', '31', 'zzz.sbt.zip',
                   *testdata_sigs)

    outfile = c.output('zzz.sbt.zip')
    assert os.path.exists(outfile)

    print(c)
    assert c.last_result.status == 0
    assert 'Finished saving SBT index, available at' in c.last_result.err

    # look internally at the zip file
    with zipfile.ZipFile(outfile) as zf:
        content = zf.namelist()
        assert len(content) == 26
        assert len([c for c in content if 'internal' in c]) == 11
        assert ".sbt.zzz/" in content
        sbts = [c for c in content if c.endswith(".sbt.json")]
        assert len(sbts) == 1
        assert sbts[0] == "zzz.sbt.json"


@utils.in_tempdir
def test_do_sourmash_index_zipfile_append(c):
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)
    half_point = int(len(testdata_sigs) / 2)
    first_half = testdata_sigs[:half_point]
    second_half = testdata_sigs[half_point:]

    print(first_half)
    print(second_half)

    # should be no overlap
    assert not set(first_half).intersection(set(second_half))

    with pytest.warns(None) as record:
        c.run_sourmash('index', '-k', '31', 'zzz.sbt.zip',
                       *first_half)
    # UserWarning is raised when there are duplicated entries in the zipfile
    assert not record, record

    outfile = c.output('zzz.sbt.zip')
    assert os.path.exists(outfile)

    print(c)
    assert c.last_result.status == 0
    assert 'Finished saving SBT index, available at' in c.last_result.err

    with pytest.warns(None) as record:
        c.run_sourmash('index', "--append", '-k', '31', 'zzz.sbt.zip',
                       *second_half)
    # UserWarning is raised when there are duplicated entries in the zipfile
    print(record)
    #assert not record, record

    print(c)
    assert c.last_result.status == 0
    assert 'Finished saving SBT index, available at' in c.last_result.err

    # look internally at the zip file
    with zipfile.ZipFile(outfile) as zf:
        content = zf.namelist()
        print(content)
        assert len(content) == 26
        assert len([c for c in content if 'internal' in c]) == 11
        assert ".sbt.zzz/" in content
        sbts = [c for c in content if c.endswith(".sbt.json")]
        assert len(sbts) == 1
        assert sbts[0] == "zzz.sbt.json"


def test_index_with_picklist(runtmp):
    # test 'sourmash index' with picklists
    gcf_sig_dir = utils.get_test_data('gather/')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    output_db = runtmp.output('thermo.sbt.zip')

    runtmp.sourmash('index', output_db, gcf_sig_dir,
                    '-k', '31', '--picklist', f"{picklist}:md5:md5")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 3 matches to 9 distinct values" in err

    # these are the different ksizes
    assert "WARNING: 6 missing picklist values." in err

    # verify:
    siglist = list(sourmash.load_file_as_signatures(output_db))
    assert len(siglist) == 3
    for ss in siglist:
        assert 'Thermotoga' in ss.name


def test_index_with_picklist_exclude(runtmp):
    # test 'sourmash index' with picklists - exclude
    gcf_sig_dir = utils.get_test_data('gather/')
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')

    output_db = runtmp.output('thermo-exclude.sbt.zip')

    runtmp.sourmash('index', output_db, gcf_sig_dir,
                    '-k', '31', '--picklist', f"{picklist}:md5:md5:exclude")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 9 matches by excluding 9 distinct values" in err

    # verify:
    siglist = list(sourmash.load_file_as_signatures(output_db))
    assert len(siglist) == 9
    for ss in siglist:
        assert 'Thermotoga' not in ss.name


def test_index_matches_search_with_picklist(runtmp):
    # test 'sourmash index' with picklists
    gcf_sig_dir = utils.get_test_data('gather/')
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')
    metag_sig = utils.get_test_data('gather/combined.sig')

    output_db = runtmp.output('thermo.sbt.zip')

    runtmp.sourmash('index', output_db, gcf_sig_dir, '-k', '21')
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    # verify:
    siglist = list(sourmash.load_file_as_signatures(output_db))
    assert len(siglist) > 3     # all signatures included...

    n_thermo = 0
    for ss in siglist:
        if 'Thermotoga' in ss.name:
            n_thermo += 1

    assert n_thermo == 3

    runtmp.sourmash('search', metag_sig, output_db, '--containment',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 3 matches to 9 distinct values" in err
    # these are the different ksizes
    assert "WARNING: 6 missing picklist values." in err

    out = runtmp.last_result.out
    print(out)
    assert "3 matches:" in out
    assert "13.1%       NC_000853.1 Thermotoga" in out
    assert "13.0%       NC_009486.1 Thermotoga" in out
    assert "12.8%       NC_011978.1 Thermotoga" in out


def test_index_matches_search_with_picklist_exclude(runtmp):
    # test 'sourmash index' with picklists - exclude
    gcf_sig_dir = utils.get_test_data('gather/')
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    picklist = utils.get_test_data('gather/thermotoga-picklist.csv')
    metag_sig = utils.get_test_data('gather/combined.sig')

    output_db = runtmp.output('thermo-exclude.sbt.zip')

    runtmp.sourmash('index', output_db, gcf_sig_dir, '-k', '21')
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    # verify:
    siglist = list(sourmash.load_file_as_signatures(output_db))
    assert len(siglist) > 3     # all signatures included...

    n_thermo = 0
    for ss in siglist:
        if 'Thermotoga' in ss.name:
            n_thermo += 1

    assert n_thermo == 3

    runtmp.sourmash('search', metag_sig, output_db, '--containment',
                    '-k', '21', '--picklist', f"{picklist}:md5:md5:exclude")

    err = runtmp.last_result.err
    print(err)
    assert "for given picklist, found 10 matches by excluding 9 distinct values" in err
    ### NTP: FIX REPORTING
    assert "WARNING: -1 missing picklist values"

    out = runtmp.last_result.out
    print(out)
    assert "10 matches; showing first 3:" in out
    assert "100.0%       -" in out
    assert "33.2%       NC_003198.1 Salmonella" in out
    assert "33.1%       NC_003197.2 Salmonella" in out


def test_gather_with_prefetch_picklist(runtmp, linear_gather):
    # test 'gather' using a picklist taken from 'sourmash prefetch' output
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    prefetch_csv = runtmp.output('prefetch-out.csv')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '-k', '21', '-o', prefetch_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "total of 12 matching signatures." in err
    assert "of 1466 distinct query hashes, 1466 were found in matches above threshold." in err

    # now, do a gather with the results
    runtmp.sourmash('gather', metag_sig, *gcf_sigs, linear_gather,
                    '-k', '21', '--picklist',
                    f'{prefetch_csv}:match_md5:md5short')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out


def test_gather_with_prefetch_picklist_2_prefetch(runtmp, linear_gather):
    # test 'gather' using a picklist taken from 'sourmash prefetch' output
    # using ::prefetch
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    prefetch_csv = runtmp.output('prefetch-out.csv')

    runtmp.sourmash('prefetch', metag_sig, *gcf_sigs,
                    '-k', '21', '-o', prefetch_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "total of 12 matching signatures." in err
    assert "of 1466 distinct query hashes, 1466 were found in matches above threshold." in err

    # now, do a gather with the results
    runtmp.sourmash('gather', metag_sig, *gcf_sigs, linear_gather,
                    '-k', '21', '--picklist',
                    f'{prefetch_csv}::prefetch')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out


def test_gather_with_prefetch_picklist_3_gather(runtmp, linear_gather):
    # test 'gather' using a picklist taken from 'sourmash gather' output,
    # using ::gather.
    # (this doesn't really do anything useful, but it's an ok test :)
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    gather_csv = runtmp.output('gather-out.csv')

    runtmp.sourmash('gather', metag_sig, *gcf_sigs,
                    '-k', '21', '-o', gather_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out

    # now, do another gather with the results
    runtmp.sourmash('gather', metag_sig, *gcf_sigs, linear_gather,
                    '-k', '21', '--picklist',
                    f'{gather_csv}::gather')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out


def test_gather_with_prefetch_picklist_3_gather_badcol(runtmp):
    # test 'gather' using a picklist taken from 'sourmash gather' output,
    # using ::gather.
    # (this doesn't really do anything useful, but it's an ok test :)
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    gather_csv = runtmp.output('gather-out.csv')

    runtmp.sourmash('gather', metag_sig, *gcf_sigs,
                    '-k', '21', '-o', gather_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out

    # now, do another gather with the results, but with a bad picklist
    # parameter
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('gather', metag_sig, *gcf_sigs,
                        '-k', '21', '--picklist',
                        f'{gather_csv}:FOO:gather')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "ERROR: could not load picklist." in err
    assert "no column name allowed for coltype 'gather'" in err


def test_gather_with_prefetch_picklist_4_manifest(runtmp, linear_gather):
    # test 'gather' using a picklist taken from 'sourmash sig manifest'
    # output, using ::manifest.
    # (this doesn't really do anything useful, but it's an ok test :)
    gather_dir = utils.get_test_data('gather/')
    metag_sig = utils.get_test_data('gather/combined.sig')
    manifest_csv = runtmp.output('manifest.csv')

    runtmp.sourmash('sig', 'manifest', gather_dir, '-o', manifest_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    # now, do a gather on the manifest
    runtmp.sourmash('gather', metag_sig, gather_dir, linear_gather,
                    '-k', '21', '--picklist',
                    f'{manifest_csv}::manifest')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 1 matches total;" in out
    assert "the recovered matches hit 100.0% of the query" in out

    # the query sig itself is in there, so :shrug: that matches at 100%
    assert "14.7 Mbp     100.0%  100.0%    -" in out


def test_gather_with_prefetch_picklist_4_manifest_excl(runtmp, linear_gather):
    # test 'gather' using a picklist taken from 'sourmash sig manifest'
    # output, using ::manifest.
    # (this doesn't really do anything useful, but it's an ok test :)
    gather_dir = utils.get_test_data('gather/')
    metag_sig = utils.get_test_data('gather/combined.sig')
    manifest_csv = runtmp.output('manifest.csv')

    runtmp.sourmash('sig', 'manifest', gather_dir, '-o', manifest_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    # now, do a gather on the manifest
    runtmp.sourmash('gather', metag_sig, gather_dir, linear_gather,
                    '-k', '21', '--picklist',
                    f'{manifest_csv}::manifest:exclude')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    # excluded everything, so nothing to match!
    assert "found 0 matches total;" in out
    assert "the recovered matches hit 0.0% of the query" in out


def test_gather_with_prefetch_picklist_5_search(runtmp):
    # test 'gather' using a picklist taken from 'sourmash prefetch' output
    # using ::prefetch
    gcf_sigs = glob.glob(utils.get_test_data('gather/GCF*.sig'))
    metag_sig = utils.get_test_data('gather/combined.sig')
    search_csv = runtmp.output('search-out.csv')

    runtmp.sourmash('search', '--containment', metag_sig, *gcf_sigs,
                    '-k', '21', '-o', search_csv)

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "12 matches; showing first 3:" in out
    assert " 33.2%       NC_003198.1 Salmonella enterica subsp." in out

    # now, do a gather with the results
    runtmp.sourmash('gather', metag_sig, *gcf_sigs,
                    '-k', '21', '--picklist',
                    f'{search_csv}::search')

    err = runtmp.last_result.err
    print(err)

    out = runtmp.last_result.out
    print(out)

    assert "found 11 matches total;" in out
    assert "the recovered matches hit 99.9% of the query" in out

    assert "4.9 Mbp       33.2%  100.0%    NC_003198.1 " in out
    assert "1.9 Mbp       13.1%  100.0%    NC_000853.1 " in out


def test_gather_scaled_1(runtmp, linear_gather, prefetch_gather):
    # test gather on a sig indexed with scaled=1
    inp = utils.get_test_data('short.fa')
    outp = runtmp.output('out.sig')

    # prepare a signature with a scaled of 1
    runtmp.sourmash('sketch', 'dna', '-p', 'scaled=1,k=31', inp, '-o', outp)

    # run with a low threshold
    runtmp.sourmash('gather', outp, outp, '--threshold-bp', '0')

    print(runtmp.last_result.out)
    print('---')
    print(runtmp.last_result.err)

    assert "1.0 kbp      100.0%  100.0%" in runtmp.last_result.out
    assert "found 1 matches total;" in runtmp.last_result.out
