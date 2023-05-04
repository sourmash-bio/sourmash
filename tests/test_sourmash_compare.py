"""
(More) sourmash compare & plot tests
"""
import os
import glob
import csv
import pprint

import sourmash_tst_utils as utils

from sourmash import sourmash_args


def test_sourmash_compare_labels_to(runtmp):
    # test compare --labels-to
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)
    assert len(testsigs) == 4

    labels_out = runtmp.output('labels.csv')

    runtmp.sourmash('compare', '-o', 'cmp', '-k', '7',
                    '--labels-to', labels_out,
                    *testsigs)

    assert os.path.exists(labels_out)

    with sourmash_args.FileInputCSV(labels_out) as r:
        assert set(r.fieldnames) == { 'order', 'md5', 'label', 'name',
                                      'filename', 'signature_file' }

        rows = list(r)

    assert len(rows) == 4
    print(rows)

    d = {}
    for row in rows:
        d[row['md5']] = row['order'], row['label'], row['name'], row['filename'], row['signature_file']

    pprint.pprint(d)
    
    r1 = d['76e45c8afb20b29a7fa022b1562b0971']
    order, label, name, filename, location = r1
    assert order == '1'
    assert label == 'genome-s12'
    assert name == 'genome-s12'
    assert filename.endswith('genome-s12.fa.gz')
    assert location.endswith('genome-s12.fa.gz.sig')

    r2 = d['93d5d09abf399740a4506310680eb62c']
    order, label, name, filename, location = r2
    assert order == '4'
    assert label == 'genome-s10+s11'
    assert name == 'genome-s10+s11'
    assert filename == '-'
    assert location.endswith('-s10+s11.sig')


def test_sourmash_plot_labels_from(runtmp):
    # test plot --labels-from
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)
    assert len(testsigs) == 4

    labels_out = runtmp.output('labels.csv')

    runtmp.sourmash('compare', '-o', 'cmp', '-k', '7',
                    '--labels-to', labels_out,
                    *testsigs)

    assert os.path.exists(labels_out)

    runtmp.sourmash('plot', 'cmp', '--labels-from', labels_out)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert 'loading labels from CSV file' in runtmp.last_result.err
    assert '(NOTE: --labels-from implies --labels)' in runtmp.last_result.err
