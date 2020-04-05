from __future__ import print_function, unicode_literals
import pytest

import sourmash
from sourmash import sourmash_args
from . import sourmash_tst_utils as utils


def make_args_selector(**kw):
    class AttrDict(dict):
        def __init__(self, *args, **kwargs):
            super(AttrDict, self).__init__(*args, **kwargs)
            self.__dict__ = self

    defaults = dict(dna=False, dayhoff=False, hp=False, protein=False,
                    ksize=None)
    defaults.update(kw)
    args = AttrDict(**defaults)
    return args


def test_search_dbl2_1():
    # test basic arg/query selection with --dna and ksize=31.
    dbl2 = sourmash_args.SearchDBLoader2(require_scaled=True)
    query_params = sourmash_args.SignatureParams([31], ['DNA'], {}, {1000})

    args_selectors = make_args_selector(moltype='DNA', ksize=31)

    dbl2.query_params = query_params
    dbl2.is_query_loaded = True

    dbl2.parse_args_selectors(args_selectors)

    assert dbl2.check_query_against_arg_selectors()
    assert query_params.ksizes == { 31 }
    assert query_params.moltypes == { 'DNA' }


def test_search_dbl2_2():
    # test failed arg/query selection with --dna and ksize=33
    # against a DNA/k=31 query.
    dbl2 = sourmash_args.SearchDBLoader2(require_scaled=True)
    query_params = sourmash_args.SignatureParams([31], ['DNA'], {}, {1000})

    # k31 != k33
    args_selectors = make_args_selector(moltype='DNA', ksize=33)

    dbl2.query_params = query_params
    dbl2.is_query_loaded = True

    dbl2.parse_args_selectors(args_selectors)

    assert not dbl2.check_query_against_arg_selectors()


def test_search_dbl2_3():
    # test failed arg/query selection with --dna and ksize=31,
    # against a protein/k=31 query.
    dbl2 = sourmash_args.SearchDBLoader2(require_scaled=True)
    query_params = sourmash_args.SignatureParams([31], ['protein'], {}, {1000})

    # protein != DNA
    args_selectors = make_args_selector(dna=True, ksize=31)

    dbl2.query_params = query_params
    dbl2.is_query_loaded = True

    dbl2.parse_args_selectors(args_selectors)

    assert not dbl2.check_query_against_arg_selectors()


@utils.in_tempdir
def test_load_dbs_and_sigs(c):
    from sourmash.sourmash_args import load_dbs_and_sigs, SearchDatabaseLoader

    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('compute', '-k', '31', testdata1)
    c.run_sourmash('compute', '-k', '31', testdata2)

    c.run_sourmash('index', 'k31-1.sbt.json', 'short.fa.sig')
    c.run_sourmash('index', 'k31-2.sbt.json', 'short2.fa.sig')

    loader = SearchDatabaseLoader([c.output('k31-1.sbt.json'),
                                   c.output('k31-2.sbt.json')],
                                  False, False)

    loader.load_all()

    assert loader.ksize == 31
    assert loader.moltype == 'DNA'
    assert loader.scaled == False


@utils.in_tempdir
def test_load_dbs_and_sigs_incompatible_ksize(c):
    from sourmash.sourmash_args import load_dbs_and_sigs, SearchDatabaseLoader

    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    c.run_sourmash('compute', '-k', '31', testdata1)
    c.run_sourmash('compute', '-k', '33', testdata2)

    c.run_sourmash('index', 'k31.sbt.json', 'short.fa.sig')
    c.run_sourmash('index', 'k33.sbt.json', 'short2.fa.sig')

    loader = SearchDatabaseLoader([c.output('k31.sbt.json'),
                                   c.output('k33.sbt.json')],
                                  False, False)

    loader.load_all()

    with pytest.raises(ValueError):
        loader.ksize
    
    assert loader.moltype == 'DNA'
    assert loader.scaled == False
