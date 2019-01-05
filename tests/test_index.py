from __future__ import print_function, unicode_literals

import sourmash
from sourmash.index import LinearIndex
from sourmash_lib.sbt import SBT, GraphFactory, Leaf
from . import sourmash_tst_utils as utils


def test_simple_index(n_children):
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1 = Leaf("a", factory())
    leaf1.data.count("AAAAA")
    leaf1.data.count("AAAAT")
    leaf1.data.count("AAAAC")

    leaf2 = Leaf("b", factory())
    leaf2.data.count("AAAAA")
    leaf2.data.count("AAAAT")
    leaf2.data.count("AAAAG")

    leaf3 = Leaf("c", factory())
    leaf3.data.count("AAAAA")
    leaf3.data.count("AAAAT")
    leaf3.data.count("CAAAA")

    leaf4 = Leaf("d", factory())
    leaf4.data.count("AAAAA")
    leaf4.data.count("CAAAA")
    leaf4.data.count("GAAAA")

    leaf5 = Leaf("e", factory())
    leaf5.data.count("AAAAA")
    leaf5.data.count("AAAAT")
    leaf5.data.count("GAAAA")

    root.insert(leaf1)
    root.insert(leaf2)
    root.insert(leaf3)
    root.insert(leaf4)
    root.insert(leaf5)

    def search_kmer(obj, seq):
        return obj.data.get(seq)

    kmers = ["AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA"]

    linear = LinearIndex()
    linear.insert(leaf1)
    linear.insert(leaf2)
    linear.insert(leaf3)
    linear.insert(leaf4)
    linear.insert(leaf5)

    for kmer in kmers:
        assert set(root.find(search_kmer, kmer)) == set(linear.find(search_kmer, kmer))

    print("-----")
    print([x.metadata for x in root.find(search_kmer, "AAAAA")])
    print([x.metadata for x in root.find(search_kmer, "AAAAT")])
    print([x.metadata for x in root.find(search_kmer, "AAAAG")])
    print([x.metadata for x in root.find(search_kmer, "CAAAA")])
    print([x.metadata for x in root.find(search_kmer, "GAAAA")])


def test_linear_index_search():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s.name for s in sr])
    assert len(sr) == 1
    assert sr[0].match_sig == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s.name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x.similarity)
    assert sr[0].match_sig == ss47
    assert sr[1].match_sig == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s.name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x.similarity)
    assert sr[0].match_sig == ss63
    assert sr[1].match_sig == ss47

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s.name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x.similarity)
    assert sr[0].match_sig == ss63
