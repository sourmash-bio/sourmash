from __future__ import print_function, unicode_literals

from glob import glob
import os

from . import sourmash_tst_utils as utils
from sourmash_lib import signature
from sourmash_lib.sbt import SBT, GraphFactory, Leaf
from sourmash_lib.sbtmh import SigLeaf, search_minhashes


SIG_FILES = [os.path.join('demo', f) for f in (
  "SRR2060939_1.sig", "SRR2060939_2.sig", "SRR2241509_1.sig",
  "SRR2255622_1.sig", "SRR453566_1.sig", "SRR453569_1.sig", "SRR453570_1.sig")
]


def test_simple():
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory)

    leaf1 = Leaf("a", factory())
    leaf1.data.count('AAAAA')
    leaf1.data.count('AAAAT')
    leaf1.data.count('AAAAC')

    leaf2 = Leaf("b", factory())
    leaf2.data.count('AAAAA')
    leaf2.data.count('AAAAT')
    leaf2.data.count('AAAAG')

    leaf3 = Leaf("c", factory())
    leaf3.data.count('AAAAA')
    leaf3.data.count('AAAAT')
    leaf3.data.count('CAAAA')

    leaf4 = Leaf("d", factory())
    leaf4.data.count('AAAAA')
    leaf4.data.count('CAAAA')
    leaf4.data.count('GAAAA')

    leaf5 = Leaf("e", factory())
    leaf5.data.count('AAAAA')
    leaf5.data.count('AAAAT')
    leaf5.data.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def search_kmer(obj, seq):
        return obj.data.get(seq)

    leaves = [leaf1, leaf2, leaf3, leaf4, leaf5 ]
    kmers = [ "AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA" ]

    def search_kmer_in_list(kmer):
        x = []
        for l in leaves:
            if l.data.get(kmer):
                x.append(l)

        return set(x)

    for kmer in kmers:
        assert set(root.find(search_kmer, kmer)) == search_kmer_in_list(kmer)

    print('-----')
    print([ x.metadata for x in root.find(search_kmer, "AAAAA") ])
    print([ x.metadata for x in root.find(search_kmer, "AAAAT") ])
    print([ x.metadata for x in root.find(search_kmer, "AAAAG") ])
    print([ x.metadata for x in root.find(search_kmer, "CAAAA") ])
    print([ x.metadata for x in root.find(search_kmer, "GAAAA") ])

def test_longer_search():
    ksize = 5
    factory = GraphFactory(ksize, 100, 3)
    root = SBT(factory)

    leaf1 = Leaf("a", factory())
    leaf1.data.count('AAAAA')
    leaf1.data.count('AAAAT')
    leaf1.data.count('AAAAC')

    leaf2 = Leaf("b", factory())
    leaf2.data.count('AAAAA')
    leaf2.data.count('AAAAT')
    leaf2.data.count('AAAAG')

    leaf3 = Leaf("c", factory())
    leaf3.data.count('AAAAA')
    leaf3.data.count('AAAAT')
    leaf3.data.count('CAAAA')

    leaf4 = Leaf("d", factory())
    leaf4.data.count('AAAAA')
    leaf4.data.count('CAAAA')
    leaf4.data.count('GAAAA')

    leaf5 = Leaf("e", factory())
    leaf5.data.count('AAAAA')
    leaf5.data.count('AAAAT')
    leaf5.data.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.data.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * (len(seq) - ksize + 1)):
            return 1
        return 0

    try1 = [ x.metadata for x in root.find(search_transcript, "AAAAT", 1.0) ]
    assert set(try1) == set([ 'a', 'b', 'c', 'e' ]), try1 # no 'd'

    try2 = [ x.metadata for x in root.find(search_transcript, "GAAAAAT", 0.6) ]
    assert set(try2) == set([ 'a', 'b', 'c', 'd', 'e' ])

    try3 = [ x.metadata for x in root.find(search_transcript, "GAAAA", 1.0) ]
    assert set(try3) == set([ 'd', 'e' ]), try3


def test_tree_v1_load():
    tree_v1 = SBT.load(utils.get_test_data('v1.sbt.json'),
                       leaf_loader=SigLeaf.load)

    tree_v2 = SBT.load(utils.get_test_data('v2.sbt.json'),
                       leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(SIG_FILES[0])
    to_search = next(signature.load_signatures(testdata1))

    results_v1 = {str(s) for s in tree_v1.find(search_minhashes,
                                               to_search, 0.1)}
    results_v2 = {str(s) for s in tree_v2.find(search_minhashes,
                                               to_search, 0.1)}

    assert results_v1 == results_v2
    assert len(results_v1) == 4


def test_tree_save_load():
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory)

    for f in SIG_FILES:
        sig = next(signature.load_signatures(utils.get_test_data(f)))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    old_result = [str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)]
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'))
        tree = SBT.load(os.path.join(location, 'demo'),
                        leaf_loader=SigLeaf.load)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = [str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)]
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_binary_nary_tree():
    factory = GraphFactory(31, 1e5, 4)
    trees = {}
    trees[2] = SBT(factory)
    trees[5] = SBT(factory, d=5)
    trees[10] = SBT(factory, d=10)

    n_leaves = 0
    for f in SIG_FILES:
        sig = next(signature.load_signatures(utils.get_test_data(f)))
        leaf = SigLeaf(os.path.basename(f), sig)
        for tree in trees.values():
            tree.add_node(leaf)
        to_search = leaf
        n_leaves += 1

    assert all([len(t.leaves()) == n_leaves for t in trees.values()])

    results = {}
    print('*' * 60)
    print("{}:".format(to_search.metadata))
    for d, tree in trees.items():
        results[d] = [str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)]
    print(*results[2], sep='\n')

    assert set(results[2]) == set(results[5])
    assert set(results[5]) == set(results[10])
