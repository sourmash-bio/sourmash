from __future__ import print_function, unicode_literals

import os

import pytest

from sourmash import load_one_signature
from sourmash.sbt import SBT, GraphFactory, Leaf, Node
from sourmash.sbtmh import (SigLeaf, search_minhashes,
                            search_minhashes_containment)
from sourmash.sbt_storage import (FSStorage, TarStorage,
                                  RedisStorage, IPFSStorage)

from . import sourmash_tst_utils as utils


def test_simple(n_children):
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

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


def test_longer_search(n_children):
    ksize = 5
    factory = GraphFactory(ksize, 100, 3)
    root = SBT(factory, d=n_children)

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

    tree_cur = SBT.load(utils.get_test_data('v4.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    results_v1 = {str(s) for s in tree_v1.find(search_minhashes_containment,
                                               to_search, 0.1)}
    results_cur = {str(s) for s in tree_cur.find(search_minhashes_containment,
                                                 to_search, 0.1)}

    assert results_v1 == results_cur
    assert len(results_v1) == 4


def test_tree_v2_load():
    tree_v2 = SBT.load(utils.get_test_data('v2.sbt.json'),
                       leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v4.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    results_v2 = {str(s) for s in tree_v2.find(search_minhashes_containment,
                                               to_search, 0.1)}
    results_cur = {str(s) for s in tree_cur.find(search_minhashes_containment,
                                                 to_search, 0.1)}

    assert results_v2 == results_cur
    assert len(results_v2) == 4


def test_tree_v3_load():
    tree_v2 = SBT.load(utils.get_test_data('v3.sbt.json'),
                       leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v4.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    results_v2 = {str(s) for s in tree_v2.find(search_minhashes_containment,
                                               to_search, 0.1)}
    results_cur = {str(s) for s in tree_cur.find(search_minhashes_containment,
                                                 to_search, 0.1)}

    assert results_v2 == results_cur
    assert len(results_v2) == 4


def test_tree_v5_load():
    tree_v2 = SBT.load(utils.get_test_data('v5.sbt.json'),
                       leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v4.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    results_v2 = {str(s) for s in tree_v2.find(search_minhashes_containment,
                                               to_search, 0.1)}
    results_cur = {str(s) for s in tree_cur.find(search_minhashes_containment,
                                                 to_search, 0.1)}

    assert results_v2 == results_cur
    assert len(results_v2) == 4


def test_tree_save_load(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    old_result = {str(s) for s in tree.find(search_minhashes,
                                            to_search.data, 0.1)}
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'))
        tree = SBT.load(os.path.join(location, 'demo'),
                        leaf_loader=SigLeaf.load)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_tree_save_load_v5(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    old_result = {str(s) for s in tree.find(search_minhashes,
                                            to_search.data, 0.1)}
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree._save_v5(os.path.join(location, 'demo'))
        tree = SBT.load(os.path.join(location, 'demo'),
                        leaf_loader=SigLeaf.load)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_search_minhashes():
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)

    to_search = next(iter(tree.leaves()))

    # this fails if 'search_minhashes' is calc containment and not similarity.
    results = tree.find(search_minhashes, to_search.data, 0.08)
    for leaf in results:
        assert to_search.data.similarity(leaf.data) >= 0.08

    print(results)


def test_binary_nary_tree():
    factory = GraphFactory(31, 1e5, 4)
    trees = {}
    trees[2] = SBT(factory)
    trees[5] = SBT(factory, d=5)
    trees[10] = SBT(factory, d=10)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        for tree in trees.values():
            tree.add_node(leaf)
        to_search = leaf
        n_leaves += 1

    assert all([len(list(t.leaves())) == n_leaves for t in trees.values()])

    results = {}
    print('*' * 60)
    print("{}:".format(to_search.metadata))
    for d, tree in trees.items():
        results[d] = {str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)}
    print(*results[2], sep='\n')

    assert results[2] == results[5]
    assert results[5] == results[10]


def test_sbt_combine(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)
    tree_1 = SBT(factory, d=n_children)
    tree_2 = SBT(factory, d=n_children)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        if n_leaves < 4:
            tree_1.add_node(leaf)
        else:
            tree_2.add_node(leaf)
        n_leaves += 1

    tree_1.combine(tree_2)

    t1_leaves = {str(l) for l in tree_1.leaves()}
    t_leaves = {str(l) for l in tree.leaves()}

    assert len(t1_leaves) == n_leaves
    assert len(t_leaves) == len(t1_leaves)
    assert t1_leaves == t_leaves

    to_search = load_one_signature(utils.get_test_data(utils.SIG_FILES[0]))
    t1_result = {str(s) for s in tree_1.find(search_minhashes,
                                             to_search, 0.1)}
    tree_result = {str(s) for s in tree.find(search_minhashes,
                                             to_search, 0.1)}
    assert t1_result == tree_result

    # TODO: save and load both trees

    # check if adding a new node will use the next empty position
    next_empty = 0
    for n, (d, _) in enumerate(tree_1):
        if n != d:
            next_empty = n
            break
    if not next_empty:
        next_empty = n + 1

    tree_1.add_node(SigLeaf(to_search.name(), to_search))
    assert tree_1.next_node == next_empty


def test_sbt_fsstorage():
    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        old_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*old_result, sep='\n')

        with FSStorage(location, '.fstree') as storage:
            tree.save(os.path.join(location, 'tree'), storage=storage)

        tree = SBT.load(os.path.join(location, 'tree'), leaf_loader=SigLeaf.load)
        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*new_result, sep='\n')

        assert old_result == new_result

        assert os.path.exists(os.path.join(location, tree.storage.subdir))
        assert os.path.exists(os.path.join(location, '.fstree'))


def test_sbt_tarstorage():
    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        old_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*old_result, sep='\n')

        with TarStorage(os.path.join(location, 'tree.tar.gz')) as storage:
            tree.save(os.path.join(location, 'tree'), storage=storage)

        with TarStorage(os.path.join(location, 'tree.tar.gz')) as storage:
            tree = SBT.load(os.path.join(location, 'tree'),
                            leaf_loader=SigLeaf.load,
                            storage=storage)

            print('*' * 60)
            print("{}:".format(to_search.metadata))
            new_result = {str(s) for s in tree.find(search_minhashes,
                                                    to_search.data, 0.1)}
            print(*new_result, sep='\n')

            assert old_result == new_result


def test_sbt_ipfsstorage():
    ipfshttpclient = pytest.importorskip('ipfshttpclient')

    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        old_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*old_result, sep='\n')

        try:
            with IPFSStorage() as storage:
                tree.save(os.path.join(location, 'tree'), storage=storage)
        except ipfshttpclient.exceptions.ConnectionError:
            pytest.xfail("ipfs not installed/functioning probably")

        with IPFSStorage() as storage:
            tree = SBT.load(os.path.join(location, 'tree'),
                            leaf_loader=SigLeaf.load,
                            storage=storage)

            print('*' * 60)
            print("{}:".format(to_search.metadata))
            new_result = {str(s) for s in tree.find(search_minhashes,
                                                    to_search.data, 0.1)}
            print(*new_result, sep='\n')

            assert old_result == new_result


def test_sbt_redisstorage():
    redis = pytest.importorskip('redis')
    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        old_result = {str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)}
        print(*old_result, sep='\n')

        try:
            with RedisStorage() as storage:
                tree.save(os.path.join(location, 'tree'), storage=storage)
        except redis.exceptions.ConnectionError:
            pytest.xfail("Couldn't connect to redis server")

        with RedisStorage() as storage:
            tree = SBT.load(os.path.join(location, 'tree'),
                            leaf_loader=SigLeaf.load,
                            storage=storage)

            print('*' * 60)
            print("{}:".format(to_search.metadata))
            new_result = {str(s) for s in tree.find(search_minhashes,
                                                    to_search.data, 0.1)}
            print(*new_result, sep='\n')

            assert old_result == new_result


def test_tree_repair():
    tree_repair = SBT.load(utils.get_test_data('leaves.sbt.json'),
                           leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v3.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    results_repair = {str(s) for s in tree_repair.find(search_minhashes,
                                                       to_search, 0.1)}
    results_cur = {str(s) for s in tree_cur.find(search_minhashes,
                                                 to_search, 0.1)}

    assert results_repair == results_cur
    assert len(results_repair) == 2


def test_tree_repair_insert():
    tree_repair = SBT.load(utils.get_test_data('leaves.sbt.json'),
                           leaf_loader=SigLeaf.load)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree_repair.add_node(leaf)

    for pos, node in tree_repair:
        # Every parent of a node must be an internal node (and not a leaf),
        # except for node 0 (the root), whose parent is None.
        if pos != 0:
            assert isinstance(tree_repair.parent(pos).node, Node)

        # Leaf nodes can't have children
        if isinstance(node, Leaf):
            assert all(c.node is None for c in tree_repair.children(pos))


def test_save_sparseness(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    old_result = {str(s) for s in tree.find(search_minhashes,
                                            to_search.data, 0.1)}
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'), sparseness=1.0)
        tree_loaded = SBT.load(os.path.join(location, 'demo'),
                               leaf_loader=SigLeaf.load)
        assert all(not isinstance(n, Node) for _, n in tree_loaded)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s) for s in tree_loaded.find(search_minhashes,
                                                       to_search.data, 0.1)}
        print(*new_result, sep='\n')

        assert old_result == new_result

        for pos, node in tree_loaded:
            # Every parent of a node must be an internal node (and not a leaf),
            # except for node 0 (the root), whose parent is None.
            if pos != 0:
                assert isinstance(tree_loaded.parent(pos).node, Node)

            # Leaf nodes can't have children
            if isinstance(node, Leaf):
                assert all(c.node is None for c in tree_loaded.children(pos))


def test_sbt_as_index_signatures():
    # test 'signatures' method from Index base class.
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'))
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'))

    tree.insert(sig47)
    tree.insert(sig63)

    xx = list(tree.signatures())
    assert len(xx) == 2

    assert sig47 in xx
    assert sig63 in xx
