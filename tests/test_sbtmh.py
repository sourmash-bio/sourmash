from sourmash import MinHash, SourmashSignature
from sourmash.sbt import GraphFactory
from sourmash.sbtmh import SigLeaf, LocalizedSBT


def test_localized_add_node(n_children, track_abundance):
    factory = GraphFactory(5, 100, 3)
    root = LocalizedSBT(factory, d=n_children)

    a = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    a.add("AAAAA")
    a.add("AAAAA")  # add k-mer twice for track abundance
    a.add('AAAAT')
    a.add('AAAAC')
    sig1 = SourmashSignature(a, name='a')
    leaf1 = SigLeaf(sig1.name(), sig1)

    b = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    b.add("AAAAA")
    b.add("AAAAC")
    b.add('AAAAT')
    b.add('AAAAG')  # Different k-mer from above, but most similar
    sig2 = SourmashSignature(b, name='b')
    leaf2 = SigLeaf(sig2.name(), sig2)

    c = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    c.add("AAAAA")
    c.add("AAAAA")  # add k-mer twice for track abundance
    c.add('AAAAT')
    c.add('GAAAA')
    sig3 = SourmashSignature(c, name='c')
    leaf3 = SigLeaf(sig3.name(), sig3)

    d = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    d.add("CAAAA")
    d.add("CAAAA")
    d.add('TAAAA')
    d.add('GAAAA')
    sig4 = SourmashSignature(d, name='d')
    leaf4 = SigLeaf(sig4.name(), sig2)

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)

    # create mapping from leaf name to node pos
    leaf_nodes = {
        sig.name: n
        for n, sig in
        root._leaves.items()
    }

    # verify most similar leaves are sharing same parent node
    assert root.parent(leaf_nodes["a"]) == root.parent(leaf_nodes["b"])
    assert root.parent(leaf_nodes["c"]) == root.parent(leaf_nodes["d"])

    raise ValueError
