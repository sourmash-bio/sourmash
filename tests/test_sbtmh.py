from sourmash import MinHash, SourmashSignature
from sourmash.sbt import GraphFactory
from sourmash.sbtmh import LocalizedSBT


def test_localized_add_node(track_abundance):
    factory = GraphFactory(5, 100, 3)
    root = LocalizedSBT(factory, track_abundance=track_abundance)

    a = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    a.add("AAAAA")
    a.add("AAAAA")  # add k-mer twice for track abundance
    a.add('AAAAT')
    a.add('AAAAC')
    sig_a = SourmashSignature(a, name='a')

    b = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    b.add("AAAAA")
    b.add("AAAAC")
    b.add('AAAAT')
    b.add('AAAAG')  # Different k-mer from above, but most similar
    sig_b = SourmashSignature(b, name='b')

    c = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    c.add("AAAAA")
    c.add("AAAAA")  # add k-mer twice for track abundance
    c.add('AAAAT')
    c.add('GAAAA')
    sig_c = SourmashSignature(c, name='c')

    d = MinHash(n=1, ksize=5, track_abundance=track_abundance)
    d.add("CAAAA")
    d.add("CAAAA")
    d.add('TAAAA')
    d.add('GAAAA')
    sig_d = SourmashSignature(d, name='d')

    # Add "b" signature in adversarial order. Most similar to "a" but added last
    root.insert(sig_a)
    root.insert(sig_c)
    root.insert(sig_d)
    root.insert(sig_b)

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
