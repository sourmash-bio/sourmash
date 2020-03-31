from sourmash import MinHash, SourmashSignature
from sourmash.sbt import GraphFactory
from sourmash.sbtmh import LocalizedSBT


def test_localized_add_node(track_abundance):
    factory = GraphFactory(5, 100, 3)
    root = LocalizedSBT(factory, track_abundance=track_abundance)

    n_hashes = 5
    a = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    a.add("AAAAA")
    a.add("AAAAA")  # add k-mer twice for track abundance
    a.add("AAAAA")  # add k-mer thrice for track abundance
    a.add('AAAAT')
    a.add('AAAAC')
    sig_a = SourmashSignature(a, name='a')

    b = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    b.add("AAAAA")
    b.add("AAAAC")
    b.add('AAAAT')
    b.add('TTTTT')
    b.add('AAAAG')  # Different k-mer from above, but most similar
    sig_b = SourmashSignature(b, name='b')

    c = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    c.add("AAAAA")
    c.add("AAAAA")  # add k-mer twice for track abundance
    c.add("AAAAA")  # add k-mer thrice for track abundance
    c.add('AAAAT')
    c.add('AAAAC')
    sig_c = SourmashSignature(c, name='c')

    d = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    d.add("CAAAA")
    d.add("CAAAA")
    d.add('TAAAA')
    d.add('GAAAA')
    d.add('CCCCC')
    sig_d = SourmashSignature(d, name='d')

    # Add "b" signature in adversarial order. Most similar to "a" but added last
    root.insert(sig_a)
    root.insert(sig_c)
    root.insert(sig_d)
    root.insert(sig_b)

    # create mapping from leaf name to node pos
    leaf_pos = {
        sig.data.name(): n
        for n, sig in
        root._leaves.items()
    }

    # Compute similarities
    from sourmash.compare import compare_all_pairs
    sigs = [sig_a, sig_b, sig_c, sig_d]
    ignore_abundance = not track_abundance
    print('\n--- ignore_abundance: {0} ---'.format(ignore_abundance))
    print(compare_all_pairs(sigs, ignore_abundance=ignore_abundance))

    # verify most similar leaves are sharing same parent node
    if track_abundance:
        assert root.parent(leaf_pos["a"]) == root.parent(leaf_pos["c"])
        assert root.parent(leaf_pos["b"]) == root.parent(leaf_pos["d"])
    else:
        assert root.parent(leaf_pos["a"]) == root.parent(leaf_pos["b"])
        assert root.parent(leaf_pos["c"]) == root.parent(leaf_pos["d"])
