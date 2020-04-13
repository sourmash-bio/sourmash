from itertools import product

from . import sourmash_tst_utils as utils
from sourmash import MinHash, SourmashSignature
from sourmash.compare import compare_all_pairs
from sourmash.sbt import GraphFactory
from sourmash.sbtmh import LocalizedSBT
from sourmash import signature as sig


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
    b.add('AAAAC')  # Same k-mer from above
    sig_b = SourmashSignature(b, name='b')

    c = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    c.add("AAAAA")
    c.add("AAAAA")  # add k-mer twice for track abundance
    c.add("AAAAG")  # add k-mer thrice for track abundance
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

    # Similarity matrices for reference
    # --- track_abundance: True (ignore_abundance: False)  similarity matrix ---
    #   a          b          c          d
    # [[1.         0.7195622  0.73043556 0.        ]
    #  [0.7195622  1.         0.68749438 0.        ]
    #  [0.73043556 0.68749438 1.         0.        ]
    #  [0.         0.         0.         1.        ]]
    # --- track_abundance: False (ignore_abundance: True) similarity matrix ---
    #   a    b    c    d
    # [[1.   1.   0.75 0.  ]
    #  [1.   1.   0.75 0.  ]
    #  [0.75 0.75 1.   0.  ]
    #  [0.   0.   0.   1.  ]]

    # Add "b" signature in adversarial order. When track_abundance=False, is most
    # similar to "a" but added last
    root.insert(sig_a)
    # Tree: (track_abundance=True and track_abundance=False)
    #     0
    #   /  \
    # a: 1  None
    root.insert(sig_c)
    # Tree: (track_abundance=True and track_abundance=False)
    #     0
    #   /  \
    # a: 1  c: 2
    root.insert(sig_d)
    # Tree: (track_abundance=True)
    #          0
    #        /  \
    #      1     d: 2
    #    /   \
    # c: 3  a: 4
    # Tree: (track_abundance=False)
    #          0
    #        /  \
    #      1     d: 2
    #    /   \
    # a: 3  c: 4
    root.insert(sig_b)
    # Tree: (track_abundance=True)
    #             0
    #         /      \
    #      1           2
    #    /   \       /   \
    # a: 3  c: 4   d: 5  b: 6
    # Tree: (track_abundance=False)
    #             0
    #         /      \
    #      1           2
    #    /   \       /   \
    # a: 3  b: 4   c: 5  d: 6

    # Make sure tree construction happened properly
    assert all(node < leaf for leaf, node in product(root._leaves, root._nodes))

    # create mapping from leaf name to node pos
    leaf_pos = {
        sig.data.name(): n
        for n, sig in
        root._leaves.items()
    }

    # Verify most similar leaves are sharing same parent node
    if track_abundance:
        # Currently leaf_pos = {'a': 3, 'd': 4, 'c': 5, 'b': 6}
        # Expected leaf_pos = {'a': 3, 'd': 5, 'c': 4, 'b': 6}
        assert root.parent(leaf_pos["a"]) == root.parent(leaf_pos["c"])
        assert root.parent(leaf_pos["b"]) == root.parent(leaf_pos["d"])
    else:
        # Currently leaf_pos = {'a': 3, 'd': 4, 'c': 5, 'b': 6}
        # Expected leaf_pos = {'a': 3, 'd': 6, 'c': 5, 'b': 4}
        assert root.parent(leaf_pos["a"]) == root.parent(leaf_pos["b"])
        assert root.parent(leaf_pos["c"]) == root.parent(leaf_pos["d"])


def test_localized_sbt_more_files():
    factory = GraphFactory(5, 100, 3)
    root = LocalizedSBT(factory, track_abundance=False)

    with utils.TempDirectory() as location:
        files = [utils.get_test_data(f) for f in utils.SIG_FILES]
        signatures = []
        for filename in files:
            loaded = sig.load_signatures(filename, ksize=31)
            signatures.extend(loaded)

        compare = compare_all_pairs(signatures, ignore_abundance=True)
        print([x.name().split('.')[0] for x in signatures])
        print(compare)

        # --- track_abundance=False, ignore_abundance=True ---
        # ['SRR2060939_1', 'SRR2060939_2', 'SRR2241509_1', 'SRR2255622_1',
        #  'SRR453566_1', 'SRR453569_1', 'SRR453570_1']
        # [[1.    0.356 0.078 0.086 0.    0.    0.   ]
        #  [0.356 1.    0.072 0.078 0.    0.    0.   ]
        #  [0.078 0.072 1.    0.074 0.    0.    0.   ]
        #  [0.086 0.078 0.074 1.    0.    0.    0.   ]
        #  [0.    0.    0.    0.    1.    0.382 0.364]
        #  [0.    0.    0.    0.    0.382 1.    0.386]
        #  [0.    0.    0.    0.    0.364 0.386 1.   ]]

        # for signature in signatures:
        #     root.insert(signature)

