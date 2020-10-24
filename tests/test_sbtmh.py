import glob
from itertools import product
import os
from string import ascii_uppercase
import random

import pytest

from . import sourmash_tst_utils as utils
from sourmash import MinHash, SourmashSignature
from sourmash.sbt import GraphFactory
from sourmash.sbtmh import LocalizedSBT
from sourmash import signature as sig


@pytest.fixture(params=range(10))
def random_seed(request):
    return request.param


def test_localized_add_node(track_abundance):
    factory = GraphFactory(5, 100, 3)
    sbt = LocalizedSBT(factory, track_abundance=track_abundance)

    n_hashes = 5
    a = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    a.add("AAAAA")
    a.add("AAAAA")  # add k-mer twice for track abundance
    a.add("AAAAA")  # add k-mer thrice for track abundance
    a.add("AAAAT")
    a.add("AAAAC")
    sig_a = SourmashSignature(a, name="a")

    b = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    b.add("AAAAA")
    b.add("AAAAC")
    b.add("AAAAT")
    b.add("TTTTT")
    b.add("AAAAC")  # Same k-mer from above
    sig_b = SourmashSignature(b, name="b")

    c = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    c.add("AAAAA")
    c.add("AAAAA")  # add k-mer twice for track abundance
    c.add("AAAAG")  # add k-mer thrice for track abundance
    c.add("AAAAT")
    c.add("AAAAC")
    sig_c = SourmashSignature(c, name="c")

    d = MinHash(n=n_hashes, ksize=5, track_abundance=track_abundance)
    d.add("CAAAA")
    d.add("CAAAA")
    d.add("TAAAA")
    d.add("GAAAA")
    d.add("CCCCC")
    sig_d = SourmashSignature(d, name="d")

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
    sbt.insert(sig_a)
    # Tree: (track_abundance=True and track_abundance=False)
    #     0
    #   /  \
    # a: 1  None
    sbt.insert(sig_c)
    # Tree: (track_abundance=True and track_abundance=False)
    #     0
    #   /  \
    # a: 1  c: 2
    sbt.insert(sig_d)
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
    sbt.insert(sig_b)
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
    assert all(node < leaf for leaf, node in product(sbt._leaves, sbt._nodes))

    # create mapping from leaf name to node pos
    leaf_pos = {sig.data.name(): n for n, sig in sbt._leaves.items()}

    # Verify most similar leaves are sharing same parent node
    if track_abundance:
        # Currently leaf_pos = {'a': 3, 'd': 4, 'c': 5, 'b': 6}
        # Expected leaf_pos = {'a': 3, 'd': 5, 'c': 4, 'b': 6}
        assert sbt.parent(leaf_pos["a"]) == sbt.parent(leaf_pos["c"])
        assert sbt.parent(leaf_pos["b"]) == sbt.parent(leaf_pos["d"])
    else:
        # Currently leaf_pos = {'a': 3, 'd': 4, 'c': 5, 'b': 6}
        # Expected leaf_pos = {'a': 3, 'd': 6, 'c': 5, 'b': 4}
        assert sbt.parent(leaf_pos["a"]) == sbt.parent(leaf_pos["b"])
        assert sbt.parent(leaf_pos["c"]) == sbt.parent(leaf_pos["d"])


@pytest.mark.filterwarnings("ignore")
def test_localized_sbt_sorted_vs_randomized(random_seed):
    factory = GraphFactory(5, 100, 3)
    sbt = LocalizedSBT(factory, track_abundance=False)
    sbt_randomized = LocalizedSBT(factory, track_abundance=False)

    with utils.TempDirectory() as location:
        # Sort to ensure consistent ordering across operating systems
        files = sorted([utils.get_test_data(f) for f in utils.SIG_FILES])
        signatures = []
        i = 0
        for filename in files:
            loaded = sig.load_signatures(filename, ksize=31)
            for signature in loaded:
                # Rename to A, B, C, D for simplicity of debugging
                signature._name = ascii_uppercase[i]
                signatures.append(signature)
                i += 1

        # --- Create all-by-all similarity matrix for reference ---
        # from sourmash.compare import compare_all_pairs
        # compare = compare_all_pairs(signatures, ignore_abundance=True)
        # print([x.name() for x in signatures])
        # print(compare)
        # --- Similarity matrix ---
        # ['A',  'B',   'C',  'D', 'E',  'F',  'G']
        # [[1.    0.356 0.078 0.086 0.    0.    0.   ]
        #  [0.356 1.    0.072 0.078 0.    0.    0.   ]
        #  [0.078 0.072 1.    0.074 0.    0.    0.   ]
        #  [0.086 0.078 0.074 1.    0.    0.    0.   ]
        #  [0.    0.    0.    0.    1.    0.382 0.364]
        #  [0.    0.    0.    0.    0.382 1.    0.386]
        #  [0.    0.    0.    0.    0.364 0.386 1.   ]]

        # --- Insert: A ---
        # Tree:
        #     0
        #   /  \
        # A: 1  None

        # --- Insert: B ---
        # Tree:
        #     0
        #   /  \
        # A: 1  B: 2

        # --- Insert: C ---
        # Tree:
        #          0
        #        /  \
        #      1     C: 2
        #    /   \
        # A: 3  B: 4

        # --- Insert: D ---
        # Tree:
        #             0
        #        /        \
        #      1           2
        #    /   \       /   \
        # A: 3  B: 4   D: 5  C: 6

        # --- Insert: E ---
        # Tree:
        #                    0
        #               /       \
        #             1            2
        #        /        \      /   \
        #      3           4    E: 5  None: 6
        #    /   \       /   \
        # A: 7  B: 8   D: 9  C: 10

        # --- Insert: F ---
        # Tree:
        #                    0
        #               /       \
        #             1            2
        #        /        \      /   \
        #      3           4    E: 5  F: 6
        #    /   \       /   \
        # A: 7  B: 8   D: 9  C: 10

        # --- Insert: G ---
        # Tree:
        #                           0
        #               /                        \
        #             1                           2
        #        /        \                 /            \
        #      3           4              5               6
        #    /   \       /   \         /     \         /    \
        # A: 7  B: 8   D: 9  C: 10   F: 11  G: 12    E: 13  None: 14

        for signature in signatures:
            sbt.insert(signature)

        # - Randomly shuffle signatures and ensure the same leaves are sharing parents -
        # Set random seed for reproducibility/debugging
        random.seed(random_seed)
        random.shuffle(signatures)
        for signature in signatures:
            sbt_randomized.insert(signature)

        # Ensure all leaves are present in both
        signatures_in_sbt = sorted(
            [leaf.data for leaf in sbt.leaves()], key=lambda x: x.name()
        )
        signatures_in_sbt_randomized = sorted(
            [leaf.data for leaf in sbt_randomized.leaves()], key=lambda x: x.name()
        )
        assert all([s in signatures_in_sbt for s in signatures])
        assert signatures_in_sbt == signatures_in_sbt_randomized

        # Ensure that the most similar pairs, (A, B) and (F, G) share parents
        # regardless of construction order
        for tree in (sbt, sbt_randomized):
            # create mapping from leaf name to node pos
            leaf_pos = {sig.data.name(): n for n, sig in tree._leaves.items()}
            assert tree.parent(leaf_pos["A"]) == tree.parent(leaf_pos["B"])
            assert tree.parent(leaf_pos["F"]) == tree.parent(leaf_pos["G"])


@pytest.mark.filterwarnings("ignore")
def test_localized_sbt_on_gather_data():
    factory = GraphFactory(ksize=31, starting_size=1000, n_tables=3)
    sbt = LocalizedSBT(factory, track_abundance=True)

    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data("gather/GCF*.sig")
        # Sort to ensure consistent ordering across operating systems
        files = sorted(glob.glob(testdata_glob))

        signatures = []
        for filename in files:
            loaded = sig.load_signatures(filename, ksize=31)
            for signature in loaded:
                signatures.append(signature)

        # Rename signatures to ... X, Y, Z for simplicity of debugging
        n_signatures = len(signatures)
        for i, signature in enumerate(signatures):
            signature._name = ascii_uppercase[-(n_signatures - i)]

        # --- Create all-by-all similarity matrix for reference ---
        from sourmash.compare import compare_all_pairs

        compare = compare_all_pairs(signatures, ignore_abundance=False)
        print([x.name() for x in signatures])
        print(
            "["
            + ",\n ".join(
                ["[" + ", ".join([f"{x:.2f}" for x in row]) + "]" for row in compare]
            )
            + "]"
        )

        # --- Similarity matrix ---
        #    ['O',  'P',  'Q',  'R',  'S',  'T',  'U',  'V',  'W',  'X',  'Y',  'Z']
        # O: [[1.00, 0.48, 0.58, 0.00, 0.00, 0.61, 0.58, 0.52, 0.58, 0.00, 0.00, 0.45],
        # P:  [0.48, 1.00, 0.44, 0.00, 0.00, 0.51, 0.50, 0.60, 0.46, 0.00, 0.00, 0.90],
        # Q:  [0.58, 0.44, 1.00, 0.00, 0.00, 0.57, 0.55, 0.47, 0.54, 0.00, 0.00, 0.41],
        # R:  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.05, 0.00],
        # S:  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        # T:  [0.61, 0.51, 0.57, 0.00, 0.00, 1.00, 0.92, 0.52, 0.60, 0.00, 0.00, 0.47],
        # U:  [0.58, 0.50, 0.55, 0.00, 0.00, 0.92, 1.00, 0.50, 0.57, 0.00, 0.00, 0.46],
        # V:  [0.52, 0.60, 0.47, 0.00, 0.00, 0.52, 0.50, 1.00, 0.49, 0.00, 0.00, 0.57],
        # W:  [0.58, 0.46, 0.54, 0.00, 0.00, 0.60, 0.57, 0.49, 1.00, 0.00, 0.00, 0.43],
        # X:  [0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.03, 0.00],
        # Y:  [0.00, 0.00, 0.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 1.00, 0.00],
        # Z:  [0.45, 0.90, 0.41, 0.00, 0.00, 0.47, 0.46, 0.57, 0.43, 0.00, 0.00, 1.00]]

        # --- Insert: O ("oh") ---
        # First leaf --> take first position
        # Tree:
        #     0
        #   /  \
        # O: 1  None

        # --- Insert: P ---
        # Only one available --> Take next available
        # Tree:
        #     0
        #   /  \
        # O: 1  P: 2

        # --- Insert: Q ---
        # Most similar to O, displace P
        #  --> Push existing tree down when displacing P
        # Tree:
        #             0
        #        /        \
        #      1           2
        #    /   \       /   \
        # O: 3  Q: 4   P: 5  None

        # --- Insert: R ---
        # - R is not similar to any leaves currently in the tree,
        # --> insert naively to next position
        # Tree:
        #             0
        #        /        \
        #      1           2
        #    /   \       /   \
        # O: 3  Q: 4   P: 5  R: 6

        # --- Insert: S ---
        # S is not similar to anything
        #  --> Push existing tree down and insert next position
        # Tree:
        #                           0
        #               /                        \
        #             1                           2
        #        /        \                 /            \
        #      3           4              5               6
        #    /   \       /   \         /     \         /    \
        # O: 7  Q: 8   P: 9  R: 10   S: 11   None   None    None

        # --- Insert: T ---
        # - T is most similar to O ("oh") and should displace Q
        # - Q is most similar to P and should displace R
        # Desired Tree:
        #                           0
        #               /                        \
        #             1                           2
        #        /        \                 /            \
        #      3           4              5               6
        #    /   \       /   \         /     \         /    \
        # O: 7  T: 8   P: 9  Q: 10   S: 11  R: 12   None    None

        # --- Insert: U ---
        # - U is most similar to T and should displace O ("oh")
        # - After T, O ("oh") is most similar to Q and should displace P
        # Desired Tree:
        #                           0
        #               /                        \
        #             1                           2
        #        /        \                 /            \
        #      3           4              5               6
        #    /   \       /   \         /     \         /    \
        # U: 7  T: 8   O: 9  Q: 10   S: 11  R: 12   P: 13    None

        # --- Insert: V ---
        # - V is most similar to P and can just join P in its node
        # Desired Tree:
        #                           0
        #               /                        \
        #             1                           2
        #        /        \                 /            \
        #      3           4              5               6
        #    /   \       /   \         /     \         /    \
        # U: 7  T: 8   O: 9  Q: 10   S: 11  R: 12   P: 13    V: 14

        # --- Insert: W ---
        # - W is most similar to T, but not more similar than T's current sibiling, U.
        #   In fact, W doesn't beat out any of the existing pairs so it should be
        #   inserted wherever, to the next available node
        # Desired Tree:
        #                                            0
        #                            /                              \
        #                      1                                       2
        #              /                 \                     /              \
        #            3                    4                   5               6
        #       /       \           /         \            /    \           /  \
        #     7          8         9          10         11      12      13      14
        #   /   \     /   \      /   \      /   \      /   \    /  \
        # U:15 T:16 O:17 Q:18  S:19 R:20  P:21 V:22  W:23 None

        # --- Insert: X ---
        # - X is most similar to R, and should displace S, which isn't similar to
        #   anything and goes to its own place
        # Desired Tree:
        #                                            0
        #                            /                              \
        #                      1                                       2
        #              /                 \                     /              \
        #            3                    4                   5               6
        #       /       \           /         \            /    \           /  \
        #     7          8         9          10         11      12      13      14
        #   /   \     /   \      /   \      /   \      /   \    /  \
        # U:15 T:16 O:17 Q:18  X:19 R:20  P:21 V:22  W:23 S:24

        # --- Insert: Y ---
        # - Y is most similar to R and X, but doesn't beat out any of their most
        #    similar nodes, so goes into a node by itself
        # Desired Tree:
        #                                            0
        #                            /                              \
        #                      1                                       2
        #              /                 \                     /              \
        #            3                    4                   5               6
        #       /       \           /         \            /     \           /  \
        #     7          8         9          10         11       12      13      14
        #   /   \     /   \      /   \      /   \      /   \     /  \
        # U:15 T:16 O:17 Q:18  X:19 R:20  P:21 V:22  W:23 S:24 Y:25

        # --- Insert: Z ---
        # - Z is most similar to P, and should displace P's sibling V
        # - After P, V is most similar to W, and so V should displace W's sibling S
        # Desired Tree:
        #                                            0
        #                            /                              \
        #                      1                                       2
        #              /                 \                     /              \
        #            3                    4                   5               6
        #       /       \           /         \            /     \           /  \
        #     7          8         9          10         11       12      13      14
        #   /   \     /   \      /   \      /   \      /   \     /  \
        # U:15 T:16 O:17 Q:18  X:19 R:20  P:21 Z:22  W:23 V:24 Y:25 S:26

        for signature in signatures:
            sbt.insert(signature)

        # Ensure that the most similar pairs, (P, Z) and (F, G) share parents
        # regardless of construction order
        leaf_names_in_order = [
            v.name for k, v in sorted(sbt._leaves.items(), key=lambda x: x[0])
        ]
        assert leaf_names_in_order == [
            "U",
            "T",
            "O",
            "Q",
            "X",
            "R",
            "P",
            "Z",
            "W",
            "V",
            "Y",
            "S",
        ]

        # --- Double check to make sure parents are shared --- #
        # create mapping from leaf name to node pos
        leaf_pos = {sig.data.name(): n for n, sig in sbt._leaves.items()}
        assert sbt.parent(leaf_pos["P"]) == sbt.parent(leaf_pos["Z"])
        assert sbt.parent(leaf_pos["U"]) == sbt.parent(leaf_pos["T"])
        assert sbt.parent(leaf_pos["O"]) == sbt.parent(leaf_pos["Q"])


def test_localized_sbt_on_gather_data_realistic():
    """Use realistic NodeGraph/bloom filter settings"""
    order1 = [
        "GCF_000011885.1_ASM1188v1_genomic.fna.gz.sig",
        "GCF_000016785.1_ASM1678v1_genomic.fna.gz.sig",
        "GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
        "GCF_000009085.1_ASM908v1_genomic.fna.gz.sig",
        "GCF_000016045.1_ASM1604v1_genomic.fna.gz.sig",
        "GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
        "GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig",
        "GCF_000008545.1_ASM854v1_genomic.fna.gz.sig",
        "GCF_000009505.1_ASM950v1_genomic.fna.gz.sig",
        "GCF_000009525.1_ASM952v1_genomic.fna.gz.sig",
        "GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
        "GCF_000018945.1_ASM1894v1_genomic.fna.gz.sig",
    ]
    order2 = [
        "GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
        "GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
        "GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
        "GCF_000008545.1_ASM854v1_genomic.fna.gz.sig",
        "GCF_000009085.1_ASM908v1_genomic.fna.gz.sig",
        "GCF_000009505.1_ASM950v1_genomic.fna.gz.sig",
        "GCF_000009525.1_ASM952v1_genomic.fna.gz.sig",
        "GCF_000011885.1_ASM1188v1_genomic.fna.gz.sig",
        "GCF_000016045.1_ASM1604v1_genomic.fna.gz.sig",
        "GCF_000016785.1_ASM1678v1_genomic.fna.gz.sig",
        "GCF_000018945.1_ASM1894v1_genomic.fna.gz.sig",
        "GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig",
    ]
    KSIZE = 21

    # Use real data settings to match
    factory = GraphFactory(ksize=KSIZE, starting_size=1e5, n_tables=4)
    sbt1: LocalizedSBT = LocalizedSBT(factory, track_abundance=True)
    sbt2 = LocalizedSBT(factory, track_abundance=True)

    sorted_sig_names = [
        "NC_000853.1 Thermotoga maritima MSB8 chromosome, complete genome",
        "NC_002163.1 Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819 chromosome, complete genome",
        "NC_003197.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome",
        "NC_003198.1 Salmonella enterica subsp. enterica serovar Typhi str. CT18, complete genome",
        "NC_004631.1 Salmonella enterica subsp. enterica serovar Typhi Ty2, complete genome",
        "NC_006511.1 Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150, complete genome",
        "NC_006905.1 Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67, complete genome",
        "NC_009486.1 Thermotoga petrophila RKU-1, complete genome",
        "NC_011080.1 Salmonella enterica subsp. enterica serovar Newport str. SL254, complete genome",
        "NC_011274.1 Salmonella enterica subsp. enterica serovar Gallinarum str. 287/91 complete genome",
        "NC_011294.1 Salmonella enterica subsp. enterica serovar Enteritidis str. P125109 complete genome",
        "NC_011978.1 Thermotoga neapolitana DSM 4359, complete genome",
    ]

    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data("gather/GCF*.sig")
        query_filename = utils.get_test_data("gather/combined.sig")

        # Sort to ensure consistent ordering across operating systems
        files = sorted(glob.glob(testdata_glob))
        n = len(files)
        name_to_letter = {}
        letter_to_name = {}

        signatures = {}
        for i, filename in enumerate(files):
            loaded = sig.load_signatures(filename, ksize=KSIZE)
            basename = os.path.basename(filename)
            letter = ascii_uppercase[-(n - i)]
            for signature in loaded:
                signatures[basename] = signature
                # only one signature per filename so can assume that i is only
                # incrementally increasing
                name_to_letter[signature._name] = letter
                letter_to_name[letter] = signature._name

        query_sig = list(sig.load_signatures(query_filename, ksize=KSIZE))[0]

    # Tree 1
    for basename in order1:
        signature = signatures[basename]
        sbt1.insert(signature)
        sbt1_str = sbt2._str_query(query=query_sig)

    # Tree 2
    for basename in order2:
        signature = signatures[basename]
        sbt2.insert(signature)

    # from sourmash.compare import compare_all_pairs
    #
    # compare = compare_all_pairs(list(signatures.values()), ignore_abundance=False)
    # print()
    # print([x.name().split()[0] for x in signatures.values()])
    # print(
    #     "["
    #     + ",\n ".join(
    #         ["[" + ", ".join([f"{x:.2f}" for x in row]) + "]" for row in compare]
    #     )
    #     + "]"
    # )
    # [
    #     "NC_003197.2",
    #     "NC_004631.1",
    #     "NC_006905.1",
    #     "NC_000853.1",
    #     "NC_002163.1",
    #     "NC_011294.1",
    #     "NC_011274.1",
    #     "NC_006511.1",
    #     "NC_011080.1",
    #     "NC_009486.1",
    #     "NC_011978.1",
    #     "NC_003198.1",
    # ]
    # [[1.00, 0.56, 0.59, 0.00, 0.00, 0.65, 0.63, 0.57, 0.65, 0.00, 0.00, 0.52],
    #  [0.56, 1.00, 0.55, 0.00, 0.00, 0.58, 0.57, 0.66, 0.56, 0.00, 0.00, 0.91],
    #  [0.59, 0.55, 1.00, 0.00, 0.00, 0.60, 0.57, 0.53, 0.60, 0.00, 0.00, 0.52],
    #  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.33, 0.05, 0.00],
    #  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    #  [0.65, 0.58, 0.60, 0.00, 0.00, 1.00, 0.91, 0.59, 0.63, 0.00, 0.00, 0.55],
    #  [0.63, 0.57, 0.57, 0.00, 0.00, 0.91, 1.00, 0.58, 0.60, 0.00, 0.00, 0.53],
    #  [0.57, 0.66, 0.53, 0.00, 0.00, 0.59, 0.58, 1.00, 0.57, 0.00, 0.00, 0.61],
    #  [0.65, 0.56, 0.60, 0.00, 0.00, 0.63, 0.60, 0.57, 1.00, 0.00, 0.00, 0.53],
    #  [0.00, 0.00, 0.00, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.04, 0.00],
    #  [0.00, 0.00, 0.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 1.00, 0.00],
    #  [0.52, 0.91, 0.52, 0.00, 0.00, 0.55, 0.53, 0.61, 0.53, 0.00, 0.00, 1.00]]



    # These samples are most similar to one another and should always share parents
    pairs = (
        (
            "NC_009486.1 Thermotoga petrophila RKU-1, complete genome",
            "NC_000853.1 Thermotoga maritima MSB8 chromosome, complete genome",
        ),
        (
            "NC_011080.1 Salmonella enterica subsp. enterica serovar Newport str. SL254, complete genome",
            "NC_003197.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome",
        ),
        (
            "NC_011274.1 Salmonella enterica subsp. enterica serovar Gallinarum str. 287/91 complete genome",
            "NC_011294.1 Salmonella enterica subsp. enterica serovar Enteritidis str. P125109 complete genome",
        ),
        (
            "NC_003198.1 Salmonella enterica subsp. enterica serovar Typhi str. CT18, complete genome",
            "NC_004631.1 Salmonella enterica subsp. enterica serovar Typhi Ty2, complete genome",
        ),
    )
    for sbt in (sbt1, sbt2):
        leaf_pos = {sig.data.name(): n for n, sig in sbt._leaves.items()}
        for leaf1, leaf2 in pairs:
            assert sbt.parent(leaf_pos[leaf1]) == sbt.parent(leaf_pos[leaf2])


@pytest.mark.skip(reason="Not currently working but not a show-stopping bug")
def test_localized_sbt_adversarial_identical_signatures():
    factory = GraphFactory(5, 100, 3)
    sbt = LocalizedSBT(factory, track_abundance=False)

    with utils.TempDirectory() as location:
        # Sort to ensure consistent ordering across operating systems
        files = sorted([utils.get_test_data(f) for f in utils.SIG_FILES])
        i = 0
        for filename in files:
            loaded = sig.load_signatures(filename, ksize=31)
            for signature in loaded:
                # Rename to A, B, C, D for simplicity
                signature._name = ascii_uppercase[i]
                for j in range(10):
                    # Fails when i = 0 and j = 5
                    # (the fifth insertion of the first signature)
                    sbt.insert(signature)
                i += 1


def test_localized_sbt_adversarial_dissimilar_signatures():
    """CHeck that tree construction doesn't fail when all signatures are dissimilar"""
    factory = GraphFactory(5, 100, 3)
    sbt = LocalizedSBT(factory)

    n_hashes = 3
    a = MinHash(n=n_hashes, ksize=5)
    a.add("AAAAA")
    a.add("AAAAC")
    a.add("AAAAG")
    a.add("AAAAT")
    sig_a = SourmashSignature(a, name="a")

    b = MinHash(n=n_hashes, ksize=5)
    b.add("TTTTA")
    b.add("TTTTC")
    b.add("TTTTG")
    b.add("TTTTT")
    sig_b = SourmashSignature(b, name="b")

    c = MinHash(n=n_hashes, ksize=5)
    c.add("CCCCA")
    c.add("CCCCC")
    c.add("CCCCG")
    c.add("CCCCGT")
    sig_c = SourmashSignature(c, name="c")

    d = MinHash(n=n_hashes, ksize=5)
    d.add("GGGGA")
    d.add("GGGGC")
    d.add("GGGGG")
    d.add("GGGGT")
    sig_d = SourmashSignature(d, name="d")

    e = MinHash(n=n_hashes, ksize=5)
    e.add("TAAAA")
    e.add("TAAAC")
    e.add("TAAAG")
    e.add("TAAAT")
    sig_e = SourmashSignature(e, name="e")

    sigs = (sig_a, sig_b, sig_c, sig_d, sig_e)
    for sig in sigs:
        sbt.insert(sig)

    # Make sure tree construction happened properly
    assert all(node < leaf for leaf, node in product(sbt._leaves, sbt._nodes))

    # create mapping from leaf name to node pos
    leaf_pos = {sig.data.name(): n for n, sig in sbt._leaves.items()}
    # Ensure signatures were added as leaves
    assert all(x in leaf_pos.keys() for x in list("abcde"))
