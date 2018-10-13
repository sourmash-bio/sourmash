from __future__ import print_function, unicode_literals

from sourmash.index import LinearIndex
from sourmash_lib.sbt import SBT, GraphFactory, Leaf


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
