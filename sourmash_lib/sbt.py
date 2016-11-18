#!/usr/bin/env python

"""
A trial implementation of sequence bloom trees, Solomon & Kingsford, 2015.

This is a simple in-memory version where all of the graphs are in
memory at once; to move it onto disk, the graphs would need to be
dynamically loaded for each query.

To try it out, do::

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    graph1 = factory.create_nodegraph()
    # ... add stuff to graph1 ...
    leaf1 = Leaf("a", graph1)
    root.add_node(leaf1)

For example, ::

    # filenames: list of fa/fq files
    # ksize: k-mer size
    # tablesizes: Bloom filter table sizes

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    for filename in filenames:
        graph = factory.create_nodegraph()
        graph.consume_fasta(filename)
        leaf = Leaf(filename, graph)
        root.add_node(leaf)

then define a search function, ::

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.graph.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * len(seq)):
            return 1
        return 0
"""

from __future__ import print_function, unicode_literals

from collections import namedtuple
import hashlib
import json
import math
import os
import random
import shutil
from tempfile import NamedTemporaryFile

import khmer
from khmer import khmer_args
from random import randint
from numpy import array


NodePos = namedtuple("NodePos", ["pos", "node"])


def GraphFactory(ksize, starting_size, n_tables):
    "Build new nodegraphs (Bloom filters) of a specific (fixed) size."

    def create_nodegraph():
        return khmer.Nodegraph(ksize, starting_size, n_tables)

    return create_nodegraph


class SBT(object):

    def __init__(self, factory, d=2):
        self.factory = factory
        self.nodes = [None]
        self.d = d

    def add_node(self, node):
        try:
            pos = self.nodes.index(None)
        except ValueError:
            # There aren't any empty positions left.
            # Extend array
            current_size = len(self.nodes)
            # TODO: this is too much, figure out the lower bound
            self.nodes += [None] * (2 * (current_size + 1) * self.d)
            pos = self.nodes.index(None)

        if pos == 0:  # empty tree
            self.nodes[0] = node
            return

        # Cases:
        # 1) parent is a Leaf (already covered)
        # 2) parent is a Node (with empty position available)
        #    - add Leaf, update parent
        # 3) parent is a Node (no position available)
        #    - this is covered by case 1
        p = self.parent(pos)
        if isinstance(p.node, Leaf):
            # Create a new internal node
            # node and parent are children of new internal node
            n = Node(self.factory, name="internal." + str(p.pos))
            self.nodes[p.pos] = n

            c1, c2, *remainder = self.children(p.pos)

            self.nodes[c1.pos] = p.node
            self.nodes[c2.pos] = node

            for child in (p.node, node):
                child.update(n)
        elif isinstance(p.node, Node):
            self.nodes[pos] = node
            node.update(p.node)

        # update all parents!
        p = self.parent(p.pos)
        while p:
            node.update(p.node)
            p = self.parent(p.pos)


    def find(self, search_fn, *args):
        matches = []
        visited, queue = set(), [0]
        while queue:
            node_p = queue.pop(0)
            node_g = self.nodes[node_p]
            if node_p not in visited and node_g is not None:
                visited.add(node_p)
                if search_fn(node_g, *args):
                    if isinstance(node_g, Leaf):
                        matches.append(node_g)
                    elif isinstance(node_g, Node):
                        queue.extend(c.pos for c in self.children(node_p))
        return matches

    def parent(self, pos):
        if pos == 0:
            return None
        p = int(math.floor((pos - 1) / self.d))
        return NodePos(p, self.nodes[p])

    def children(self, pos):
        return [self.child(pos, c) for c in range(self.d)]

    def child(self, parent, pos):
        cd = self.d * parent + pos + 1
        return NodePos(cd, self.nodes[cd])

    def save(self, tag):
        dirname = '.sbt.' + tag

        if not os.path.exists(dirname):
            os.makedirs(dirname)

        structure = []
        for node in self.nodes:
            if node is None:
                structure.append(None)
                continue

            data = {
                'filename': os.path.join('.sbt.' + tag,
                                         '.'.join([tag, node.name, 'sbt'])),
                'name': node.name
            }
            if isinstance(node, Leaf):
                data['metadata'] = node.metadata

            node.save(data['filename'])
            structure.append(data)

        fn = tag + '.sbt.json'
        with open(fn, 'w') as fp:
            json.dump(structure, fp)

        return fn

    @staticmethod
    def load(sbt_fn, leaf_loader=None):
        if leaf_loader is None:
            leaf_loader = Leaf.load

        with open(sbt_fn) as fp:
            nodes = json.load(fp)

        if nodes[0] is None:
            # TODO error!
            raise ValueError("Empty tree!")

        sbt_nodes = []

        ksize, tablesize, ntables, _, _, _ = khmer.extract_nodegraph_info(nodes[0]['filename'])
        factory = GraphFactory(ksize, tablesize, ntables)

        for node in nodes:
            if node is None:
                sbt_nodes.append(None)
                continue

            if 'internal' in node['filename']:
                node['factory'] = factory
                new_node = Node.load(node)
            else:
                new_node = leaf_loader(node)

            sbt_nodes.append(new_node)

        tree = SBT(factory)
        tree.nodes = sbt_nodes

        return tree

    def print_dot(self):
        print("""
        digraph G {
        nodesep=0.3;
        ranksep=0.2;
        margin=0.1;
        node [shape=ellipse];
        edge [arrowsize=0.8];
        """)

        for i, node in enumerate(self.nodes):
            if node is None:
                continue

            p = self.parent(i)
            if p is not None:
                if isinstance(node, Leaf):
                    print('"', p.pos, '"', '->', '"', node.name, '";')
                else:
                    print('"', p.pos, '"', '->', '"', i, '";')
        print("}")

    def print(self):
        visited, stack = set(), [0]
        while stack:
            node_p = stack.pop()
            node_g = self.nodes[node_p]
            if node_p not in visited and node_g is not None:
                visited.add(node_p)
                depth = int(math.floor(math.log(node_p + 1, self.d)))
                print(" " * 4 * depth, node_g)
                if isinstance(node_g, Node):
                    stack.extend(c.pos for c in self.children(node_p)
                                       if c.pos not in visited)

class Node(object):
    "Internal node of SBT; has 0, 1, or 2 children."

    def __init__(self, factory, name=None):
        self.data = factory()
        self.name = name

    def __str__(self):
        return '*Node:{name} [occupied: {nb}, fpr: {fpr:.2}]'.format(
                name=self.name, nb=self.data.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    def save(self, filename):
        self.data.save(filename)

    @staticmethod
    def load(info):
        new_node = Node(info['factory'], name=info['name'])
        new_node.data = khmer.load_nodegraph(info['filename'])
        return new_node


class Leaf(object):
    def __init__(self, metadata, data, name=None):
        self.metadata = metadata
        if name is None:
            name = metadata
        self.name = name
        self.data = data

    def __str__(self):
        return '**Leaf:{name} [occupied: {nb}, fpr: {fpr:.2}] -> {metadata}'.format(
                name=self.name, metadata=self.metadata,
                nb=self.data.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    def save(self, filename):
        self.data.save(filename)

    def update(self, parent):
        parent.data.update(self.data)

    @staticmethod
    def load(info):
        data = khmer.load_nodegraph(info['filename'])
        return Leaf(info['metadata'], data, name=info['name'])


def filter_distance( filter_a, filter_b, n=1000 ) :
    """
    Compute a heuristic distance per bit between two Bloom
    filters.

    filter_a : First filter
    filter_b : Second filter
    n        : Number of positions to compare (in groups of 8)
    """
    A = filter_a.graph.get_raw_tables()
    B = filter_b.graph.get_raw_tables()
    distance = 0
    for q,p in zip( A, B ) :
        a = array( q, copy=False )
        b = array( p, copy=False )
        for i in map( lambda x : randint( 0, len(a) ), range(n) ) :
            distance += sum( map( int, [ not bool((a[i]>>j)&1)
                                           ^ bool((b[i]>>j)&1)
                                         for j in range(8) ] ) )
    return distance / ( 8.0 * len(A) * n )


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
