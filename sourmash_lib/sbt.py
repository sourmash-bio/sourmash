#!/usr/bin/env python
"""
A trial implementation of sequence bloom trees, Solomon & Kingsford, 2015.

This is a simple in-memory version where all of the graphs are in
memory at once; to move it onto disk, the graphs would need to be
dynamically loaded for each query.

To try it out, do::

    factory = GraphFactory(ksize, tablesizes, n_tables)
    root = Node(factory)

    graph1 = factory()
    # ... add stuff to graph1 ...
    leaf1 = Leaf("a", graph1)
    root.add_node(leaf1)

For example, ::

    # filenames: list of fa/fq files
    # ksize: k-mer size
    # tablesizes: Bloom filter table sizes
    # n_tables: Number of tables

    factory = GraphFactory(ksize, tablesizes, n_tables)
    root = Node(factory)

    for filename in filenames:
        graph = factory()
        graph.consume_fasta(filename)
        leaf = Leaf(filename, graph)
        root.add_node(leaf)

then define a search function, ::

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.data.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * len(seq)):
            return 1
        return 0
"""

from __future__ import print_function, unicode_literals, division

from collections import namedtuple, Mapping
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

    def new_node_pos(self, node):
        try:
            pos = self.nodes.index(None)
        except ValueError:
            # There aren't any empty positions left.
            # Extend array
            height = math.floor(math.log(len(self.nodes), self.d)) + 1
            self.nodes += [None] * int(self.d ** height)
            pos = self.nodes.index(None)
        return pos

    def add_node(self, node):
        pos = self.new_node_pos(node)

        if pos == 0:  # empty tree; initialize w/node.
            n = Node(self.factory, name="internal." + str(pos))
            self.nodes[0] = n
            pos = self.new_node_pos(node)

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

            c1, c2 = self.children(p.pos)[:2]

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

    def find(self, search_fn, *args, **kwargs):
        matches = []
        visited, queue = set(), [0]
        while queue:
            node_p = queue.pop(0)
            node_g = self.nodes[node_p]
            if node_g is None:
                continue

            if node_p not in visited:
                visited.add(node_p)
                if search_fn(node_g, *args):
                    if isinstance(node_g, Leaf):
                        matches.append(node_g)
                    elif isinstance(node_g, Node):
                        if kwargs.get('dfs', True):  # defaults search to dfs
                            for c in self.children(node_p):
                                queue.insert(0, c.pos)
                        else: # bfs
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
        version = 2
        basetag = os.path.basename(tag)
        dirprefix = os.path.dirname(tag)
        dirname = os.path.join(dirprefix, '.sbt.' + basetag)

        if not os.path.exists(dirname):
            os.makedirs(dirname)

        info = {}
        info['d'] = self.d
        info['version'] = version

        structure = {}
        for i, node in iter(self):
            if node is None:
                structure[i] = None
                continue

            basename = os.path.basename(node.name)
            data = {
                'filename': os.path.join('.sbt.' + basetag,
                                         '.'.join([basetag, basename, 'sbt'])),
                'name': node.name
            }
            if isinstance(node, Leaf):
                data['metadata'] = node.metadata

            node.save(os.path.join(dirprefix, data['filename']))
            structure[i] = data

        fn = tag + '.sbt.json'
        info['nodes'] = structure
        with open(fn, 'w') as fp:
            json.dump(info, fp)

        return fn

    @classmethod
    def load(cls, sbt_name, leaf_loader=None):
        dirname = os.path.dirname(sbt_name)
        sbt_name = os.path.basename(sbt_name)

        loaders = {
            1: cls._load_v1,
            2: cls._load_v2,
        }

        # @CTB hack: check to make sure khmer Nodegraph supports the
        # correct methods.
        x = khmer.Nodegraph(1, 1, 1)
        try:
            x.count(10)
        except TypeError:
            raise Exception("khmer version is too old; need >= 2.1.")

        if leaf_loader is None:
            leaf_loader = Leaf.load

        sbt_fn = sbt_name
        if not sbt_fn.endswith('.sbt.json'):
            sbt_fn = sbt_fn + '.sbt.json'
        with open(os.path.join(dirname, sbt_fn)) as fp:
            jnodes = json.load(fp)

        version = 1
        if isinstance(jnodes, Mapping):
            version = jnodes['version']

        return loaders[version](jnodes, leaf_loader, dirname)

    @staticmethod
    def _load_v1(jnodes, leaf_loader, dirname):

        if jnodes[0] is None:
            # TODO error!
            raise ValueError("Empty tree!")

        sbt_nodes = []

        sample_bf = os.path.join(dirname, jnodes[0]['filename'])
        ksize, tablesize, ntables = khmer.extract_nodegraph_info(sample_bf)[:3]
        factory = GraphFactory(ksize, tablesize, ntables)

        for jnode in jnodes:
            if jnode is None:
                sbt_nodes.append(None)
                continue

            if 'internal' in jnode['filename']:
                jnode['factory'] = factory
                sbt_node = Node.load(jnode, dirname)
            else:
                sbt_node = leaf_loader(jnode, dirname)

            sbt_nodes.append(sbt_node)

        tree = SBT(factory)
        tree.nodes = sbt_nodes

        return tree

    @classmethod
    def _load_v2(cls, info, leaf_loader, dirname):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if nodes[0] is None:
            raise ValueError("Empty tree!")

        sbt_nodes = []

        sample_bf = os.path.join(dirname, nodes[0]['filename'])
        k, size, ntables = khmer.extract_nodegraph_info(sample_bf)[:3]
        factory = GraphFactory(k, size, ntables)

        for i, node in sorted(nodes.items()):
            if node is None:
                sbt_nodes.append(None)
                continue

            if 'internal' in node['filename']:
                node['factory'] = factory
                sbt_node = Node.load(node, dirname)
            else:
                sbt_node = leaf_loader(node, dirname)

            sbt_nodes.append(sbt_node)

        tree = cls(factory, d=info['d'])
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

        for i, node in iter(self):
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

    def __iter__(self):
        for i, node in enumerate(self.nodes):
            yield (i, node)


class Node(object):
    "Internal node of SBT."

    def __init__(self, factory, name=None, fullpath=None):
        self.name = name
        self._factory = factory
        self._data = None
        self._filename = fullpath

    def __str__(self):
        return '*Node:{name} [occupied: {nb}, fpr: {fpr:.2}]'.format(
                name=self.name, nb=self.data.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    def save(self, filename):
        self.data.save(filename)

    @property
    def data(self):
        if self._data is None:
            if self._filename is None:
                self._data = self._factory()
            else:
                self._data = khmer.load_nodegraph(self._filename)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    @staticmethod
    def load(info, dirname):
        filename = os.path.join(dirname, info['filename'])
        new_node = Node(info['factory'], name=info['name'], fullpath=filename)
        return new_node


class Leaf(object):
    def __init__(self, metadata, data=None, name=None, fullpath=None):
        self.metadata = metadata
        if name is None:
            name = metadata
        self.name = name
        self._data = data
        self._filename = fullpath

    def __str__(self):
        return '**Leaf:{name} [occupied: {nb}, fpr: {fpr:.2}] -> {metadata}'.format(
                name=self.name, metadata=self.metadata,
                nb=self.data.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    @property
    def data(self):
        if self._data is None:
            # TODO: what if self._filename is None?
            self._data = khmer.load_nodegraph(self._filename)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    def save(self, filename):
        self.data.save(filename)

    def update(self, parent):
        parent.data.update(self.data)

    @classmethod
    def load(cls, info, dirname):
        filename = os.path.join(dirname, info['filename'])
        return cls(info['metadata'], name=info['name'], fullpath=filename)


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
