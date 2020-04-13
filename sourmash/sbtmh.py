from __future__ import print_function
from __future__ import division

import math
from io import BytesIO, TextIOWrapper
import sys

from .sbt import Leaf, SBT, GraphFactory
from . import signature


def load_sbt_index(filename, print_version_warning=True):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning)


def create_sbt_index(bloom_filter_size=1e5, n_children=2):
    "Create an empty SBT index."
    factory = GraphFactory(1, bloom_filter_size, 4)
    tree = SBT(factory, d=n_children)
    return tree


def search_sbt_index(tree, query, threshold):
    """\
    Search an SBT index `tree` with signature `query` for matches above
    `threshold`.

    Usage:

        for match_sig, similarity in search_sbt_index(tree, query, threshold):
           ...
    """
    for leaf in tree.find(search_minhashes, query, threshold, unload_data=True):
        similarity = query.similarity(leaf.data)
        yield leaf.data, similarity


class SigLeaf(Leaf):
    def __str__(self):
        return '**Leaf:{name} -> {metadata}'.format(
                name=self.name, metadata=self.metadata)

    def save(self, path):
        # this is here only for triggering the property load
        # before we reopen the file (and overwrite the previous
        # content...)
        self.data

        buf = BytesIO()
        with TextIOWrapper(buf) as out:
            signature.save_signatures([self.data], out)
            out.flush()
            return self.storage.save(path, buf.getvalue())

    def update(self, parent):
        mh = self.data.minhash
        for v in mh.get_mins():
            parent.data.count(v)
        min_n_below = parent.metadata.get('min_n_below', sys.maxsize)
        min_n_below = min(len(mh), min_n_below)

        if min_n_below == 0:
            min_n_below = 1

        parent.metadata['min_n_below'] = min_n_below

    @property
    def data(self):
        if self._data is None:
            buf = BytesIO(self.storage.load(self._path))
            self._data = signature.load_one_signature(buf)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data


### Search functionality.

def _max_jaccard_underneath_internal_node(node, hashes):
    """\
    calculate the maximum possibility similarity score below
    this node, based on the number of matches in 'hashes' at this node,
    divided by the smallest minhash size below this node.

    This should yield be an upper bound on the Jaccard similarity
    for any signature below this point.
    """
    if len(hashes) == 0:
        return 0.0

    # count the maximum number of hash matches beneath this node
    get = node.data.get
    matches = sum(1 for value in hashes if get(value))

    # get the size of the smallest collection of hashes below this point
    min_n_below = node.metadata.get('min_n_below', -1)

    if min_n_below == -1:
        raise Exception('cannot do similarity search on this SBT; need to rebuild.')

    # max of numerator divided by min of denominator => max Jaccard
    max_score = float(matches) / min_n_below

    return max_score


def search_minhashes(node, sig, threshold, results=None):
    """\
    Default tree search function, searching for best Jaccard similarity.
    """
    mins = sig.minhash.get_mins()
    score = 0

    if isinstance(node, SigLeaf):
        score = node.data.minhash.similarity(sig.minhash)
    else:  # Node minhash comparison
        score = _max_jaccard_underneath_internal_node(node, mins)

    if results is not None:
        results[node.name] = score

    if score >= threshold:
        return 1

    return 0


class SearchMinHashesFindBest(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        mins = sig.minhash.get_mins()
        score = 0

        if isinstance(node, SigLeaf):
            score = node.data.minhash.similarity(sig.minhash)
        else:  # internal object, not leaf.
            score = _max_jaccard_underneath_internal_node(node, mins)

        if results is not None:
            results[node.name] = score

        if score >= threshold:
            # have we done better than this elsewhere? if yes, truncate.
            if score > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = score
                return 1

        return 0


def search_minhashes_containment(node, sig, threshold, results=None, downsample=True):
    mins = sig.minhash.get_mins()

    if isinstance(node, SigLeaf):
        matches = node.data.minhash.count_common(sig.minhash, downsample)
    else:  # Node or Leaf, Nodegraph by minhash comparison
        get = node.data.get
        matches = sum(1 for value in mins if get(value))

    if results is not None:
        results[node.name] = float(matches) / len(mins)

    if len(mins) and float(matches) / len(mins) >= threshold:
        return 1
    return 0


class GatherMinHashes(object):
    def __init__(self):
        self.best_match = 0

    def search(self, node, query, threshold, results=None):
        score = 0
        if not len(query.minhash):
            return 0

        if isinstance(node, SigLeaf):
            matches = query.minhash.count_common(node.data.minhash, True)
        else:  # Nodegraph by minhash comparison
            mins = query.minhash.get_mins()
            get = node.data.get
            matches = sum(1 for value in mins if get(value))

        score = float(matches) / len(query.minhash)

        # store results if we have passed in an appropriate dictionary
        if results is not None:
            results[node.name] = score

        # have we done better than this? if no, truncate searches below.
        if score >= self.best_match:
            # update best if it's a leaf node...
            if isinstance(node, SigLeaf):
                self.best_match = score
            return 1

        return 0


### --- Nearest Neighbor Functionality ---

def get_leaves(tree):
    for i, node in tree.nodes.items():
        if isinstance(node, SigLeaf) or isinstance(node, Leaf):
            yield i, node


def nearest_neighbors(tree, n_neighbors, ignore_abundance, downsample,
                      min_similarity=0.0):
    n_parent_levels = math.log2(n_neighbors) + 1

    for position1, leaf1 in get_leaves(tree):
        n = 1
        upper_internal_node = tree.parent(position1)
        while n < n_parent_levels:
            upper_internal_node = tree.parent(upper_internal_node.pos)
            n += 1
        #         print("upper_internal_node:", upper_internal_node)
        leaves = tree.get_leaves_under(upper_internal_node.pos)

        similarities = []
        for leaf2 in leaves:
            if leaf2 == leaf1:
                continue
            similarity = leaf1.data.similarity(leaf2.data,
                                              ignore_abundance=ignore_abundance,
                                              downsample=downsample)
            if similarity > min_similarity:
                similarities.append(
                    [leaf1.data.name(), leaf2.data.name(), similarity])
        adjacent = sorted(similarities, key=lambda x: x[1])[-n_neighbors:]
        for adjacency in adjacent:
            yield adjacency


def adjacency_to_knn(adjacencies, tree):
    leaf_to_index = dict(
        (node.data.name(), i) for i, node in enumerate(get_leaves(tree)))
    index_to_leaf = dict(zip(leaf_to_index.values(), leaf_to_index.keys()))
    knn_indices = []
    knn_dists = []

    for u, items in itertools.groupby(adjacencies, key=lambda x: x[0]):
        knn_indices_line = []
        knn_dists_line = []
        for u, v, similarity in items:
            knn_indices_line.append(leaf_to_index[v])
            knn_dists_line.append(1 - similarity)
        knn_indices.append(knn_indices_line)
        knn_dists.append(knn_dists_line)
