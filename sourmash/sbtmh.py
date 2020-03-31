from __future__ import division
from __future__ import print_function

import sys
from io import BytesIO, TextIOWrapper

from . import signature
from .sbt import Leaf, SBT, GraphFactory


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


class LocalizedSBT(SBT):
    """A Sequence Bloom Tree implementation which guarantees new leaves are plaaced in

    The default node is a Bloom Filter (like the original implementation),
    and the leaves are MinHash leaf class (in the sourmash.sbtmh.SigLeaf class)

    Parameters
    ----------
    factory: Factory
        Callable for generating new datastores for internal nodes.
    d: int
        Number of children for each internal node. Defaults to 2 (a binary tree)
    storage: Storage, default: None
        A Storage is any place where we can save and load data for the nodes.
        If set to None, will use a FSStorage.

    Notes
    -----
    We use two dicts to store the tree structure: One for the internal nodes,
    and another for the leaves (datasets).
    """

    def __init__(self, factory, d=2, storage=None, track_abundance=False,
                 do_containment=False):
        SBT.__init__(self, factory=factory, d=d, storage=storage)
        self.track_abundance = track_abundance
        self.ignore_abundance = not self.track_abundance
        self.do_containment = do_containment

    def new_node_pos(self, node):
        if not self._nodes:
            self.next_node = 1
            return 0

        if not self._leaves:
            self.next_node = 2
            return 1

        # Not an empty tree, can search
        most_similar_node = self.search(node.data, threshold=0, best_only=1,
                                        ignore_abundance=self.ignore_abundance,
                                        do_containment=self.do_containment)

        # Get parent of the most similar node
        localized_parent = self.parent(most_similar_node[0])

        # If the parent has one child: easy, insert the new child here
        children = self.children(localized_parent.pos)
        if children[1].node is None:
            return children[1].pos
        else:
            # If parent has two children, check if the other child is more similar to
            # the most_similar_node --> then no displacement is necessary
            children[1]

        # If the parent has two children .. then need to find a new place for the
        # displaced child

        # Insert node under this parent
        import pdb;
        pdb.set_trace()
        min_leaf = min(self._leaves.keys())

        next_internal_node = None
        if self.next_node <= min_leaf:
            for i in range(min_leaf):
                if all((i not in self._nodes,
                        i not in self._leaves,
                        i not in self._missing_nodes)):
                    next_internal_node = i
                    break

        if next_internal_node is None:
            self.next_node = max(self._leaves.keys()) + 1
        else:
            self.next_node = next_internal_node

        return self.next_node
