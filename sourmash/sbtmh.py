from __future__ import print_function
from __future__ import division

from io import BytesIO, TextIOWrapper
import sys

from .sbt import Leaf, Node, SBT, GraphFactory
from . import signature
from .logging import notify



def load_sbt_index(filename, print_version_warning=True):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning)


def create_sbt_index(bloom_filter_size=1e5, n_children=2, not_localized=False):
    "Create an empty SBT index."
    factory = GraphFactory(1, bloom_filter_size, 4)
    if not_localized or n_children != 2:
        if n_children != 2:
            notify("n_children != 2 in tree, defaulting to non-localized sequence "
                   "bloom tree (SBT) index")
        tree = SBT(factory, d=n_children)
    else:
        tree = LocalizedSBT(factory, d=n_children)
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


class LocalizedSBT(SBT):
    """An SBT implementation which inserts new leaves next to most similar existing leaf

    In this Sequence Bloom Tree (SBT) implementation, the default node is a Bloom Filter
    (like the original implementation), and the leaves are MinHash leaf class (in the
    sourmash.sbtmh.SigLeaf class). This is called "localized" because leaves that are
    similar to one another are constrained to sharing parents.

    Parameters
    ----------
    factory: Factory
        Callable for generating new datastores for internal nodes.
    d: int
        Number of children for each internal node. Defaults to 2 (a binary tree).
        Currently only implemented for d=2, binary trees
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
        if d != 2:
            raise NotImplementedError("LocalizedSBT is only implemented for when the "
                                      "number of children is 2, d=2")
        super().__init__(factory=factory, d=d, storage=storage)
        self.track_abundance = track_abundance
        self.ignore_abundance = not self.track_abundance
        self.do_containment = do_containment

    def find_most_similar_leaf(self, node):
            search_results = self.search(
                node.data,
                threshold=sys.float_info.epsilon,
                best_only=True,
                ignore_abundance=self.ignore_abundance,
                do_containment=self.do_containment,
                return_leaf=True
            )
            if len(search_results) == 1:
                best_result = search_results.pop()
            elif search_results:
                # Use the computed similarity to pick the best result
                # Note: if there are ties, this takes the first one (I think)
                best_result = max(search_results, key=lambda x: x[0])
            else:
                # no similarity overlap found; search_results empty
                best_result = None

            return best_result

    def get_sibling_of_similar_leaf(self, children, most_similar_leaf):
        # if most similar node has two children already, return node
        # of least similar child (displaced)
        # Get the leaf information of the other child
        if most_similar_leaf == children[0].node:
            other_child = children[1]
        elif most_similar_leaf == children[1].node:
            other_child = children[0]
        else:
            raise ValueError(
                "Neither children in node show up as most similar"
                " leaf. Something weird happened in search."
            )

        return other_child

    def get_child_nodes(self, children):
        return [c.node for c in children]

    def check_if_all_sigleafs(self, child_nodes):
        return all(
            isinstance(x, SigLeaf)
            for x in child_nodes
        )

    def get_siblings(self, grandparent, most_similar_parent):
        return [
            x for x in self.children(grandparent.pos)
            if x != most_similar_parent
        ]

    def new_node_pos(self, node):
        if not self._nodes:
            self.next_node = 1
            return 0

        if not self._leaves:
            self.next_node = 2
            return 1
        # Not an empty tree, can search

        # TODO: There is probably a way better way to write this logic - @olgabot
        if isinstance(node, SigLeaf):
            search_results = self.search(
                node.data, threshold=sys.float_info.epsilon, best_only=True,
                ignore_abundance=self.ignore_abundance,
                do_containment=self.do_containment, return_leaf=True)
            if len(search_results) == 1:
                best_result = search_results.pop()
            elif search_results:
                # Use the computed similarity to pick the best result
                # Note: if there are ties, this takes the first one (I think)
                best_result = max(search_results, key=lambda x: x[0])
            else:
                self.next_node = self._insert_next_position(self.next_node)
                return self.next_node

            new_leaf_similarity, most_similar_leaf, most_similar_pos = best_result

            # Get parent of the most similar node
            most_similar_parent = self.parent(most_similar_pos)

            # If the parent has one child: easy, insert the new child here
            children = self.children(most_similar_parent.pos)
            if children[1].node is None:
                # Use the default next node position
                self.next_node = self._insert_next_position(self.next_node)
            else:
                # If parent has two children, check if the other child is more similar
                # to the most_similar_leaf --> then no displacement is necessary
                other_child = self.get_sibling_of_similar_leaf(children, most_similar_leaf)
                child_nodes = self.get_child_nodes(children)
                all_leaves = self.check_if_all_sigleafs(child_nodes)

                if all_leaves:
                    child_similarity = most_similar_leaf.data.similarity(
                        other_child.node.data, ignore_abundance=self.ignore_abundance)

                    if new_leaf_similarity > child_similarity:
                        # New leaf is *more* similar than the existing child
                        # --> displace existing child

                        # Get this child's displaced position
                        displaced_position = other_child.pos

                        # Place the less similar child in the neighboring node
                        grandparent = self.parent(most_similar_parent.pos)
                        parent_sibling = [x for x in self.children(grandparent.pos)
                                          if x != most_similar_parent][0]
                        self.insert_new_internal_node_with_children(other_child.node,
                                                                    parent_sibling)
                        # Remove the old location
                        del self._leaves[displaced_position]
                        return displaced_position
                    else:
                        self.next_node = self._insert_next_position(self.next_node)
                else:
                    # One of the children is a Node rather than a SigLeaf --> replace
                    # the node with the SigLeaf
                    self.next_node = self._insert_next_position(self.next_node)


        else:
            self.next_node = self._insert_next_position(self.next_node)

        return self.next_node

    def _insert_next_position(self, next_node):
        # New leaf is *less* similar than the existing child
        # --> Create new adjacent parent as done previously
        min_leaf = min(self._leaves.keys())
        next_internal_node = None
        if next_node <= min_leaf:
            for i in range(min_leaf):
                if all((i not in self._nodes,
                        i not in self._leaves,
                        i not in self._missing_nodes)):
                    next_internal_node = i
                    break
        if next_internal_node is None:
            next_node = max(self._leaves.keys()) + 1
        else:
            next_node = next_internal_node

        return next_node

    def add_node(self, node):
        pos = self.new_node_pos(node)

        if pos == 0:  # empty tree; initialize w/node.
            new_internal_node = Node(self.factory, name="internal." + str(pos))
            self._nodes[0] = new_internal_node
            pos = self.new_node_pos(node)

        # Cases:
        # 1) parent is a Leaf (already covered)
        # 2) parent is a Node (with empty position available)
        #    - add Leaf, update parent
        # 3) parent is a Node (no position available)
        #    - this is covered by case 1
        # 4) parent is None
        #    this can happen with d != 2, in this case create the parent node
        parent = self.parent(pos)
        if isinstance(parent.node, Leaf):
            # Create a new internal node
            # node and parent are children of new internal node
            # Move children of grandparent internal node to new inserted internal node
            # below, and bump new node up to grandparent node
            grandparent = self.parent(parent.pos)
            grandparent_children = self.children(grandparent.pos)
            if grandparent_children[0] == parent:
                parent_sibling = grandparent_children[1]
            else:
                parent_sibling = grandparent_children[0]
            if isinstance(parent_sibling.node, Node):
                self.insert_new_internal_node_with_children(node, parent)
            else:
                self.relocate_children_to_new_internal_node(grandparent, node, parent,
                                                            parent_sibling)
        elif isinstance(parent.node, Node):
            self._leaves[pos] = node
            node.update(parent.node)
        elif parent.node is None:
            new_internal_node = Node(self.factory, name="internal." + str(parent.pos))
            self._nodes[parent.pos] = new_internal_node
            c1 = self.children(parent.pos)[0]
            self._leaves[c1.pos] = node
            node.update(new_internal_node)

        # update all parents!
        parent = self.parent(parent.pos)
        while parent:
            self._rebuild_node(parent.pos)
            node.update(self._nodes[parent.pos])
            parent = self.parent(parent.pos)

    def relocate_children_to_new_internal_node(self, grandparent, node, parent,
                                               parent_sibling):
        new_internal_node = Node(self.factory, name="internal." + str(parent.pos))
        self._nodes[parent.pos] = new_internal_node
        c1, c2 = self.children(parent.pos)
        # Update new internal node
        self._leaves[c1.pos] = parent.node
        self._leaves[c2.pos] = parent_sibling.node
        del self._leaves[parent.pos]
        # Swap position of previous parent's sibling
        self._leaves[parent_sibling.pos] = node
        # Set the parent of the new node as the grandparent
        node.update(grandparent.node)
        # Update data for the moved parent and parent's sibling
        for child in (parent.node, parent_sibling.node):
            child.update(new_internal_node)


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
        """May return a list of matches under a node --> doesn't return a single sig"""
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
