from __future__ import print_function
from __future__ import division

from io import BytesIO
import math
import sys

from .sbt import Leaf, Node, SBT, GraphFactory
from . import signature


def load_sbt_index(filename, print_version_warning=True):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning)


def create_sbt_index(bloom_filter_size=1e5, n_children=2):
    "Create an empty SBT index."
    factory = GraphFactory(1, bloom_filter_size, 4)
    if n_children != 2:
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

        buf = signature.save_signatures([self.data], compression=1)
        return self.storage.save(path, buf)

    def update(self, parent):
        mh = self.data.minhash
        parent.data.update(mh)
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

        # When recursively inserting new leaves, ignore these.
        # Set and reset after doing recursion
        self._recursion_ignore_leaves = set([])

    def find_most_similar_leaf(self, node_to_add):
        search_results = self.search(
            node_to_add.data,
            threshold=sys.float_info.epsilon,
            best_only=True,
            ignore_abundance=self.ignore_abundance,
            do_containment=self.do_containment,
            return_leaf=True,
            ignore_empty=True,
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

        if isinstance(node, SigLeaf):
            self.next_node = self.get_best_leaf_position(node)
        else:
            self.next_node = self._insert_next_position(self.next_node)

        return self.next_node

    def get_best_leaf_position(self, node_to_add):
        best_result = self.find_most_similar_leaf(node_to_add)
        if best_result is None:
            potential_next_node = self._insert_next_position(self.next_node)
            parent_potential_next_node = self.parent(potential_next_node)
            if parent_potential_next_node.pos in self._leaves:
                # No best position for this node. Push current tree down as-is, and
                # insert this new leaf in adjacent neighbors
                next_node = self._push_existing_tree_down()
            else:
                # There's an available empty leaf! Put this new one here
                next_node = potential_next_node
            return next_node

        new_leaf_similarity, most_similar_leaf, most_similar_pos = best_result

        # Get parent of the most similar node
        most_similar_parent = self.parent(most_similar_pos)

        # If the parent has one child: easy, insert the new child here
        children_of_most_similar_parent = self.children(most_similar_parent.pos)
        if children_of_most_similar_parent[1].node is None:
            # Use the default next node position
            next_node = self._insert_next_position(self.next_node)
        else:
            next_node = self.maybe_displace_child(
                new_leaf_similarity,
                children_of_most_similar_parent,
                most_similar_leaf)
        return next_node

    def maybe_displace_child(self, new_leaf_similarity, children, most_similar_leaf):
        # If parent has two children, check if the other child is more similar
        # to the most_similar_leaf --> then no displacement is necessary
        sibling_of_best_match = self.get_sibling_of_similar_leaf(children, most_similar_leaf)
        child_nodes = self.get_child_nodes(children)
        all_leaves = self.check_if_all_sigleafs(child_nodes)
        # Check if any currently created parental nodes have an empty
        # leaf child slot
        exists_free_leaf = self.check_exists_free_leaf()

        if all_leaves:
            child_similarity = most_similar_leaf.data.similarity(
                sibling_of_best_match.node.data, ignore_abundance=self.ignore_abundance)

            if new_leaf_similarity > child_similarity:
                next_node = self.displace_child_with_new_leaf(
                    exists_free_leaf,
                    sibling_of_best_match,
                    most_similar_leaf)
            else:
                next_node = self.insert_dissimilar_leaf(exists_free_leaf)
        else:
            # One of the children is a Node rather than a SigLeaf --> replace
            # the node with the SigLeaf
            next_node = self._insert_next_position(self.next_node)
        return next_node

    def insert_dissimilar_leaf(self, exists_free_leaf):
        # This node isn't similar to anything
        if exists_free_leaf:
            # If there is a next available spot, take it
            next_node = self._insert_next_position(self.next_node)
        else:
            # Otherwise, keep the current structure as-is and insert a new node,
            # bottom-up
            next_node = self._push_existing_tree_down()
        return next_node

    def displace_child_with_new_leaf(self, exists_free_leaf, node_to_displace,
                                     most_similar_leaf):
        # New leaf is *more* similar than the existing child
        # --> displace existing child
        # Get this child's displaced position
        displaced_position = node_to_displace.pos
        if exists_free_leaf:
            # Mask both the node to displace and the most similar leaf away
            node_to_displace.node._ignore_in_search = True
            most_similar_leaf._ignore_in_search = True
            node_to_displace_new_pos = self.get_best_leaf_position(
                node_to_displace.node)
            # self._leaves[node_to_displace_new_pos] = node_to_displace.node

            # Unset ignoring after recursion is done
            for leaf in self.leaves():
                try:
                    leaf._ignore_in_search = False
                except AttributeError:
                    # Empty leaf --> NodePos object with no data
                    pass

            # displaced_node_new_position = self._insert_next_position(self.next_node)
            displaced_node_new_parent = self.parent(node_to_displace_new_pos)
        # elif has_grandparent and exists_free_leaf:
        #     displaced_node_new_parent = [x for x in self.children(grandparent.pos)
        #                              if x != most_similar_parent][0]
        else:
            displaced_node_new_position = self._push_existing_tree_down()
            # Update displaced position to where the child was in the new
            # position when the tree was pushed down
            for pos, leaf in self._leaves.items():
                if leaf == node_to_displace.node:
                    displaced_position = pos
                    break
            displaced_node_new_parent = self.parent(displaced_node_new_position)
        self.insert_new_internal_node_with_children(node_to_displace.node,
                                                    displaced_node_new_parent.pos)
        del self._leaves[displaced_position]
        return displaced_position

    def check_exists_free_leaf(self):
        potential_next_position = self._insert_next_position(self.next_node)
        potential_next_position_parent = self.parent(potential_next_position)
        if potential_next_position_parent.pos in self._nodes:
            exists_free_leaf = True
        else:
            exists_free_leaf = False
        return exists_free_leaf

    def _push_existing_tree_down(self):
        current_tree_depth = int(math.floor(math.log2(max(self._leaves))))
        add_to_leaves = 2 ** current_tree_depth
        add_to_nodes = 2 ** (current_tree_depth - 1)

        new_leaves = {}
        for i, leaf in self._leaves.items():
            new_leaves[i + add_to_leaves] = leaf
        self._leaves = new_leaves

        new_nodes = {}
        for i, (pos, node) in enumerate(self._nodes.items()):
            new_position = pos + add_to_nodes
            node.name = "internal." + str(new_position)
            new_nodes[new_position] = node

        self._nodes = new_nodes
        # Rebuild new nodes in case they don't exist
        for i in range(min(self._leaves.keys())):
            self._rebuild_internal_nodes_bottom_up(i)
        next_node = max(new_leaves) + 1
        return next_node

    def _rebuild_internal_nodes_bottom_up(self, pos=None):
        """Rebuilds internal nodes (if it is not present), recursively up the tree

        Parameters
        ----------
        pos: int
            node to be rebuild. Any internal node *above* it will be rebuild too.
            If you want to rebuild all missing internal nodes you can use pos=None
            (the default).
        """
        if pos is None:
            pos = min(self._leaves.keys()) - 1

        node = self._nodes.get(pos, None)
        if node is not None:
            # this node was already build, skip

            # But make sure its parents are built
            self._rebuild_from_position(pos)
            return

        node = Node(self.factory, name="internal.{}".format(pos))
        self._nodes[pos] = node
        self.update_children_nodes_below_position(node, pos)

        self._rebuild_from_position(pos)

    def update_children_nodes_below_position(self, node, pos):
        for c in self.children(pos):
            if c.node is not None:
                # Node may be empty because of bottom up building
                c.node.update(node)
            if c.pos in self._missing_nodes:
                node = Node(self.factory, name="internal.{}".format(c.pos))
                self._nodes[c.pos] = node
                self.update_children_nodes_below_position(node, c.pos)

    def _rebuild_from_position(self, pos):
        """Check parent of current nodes"""
        parent = self.parent(pos)
        # Rebuild all the way to the top!
        if parent is not None:
            self._rebuild_internal_nodes_bottom_up(parent.pos)

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
                self.insert_new_internal_node_with_children(node, parent.pos,
                                                            parent.node)
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

        # Rebuild all nodes, starting from the node previous to existing ones
        self._rebuild_internal_nodes_bottom_up()

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

    def insert_new_internal_node_with_children(self, node, parent_pos,
                                               parent_node=None):
        """Create a new internal node at parent.pos with child, 'node'"""
        # Create a new internal node
        # node and parent are children of new internal node
        if parent_pos not in self._nodes:
            n = Node(self.factory, name="internal." + str(parent_pos))
            self._nodes[parent_pos] = n
        else:
            n = self._nodes[parent_pos]
        c1, c2 = self.children(parent_pos)[:2]
        if parent_node is None:
            # Don't split this node into a parent - just keep it empty
            if c1.node is None:
                # First position is empty
                self._leaves[c1.pos] = node
            else:
                # First position is taken
                # The other position is free
                self._leaves[c2.pos] = node
            node.update(n)
        if parent_node is not None:
            self._leaves[c1.pos] = parent_node
            self._leaves[c2.pos] = node
            del self._leaves[parent_pos]
            for child in (parent_node, node):
                child.update(n)

### Search functionality.

def _max_jaccard_underneath_internal_node(node, mh, ignore_empty=False):
    """\
    calculate the maximum possibility similarity score below
    this node, based on the number of matches in 'hashes' at this node,
    divided by the smallest minhash size below this node.

    This should yield be an upper bound on the Jaccard similarity
    for any signature below this point.

    ignore_empty : bool
        If true, then when an internal node is empty and has no children, don't throw
        an error. Needed for localized SBTs, which many have empty internal nodes in the
        process of building
    """
    if len(mh) == 0:
        return 0.0

    # count the maximum number of hash matches beneath this node
    matches = node.data.matches(mh)

    # get the size of the smallest collection of hashes below this point
    default = 0 if ignore_empty else -1
    min_n_below = node.metadata.get('min_n_below', default)

    if min_n_below == -1:
        raise Exception('cannot do similarity search on this SBT; need to rebuild.')
    if min_n_below == 0:
        return min_n_below

    # max of numerator divided by min of denominator => max Jaccard
    max_score = float(matches) / min_n_below

    return max_score


def search_minhashes(node, sig, threshold, results=None, ignore_empty=False):
    """\
    Default tree search function, searching for best Jaccard similarity.
    """
    sig_mh = sig.minhash
    score = 0

    if isinstance(node, SigLeaf):
        score = node.data.minhash.similarity(sig_mh)
    else:  # Node minhash comparison
        score = _max_jaccard_underneath_internal_node(node, sig_mh,
                                                      ignore_empty=ignore_empty)

    if results is not None:
        results[node.name] = score

    if score >= threshold:
        return 1

    return 0


class SearchMinHashesFindBest(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None, ignore_empty=False):
        """May return a list of matches under a node --> doesn't return a single sig"""
        sig_mh = sig.minhash
        score = 0

        if isinstance(node, SigLeaf):
            score = node.data.minhash.similarity(sig_mh)
        else:  # internal object, not leaf.
            score = _max_jaccard_underneath_internal_node(node, sig_mh,
                                                          ignore_empty=ignore_empty)

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


def search_minhashes_containment(node, sig, threshold, results=None, downsample=True,
                                 ignore_empty=None):
    """Find minhashes contained within

    Parameters
    -----------
    ignore_empty: bool
        Needed for compatibility with _max_jaccard_underneath_internal_node
    """
    mh = sig.minhash

    if isinstance(node, SigLeaf):
        matches = node.data.minhash.count_common(mh, downsample)
    else:  # Node or Leaf, Nodegraph by minhash comparison
        matches = node.data.matches(mh)

    if results is not None:
        results[node.name] = float(matches) / len(mh)

    if len(mh) and float(matches) / len(mh) >= threshold:
        return 1
    return 0


class GatherMinHashes(object):
    def __init__(self):
        self.best_match = 0

    def search(self, node, query, threshold, results=None):
        mh = query.minhash
        if not len(mh):
            return 0

        if isinstance(node, SigLeaf):
            matches = mh.count_common(node.data.minhash, True)
        else:  # Nodegraph by minhash comparison
            matches = node.data.matches(mh)

        if not matches:
            return 0

        score = float(matches) / len(mh)

        if score < threshold:
            return 0

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
