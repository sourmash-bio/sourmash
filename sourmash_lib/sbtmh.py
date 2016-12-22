from __future__ import print_function
from __future__ import division
from glob import glob
import os

from sourmash_lib import signature
from sourmash_lib import sourmash_tst_utils as utils
from .sbt import SBT, GraphFactory, Leaf


class SigLeaf(Leaf):
    def __str__(self):
        return '**Leaf:{name} -> {metadata}'.format(
                name=self.name, metadata=self.metadata)

    def save(self, filename):
        from sourmash_lib import signature
        with open(filename, 'w') as fp:
            signature.save_signatures([self.data], fp)

    def update(self, parent):
        for v in self.data.estimator.mh.get_mins():
            parent.data.count(v)

    @staticmethod
    def load(info, dirname):
        from sourmash_lib import signature

        filename = os.path.join(dirname, info['filename'])
        with open(filename, 'r') as fp:
            data = next(signature.load_signatures(fp))
        return SigLeaf(info['metadata'], data, name=info['name'])


def search_minhashes(node, sig, threshold, results=None):
    mins = sig.estimator.mh.get_mins()

    if isinstance(node, SigLeaf):
        matches = node.data.estimator.count_common(sig.estimator)
    else:  # Node or Leaf, Nodegraph by minhash comparison
        matches = sum(1 for value in mins if node.data.get(value))

    if results is not None:
        results[node.name] = float(matches) / len(mins)

    if len(mins) and float(matches) / len(mins) >= threshold:
        return 1
    return 0


class SearchMinHashesFindBest(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        mins = sig.estimator.mh.get_mins()

        if isinstance(node, SigLeaf):
            matches = node.data.estimator.count_common(sig.estimator)
        else:  # Node or Leaf, Nodegraph by minhash comparison
            matches = sum(1 for value in mins if node.data.get(value))

        score = 0
        if len(mins):
            score = float(matches) / len(mins)

        if results is not None:
            results[node.name] = score

        if score >= threshold:
            # have we done better than this? if yes, truncate.
            if float(matches) / len(mins) > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = float(matches) / len(mins)
                return 1

        return 0


def test_tree_save_load():
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory)

    for f in glob("demo/*.sig"):
        with open(f, 'r') as data:
            sig = next(signature.load_signatures(data))
        leaf = SigLeaf(os.path.basename(f), sig)
        try:
            tree.add_node(leaf)
        except TypeError:
            raise Exception("check your version of khmer - you need khmer>2.1")
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    old_result = [str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)]
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'))
        tree = SBT.load(os.path.join(location, 'demo'),
                        leaf_loader=SigLeaf.load)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = [str(s) for s in tree.find(search_minhashes,
                                                to_search.data, 0.1)]
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_binary_nary_tree():
    factory = GraphFactory(31, 1e5, 4)
    trees = {}
    trees[2] = SBT(factory)
    trees[5] = SBT(factory, d=5)
    trees[10] = SBT(factory, d=10)

    for f in glob("demo/*.sig"):
        with open(f, 'r') as data:
            sig = next(signature.load_signatures(data))
        leaf = SigLeaf(os.path.basename(f), sig)
        for tree in trees.values():
            tree.add_node(leaf)
        to_search = leaf

    results = {}
    print('*' * 60)
    print("{}:".format(to_search.metadata))
    for d, tree in trees.items():
        results[d] = [str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)]
    print(*results[2], sep='\n')

    assert set(results[2]) == set(results[5])
    assert set(results[5]) == set(results[10])

if __name__ == "__main__":
    test_binary_nary_tree()
