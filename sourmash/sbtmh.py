from __future__ import print_function
from __future__ import division

from .sbt import Leaf
from . import MinHash, signature


class SigLeaf(Leaf):
    def __str__(self):
        return '**Leaf:{name} -> {metadata}'.format(
                name=self.name, metadata=self.metadata)

    def save(self, filename):
        # this is here only for triggering the property load
        # before we reopen the file (and overwrite the previous
        # content...)
        self.data

        with open(filename, 'w') as fp:
            signature.save_signatures([self.data], fp)

    def update(self, parent):
        for v in self.data.estimator.get_mins():
            parent.data.count(v)

    @property
    def data(self):
        if self._data is None:
            it = signature.load_signatures(self._filename)
            self._data, = list(it)              # should only be one signature
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data


def search_minhashes(node, sig, threshold, results=None):
    mins = sig.estimator.get_mins()

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
        mins = sig.estimator.get_mins()

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


class SearchMinHashesFindBestIgnoreMaxHash(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        mins = sig.estimator.get_mins()

        if isinstance(node, SigLeaf):
            old_est = node.data.estimator
            E = MinHash(ksize=old_est.ksize, n=old_est.num)
            for m in old_est.get_mins():
                E.add_hash(m)

            matches = E.count_common(sig.estimator)
        else:  # Node or Leaf, Nodegraph by minhash comparison
            matches = sum(1 for value in mins if node.data.get(value))

        score = 0
        if not len(mins):
            return 0

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
