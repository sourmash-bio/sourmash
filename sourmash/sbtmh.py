from __future__ import print_function
from __future__ import division

from io import BytesIO, TextIOWrapper

from .sbt import Leaf, SBT, GraphFactory
from . import signature


def load_sbt_index(filename):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load)


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
    for leaf in tree.find(search_minhashes, query, threshold):
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
        for v in self.data.minhash.get_mins():
            parent.data.count(v)
        max_n_below = parent.metadata.get('max_n_below', 0)
        max_n_below = max(len(self.data.minhash.get_mins()),
                          max_n_below)
        parent.metadata['max_n_below'] = max_n_below

    @property
    def data(self):
        if self._data is None:
            buf = BytesIO(self.storage.load(self._path))
            with TextIOWrapper(buf) as data:
                self._data = signature.load_one_signature(data)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data


def search_minhashes(node, sig, threshold, results=None, downsample=True):
    mins = sig.minhash.get_mins()
    score = 0

    if isinstance(node, SigLeaf):
        try:
            score = node.data.minhash.similarity(sig.minhash)
        except Exception as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = sig.minhash.downsample_max_hash(node.data.minhash)
                yy = node.data.minhash.downsample_max_hash(sig.minhash)

                score = yy.similarity(xx)
            else:
                raise

    else:  # Node or Leaf, Nodegraph by minhash comparison
        if len(mins):
            matches = sum(1 for value in mins if node.data.get(value))
            max_mins = node.metadata.get('max_n_below', -1)
            if max_mins == -1:
                raise Exception('cannot do similarity search on this SBT; need to rebuild.')
            score = float(matches) / max_mins

    if results is not None:
        results[node.name] = score

    if score >= threshold:
        return 1

    return 0


class SearchMinHashesFindBest(object):
    def __init__(self, downsample=True):
        self.best_match = 0.
        self.downsample = downsample

    def search(self, node, sig, threshold, results=None):
        mins = sig.minhash.get_mins()
        score = 0

        if isinstance(node, SigLeaf):
            try:
                score = node.data.minhash.similarity(sig.minhash)
            except Exception as e:
                if 'mismatch in max_hash' in str(e) and self.downsample:
                    xx = sig.minhash.downsample_max_hash(node.data.minhash)
                    yy = node.data.minhash.downsample_max_hash(sig.minhash)

                    score = yy.similarity(xx)
                else:
                    raise
        else:  # internal object, not leaf.
            if len(mins):
                matches = sum(1 for value in mins if node.data.get(value))
                max_mins = node.metadata.get('max_n_below', -1)
                if max_mins == -1:
                    raise Exception('cannot do similarity search on this SBT; need to rebuild.')
                score = float(matches) / max_mins

        if results is not None:
            results[node.name] = score

        if score >= threshold:
            # have we done better than this? if yes, truncate.
            if score > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = score
                return 1

        return 0


def search_minhashes_containment(node, sig, threshold,
                                 results=None, downsample=True):
    mins = sig.minhash.get_mins()

    if isinstance(node, SigLeaf):
        try:
            matches = node.data.minhash.count_common(sig.minhash)
        except Exception as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = sig.minhash.downsample_max_hash(node.data.minhash)
                yy = node.data.minhash.downsample_max_hash(sig.minhash)

                matches = yy.count_common(xx)
            else:
                raise

    else:  # Node or Leaf, Nodegraph by minhash comparison
        matches = sum(1 for value in mins if node.data.get(value))

    if results is not None:
        results[node.name] = float(matches) / len(mins)

    if len(mins) and float(matches) / len(mins) >= threshold:
        return 1
    return 0


class SearchMinHashesFindBestIgnoreMaxHash(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        mins = sig.minhash.get_mins()

        if isinstance(node, SigLeaf):
            max_scaled = max(node.data.minhash.scaled, sig.minhash.scaled)

            mh1 = node.data.minhash.downsample_scaled(max_scaled)
            mh2 = sig.minhash.downsample_scaled(max_scaled)
            matches = mh1.count_common(mh2)
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
