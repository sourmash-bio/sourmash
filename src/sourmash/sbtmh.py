from io import BytesIO
import sys

from .sbt import Leaf, SBT, GraphFactory
from . import signature


def load_sbt_index(filename, *, print_version_warning=True, cache_size=None):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning,
                    cache_size=cache_size)


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
    for (score, match, _) in tree.search(query, threshold=threshold,
                                         unload_data=True):
        yield match, score


class SigLeaf(Leaf):
    def __str__(self):
        return '**Leaf:{name} -> {metadata}'.format(
                name=self.name, metadata=self.metadata)

    def make_manifest_row(self, loc):
        from .index import CollectionManifest
        row = CollectionManifest.make_manifest_row(self.data,
                                                   loc,
                                                   include_signature=0)
        return row

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
